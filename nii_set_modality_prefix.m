function nii_set_modality_prefix(pth, isRename)
% use .json to set modality, e.g. T1 scans named 'T1_filename'
% pth : parent folder to search for images
% isRename : to change filenames (otherwise just reports filenames)
% Examples
% dcm2modality('/Users/rorden/tst/mat2bids/tgz', true);
% 
%pth = '/Volumes/TbStick/nifti/20180628_114055'
%pth = '/Volumes/TbStick/nifti/20110317_135407'
if ~exist('pth','var')
    pth = uigetdir(pwd);
end
if ~exist('isRename','var')
    isRename = true;
end
setModalitySub(pth, isRename);
%end dcm2modality(pth);

function setModalitySub(pth, isRename)
nams = subFileSub(pth, '.json');
fixBrokenModalities = false;
for i = 1:numel(nams)
    njson = nams{i}
    if njson(1) == '.', continue; end
    modality = [];
    [~, n] = fileparts(njson); %remove extension
    n = fullfile(pth,n);
    njson = fullfile(pth,njson);
    modality = modalitySub(njson);
    oldmodality = strsplit(nams{i},'_');
    
    if ~isempty(oldmodality) && strcmp(oldmodality{1},modality), continue; end; %already set
    if fixBrokenModalities
        oldmodality = oldmodality{1};
        if max(strcmpi({'FA','TRACE','ADC','swi','angio','fmapSErev','DTIrev','fmapSE','DTI','fmap', 'fMRI','Rest','ASL','T2','T1','T1w','func_rest','func_task','fmap','SWI','PCASL','dwi_AP','dwi_rev','FLAIR','anat_T1w_e1'},oldmodality)) < 1, 
            %error('Unknown modality %s', oldmodality);
            continue; 
        end
        if strcmp(oldmodality,modality), continue; end
        fprintf('%s -> %s: %s\n', oldmodality, modality, njson);
        renameSub(njson, modality, oldmodality);
        continue;
    end
    fprintf('>> %d/%d %s = %s\n', i, numel(nams), modality, nams{i});
    %if isempty(modality), error('Unknown modality %s\n', njson); end
    if ~isRename, continue; end
    renameSub(njson, modality);
    modality
end
%next: recursively search sub folders
nams = subFileSub(pth);
for i = 1:numel(nams)
    setModalitySub(fullfile(pth,nams{i}), isRename)
end
%end dcm2modality()

function renameSub(fnm, modality, oldmodality)
[pth,nam] = fileparts(fnm);
modality = [modality,'_'];
if startsWith(nam,modality), return; end
%fprintf('>> %s %s\n', nam, modality);
outnam = nam;
if exist('oldmodality','var')
    outnam = nam((numel(oldmodality)+2):end);
end
x = { '.json','.nii','.nii.gz','.bvec','.bval'};
for i = 1: numel(x)
    nx = [nam, x{i}];
    inx = fullfile(pth, nx);
    if ~exist(inx, 'file'); continue; end
    outx = fullfile(pth, [modality, outnam, x{i}]);
    if strcmpi(inx, outx)
        tmp = fullfile(pth,'tempMat');
        movefile(inx, tmp);
        movefile(tmp, outx);
        continue;
    end
    %fprintf('%s -> %s\n', inx, outx);
    movefile(inx, outx);    
end

function modality = modalitySub(njson)
modality = [];
if ~exist(njson, 'file'), return; end
[json] = LoadJsonSub(njson);
seqNam = [];
if isfield(json,'SequenceName')
    seqNam = json.('SequenceName');
end
isPhaseRev = false;
if isfield(json,'PhaseEncodingDirection')
    isPhaseRev = json.('PhaseEncodingDirection') == 'j';
end
seqVar = [];
if isfield(json,'SequenceVariant')
    seqVar = json.('SequenceVariant');
end
TE = nan;
if isfield(json,'EchoTime')
    TE = json.('EchoTime');
end
protNam = [];
if isfield(json,'ProtocolName')
    protNam = json.('ProtocolName');
end
imageType = [];
if isfield(json,'ImageType')
    imageType = json.('ImageType');
    
end
%check values
if contains(seqNam, '_fm2d2r')
    modality = 'fmap';
    if isPhaseRev, modality = [modality, 'rev']; end
    return;
end
if contains(seqNam, '_tfl3d1') || isequal(seqVar, 'SP_MP')
    modality = 'T1';
    return;
end
if contains(seqNam, 'ep_b')
    if ~isempty(find(strcmp(imageType,'FA')))
        modality = 'FA';
    elseif ~isempty(find(strcmp(imageType,'ADC')))
        modality = 'ADC';   
    elseif max(contains(imageType,'TRACE')) %TRACE or TRACEW
        modality = 'TRACE';
    else
        modality = 'DTI';
    end
    %njson
    %imageType
    %modality
    if isPhaseRev, modality = [modality, 'rev']; end
    return;
end
if contains(seqVar, 'TOF_SP')
    modality = 'angio';
    return;
end
if contains(seqNam, '_fl3d1r')  %this check MUST be done after angio detected!
    modality = 'swi';
    return;
end

if contains(protNam, 'anat_SWI')  %this check MUST be done after angio detected!
    modality = 'swi';
    return;
end

if contains(seqNam, 'epse2d') %this check MUST be done after DWI detected!
    modality = 'fmapSE';
    if isPhaseRev, modality = [modality, 'rev']; end
    return;
end
if contains(seqVar, 'SK_SP_MP') && contains(seqNam, 'ir')
    modality = 'FLAIR';
    return;
end
if isequal(seqVar, 'SP') && contains(seqNam, 'fl3d1r')
    modality = 'FLAIR';
    return;
end
if isequal(seqVar, 'SP_MP') && contains(seqNam, '_tfl3d1')
    modality = 'FLAIR';
    return;
end

if isequal(seqVar, 'SP') && contains(seqNam, 'fm2d2r')
    modality = 'fmap';
    return;
end
if  (TE >= 0.03) && (TE < 0.04) && (isequal(seqVar, 'SK') || isequal(seqVar, 'SK_SS'))
    if isempty(regexpi(protNam, 'REST'))
        modality = 'fMRI';
    else
        modality = 'Rest';
    end
    if isPhaseRev, modality = [modality, 'rev']; end
    return;
end
if (isequal(seqVar, 'SK') || isequal(seqVar, 'SK_OSP'))  && (TE < 0.03)
    modality = 'ASL';
    return;
end
if isequal(seqVar, 'SK') && contains(seqNam, 'epfid2d1')  && (TE > 0.04)
    modality = 'ASL';
    return;
end
if contains(seqVar, 'SK') && (TE > 0.078)
    modality = 'T2';
    return;
end
if contains(seqVar, 'SK_SP_MP') && contains(seqNam ,'tfl_me3d5_16ns') && ~isfield(json,'EchoNumber')
    modality = 'T1';
    return;
end



fprintf('Unable to classify Var="%s" Nam="%s" TE=%g\n', seqVar, seqNam, TE);
%end modalitySub()

function [json, fnm] = LoadJsonSub(fnm)
p = fileparts(which(mfilename));
if ~exist('fnm','var'), fnm = fullfile(p, 'config.json'); end;
if ~exist(fnm,'file')
    fprintf('Creating empty JSON: Unable to find %s\n', fnm);
    json = [];
    return;
end
fid = fopen(fnm);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
if exist('jsondecode', 'builtin')
    json = jsondecode(str);
else
    js = parse_json(str); %https://www.mathworks.com/matlabcentral/fileexchange/20565-json-parser
    json = js{1};
end
%end LoadJsonSub()

function nameFiles=subFileSub(pathFolder, ext)
%return all folders excep those that begin with '.' or '_'
%nameFiles=subFileSub(pth, '.nii') %names of .nii
%nameFiles=subFileSub(pth, {'.nii', '.nii.gz'}) %filenames of either .nii or .gz
%nameFiles=subFileSub(pth) %return directories
if ~exist('ext','var') || isempty(ext) 
    d = dir(pathFolder);
    isub = [d(:).isdir];
    nameFiles = {d(isub).name}';
elseif iscell(ext) && numel(ext) > 1
    nameFiles = []; 
    for i = 1: numel(ext)
         nameFiles= [nameFiles; subFileSub(pathFolder, ext{i})];
     end
     return;
else
    d = dir(fullfile(pathFolder,['*' ext]));
    isub = ~[d(:).isdir];
    nameFiles = {d(isub).name}';    
end
nameFiles=nameFiles(~startsWith(nameFiles(:),'_'));
nameFiles=nameFiles(~startsWith(nameFiles(:),'.'));
%end subFileSub()