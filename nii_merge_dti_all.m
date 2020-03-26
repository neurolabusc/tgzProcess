function nii_merge_dti_all(pth)
%search folder and subfolders, combine all matching 'DTIrev_'/'DTI_' files in same folder
% pth : parent folder to search
%Examples
% nii_merge_dti_all('/Users/rorden/tst/mat2bids/tgz');

%pth = '/Users/rorden/tst/org/m';
if ~exist('pth','var')
    pth = uigetdir(pwd);
end
mergeDtiSub(pth);
%end nii_merge_dti_all()

function mergeDtiSub(pth)
fnms=subFileSub(pth, {'DTI_*.nii', 'DTI_*.nii.gz'});
combineSub(pth, fnms);
fnms=subFileSub(pth, {'DTIrev_*.nii', 'DTIrev_*.nii.gz'});
combineSub(pth, fnms);
%recursively check sub folders
pths=subFileSub(pth);
for i = numel(pths) : -1 : 1
    mergeDtiSub(fullfile(pth, pths{i}));
end
%end mergeDtiSub()

function combineSub(pth, fnms)
if numel(fnms) < 2, return; end
%next optional: exclude images we have already processed
f = [];
for i = 1 : numel(fnms)
    [~, nam] = fsl_filepartsSub(fnms{i});
    if endsWith(nam,'_Cat'), continue; end
    f = [f; fnms(i)]; %#ok<AGROW>
end
fnms = f;
if numel(fnms) < 2, return; end
%load all images to memory
niis = [];
for i = 1 : numel(fnms)
    fnm = fullfile(pth, fnms{i});
    niis = [niis, nii_tool('load', fnm)]; %#ok<AGROW>
end
%assume images with same X,Y,Z dims should be stacked
isMatched = zeros(numel(fnms),1);
for i = 1 : (numel(fnms)-1)
    for j = i : numel(fnms)
        if (isMatched(j)), continue; end
        if isequal(niis(i).hdr.dim(2:4), niis(j).hdr.dim(2:4))
            isMatched(j) = i;
        end
    end
end
for i = 1: max(isMatched)
    idx = find(isMatched == i);
    if numel(idx) < 2, continue; end
    nii.hdr = niis(idx(1)).hdr;
    nii.img = [];
    bvalCat = [];
    bvecCat = [];
    teCat = [];
    [pth, nam, ext] = fsl_filepartsSub(niis(idx(1)).hdr.file_name);
    for j = 1: numel(idx)
        nii.img = cat(4,nii.img, niis(idx(j)).img);
        %getBvalBvecSub
        te = getEchoTimeSub(niis(idx(j)).hdr.file_name);
        teCat = [teCat; te]; %#ok<AGROW>
        [bval, bvec] = getBvalBvecSub (niis(idx(j)).hdr.file_name, size(niis(idx(j)).img, 4));
        bvecCat = [bvecCat; bvec]; %#ok<AGROW>
        bvalCat = [bvalCat; bval]; %#ok<AGROW>
        jsonNam = fullfile(pth, [nam, '.json']);
    end
    nii.hdr.dim(5) = size(nii.img, 4);
    if max(bvalCat(:)) < 200
        fprintf('%d volumes all with B=0 %s\n', nii.hdr.dim(5), nam);
        fid = fopen('Bzero.txt', 'a');
        fprintf(fid,'%s\n', nam);
        fclose(fid);
        continue;
    end
    fprintf('%d volumes in merged %s\n', nii.hdr.dim(5), nam);
    %exit here unless we want to apply these ...
    % continue;
    for j = 1: numel(idx)
        if exist(jsonNam, 'file')
           copyfile(jsonNam, fullfile(pth, [nam, '_Cat.json'])); 
        end
        renameSub(niis(idx(j)).hdr.file_name,'_');
    end
    nii_tool('save', nii, fullfile(pth, [nam, '_Cat',ext]));
    dlmwrite(fullfile(pth,[ nam, '_Cat.bval']), bvalCat','delimiter','\t');
    dlmwrite(fullfile(pth,[nam, '_Cat.bvec']), bvecCat','delimiter','\t');
    teCat = sort(unique(teCat));
    if numel(teCat) > 1 
        fprintf('warning %d echo times: correct if you calculate Kurtosis %s\n', numel(teCat), nam);
        fid = fopen(fullfile(pth, [nam, '_Cat.te']),'w');
        fprintf(fid,'{\n\t"EchoTime": [');
        fprintf(fid,' %0.4f',teCat(1));
        for j = 2: numel(teCat)
            fprintf(fid,', %0.4f',teCat(j));
        end
        fprintf(fid,' ]\n}\n');
        fclose(fid);       
    end %echo times vary
end %for matched
%end combineSub()

function renameSub(fnm,prefix)
%add prefix to image, bvec, bval names
[pth,nam] = fsl_filepartsSub(fnm);
outnam = [prefix, nam];
x = { '.json','.nii','.nii.gz','.bvec','.bval'};
fprintf(' + %s\n', nam);
for i = 1: numel(x)
    inx = fullfile(pth, [nam, x{i}]);
    if ~exist(inx, 'file'); continue; end
    outx = fullfile(pth, [outnam, x{i}]);
    movefile(inx, outx);
end

function [bval, bvec] = getBvalBvecSub (imgName, nVol)
%read .bval file and return indices for B0 volumes
[pth, nam] = fsl_filepartsSub(imgName);
nameVal = fullfile(pth,[nam '.bval']); %name for b-values
nameVec = fullfile(pth,[nam  '.bvec']); %name for b-vectors
%if (exist(imgName, 'file') == 0) , fprintf('Unable to find required image %s\n',imgName); return; end;
if ( (exist(nameVal, 'file') == 0) || (exist(nameVec, 'file') == 0) )
    fprintf('Assumg B=0: Unable to find required DTI files %s and %s\n',nameVec,nameVal);
    bval = zeros(nVol,1);
    bvec = zeros(nVol,3);
    return;
end
%bval = importdata(nameVal); %<- does not work with Matlab 2014b on Linux
%read b-values
fileID = fopen(nameVal);
bval = cell2mat( textscan(fileID,'%d'));
fclose(fileID);
%read b-vectors
fileID = fopen(nameVec);
bvec = cell2mat( textscan(fileID,'%f'));
fclose(fileID);
if mod(numel(bvec),3) ~= 0
    error('Error: number of bvecs must be divisible by three. Found %d bvecs in %s', numel(bvec), nameVec);
end
bvec = reshape(bvec,numel(bvec)/3,3);
%end getB0vols()

function [pth, nam, ext] = fsl_filepartsSub(fileName)
% a.nii.gz has the extension ".nii.gz" not ".nii"
[pth, nam, ext] = fileparts(fileName);
if (length(ext)==3)  && min((ext=='.gz')==1)
	[pth, nam, ext2] = fileparts( fullfile(pth, nam)); %remove .nii for .nii.gz
    ext = [ext2 ext];
end
%end fsl_filepartsSub()

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
         nameFiles= [nameFiles; subFileSub(pathFolder, ext{i})]; %#ok<AGROW>
     end
     return;
else
    if contains(ext,'*')
        d = dir(fullfile(pathFolder, ext));
    else
        d = dir(fullfile(pathFolder,['*' ext]));
    end
    isub = ~[d(:).isdir];
    nameFiles = {d(isub).name}';
end
nameFiles=nameFiles(~startsWith(nameFiles(:),'_'));
nameFiles=nameFiles(~startsWith(nameFiles(:),'.'));
nameFiles = sort(nameFiles);
%end subFileSub()

function te = getEchoTimeSub(imgName)
te = 0;
[pth, nam] = fsl_filepartsSub(imgName);
nameJson = fullfile(pth,[nam '.json']); %name for b-values
if ~exist(nameJson, 'file'), return; end
json = LoadJsonSub(nameJson);
if isfield(json,'EchoTime')
    te = json.('EchoTime');
end
%end getEchoTimeSub()

function [json, fnm] = LoadJsonSub(fnm)
p = fileparts(which(mfilename));
if ~exist('fnm','var'), fnm = fullfile(p, 'config.json'); end
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