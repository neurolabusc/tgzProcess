function nii_deface_all(pth)

%Deface all 'T1_*.nii','T2_*.nii','T1_*.nii.gz','T2_*.nii.gz' in path and its subfolders
% pth : parent directory that contains T1/T2 scans
%Example
% nii_deface_all('/Users/rorden/tst/mat2bids/Master_DB')

%pth = '/Volumes/TbStick/nifti'
%pth = '/Users/rorden/tst/org/todo'

if ~exist('pth','var')
    pth = uigetdir(pwd);
end
defaceAnatSub(pth);
%end nii_deface_all()

function defaceAnatSub(pth)
niis=subFileSub(pth, {'.nii', '.nii.gz'});
anats = [];
for i = 1 : numel(niis)
    fnm = niis{i};
    if ~startsWith(fnm,'T1_') && ~startsWith(fnm,'FLAIR_') && ~startsWith(fnm,'T2_'), continue; end;
    fnm = fullfile(pth,fnm);
    %fprintf('%s\n', fnm);
    anats = [anats {fnm}];
    %T1 is stored in anats{2} for ABC study
    
end

%ROGER commented out these lines temporarily dont need epilog for all participants
%nm = pth(size(pth,2)-3:size(pth,2));
%copyfile(anats{2}, ['../T1s/' nm '_T1.nii.gz'])
    
if ~isempty(anats)
   nii_deface_crop(anats);
end
pths=subFileSub(pth);
%for i = 1 : numel(pths)
for i =  numel(pths) : -1 : 1
    fprintf('%d/%d\n', i, numel(pths));
    defaceAnatSub(fullfile(pth, pths{i}));
end
%end defaceAnatSub()

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
nameFiles = sort(nameFiles);
%end subFileSub()