function nii_match_z (basepth)
%this script was designed to ensure session dates (201601202) match participants (M2002)
%  since data was acquired anonymized across different studies, we need to ensure 
%  the correct scan is associated with the correct image
%Generate one representative bitmap for each scanning session
%
%Images are basepth/Participant/Session/*.nii
%For example:
% M2002/
%    201601202/
%       T1.nii
%       fMRI.nii
% M2012/
%    201801202/
%       T1.nii
%       fMRI.nii
%    201901202/
%       DTI.nii
%       fMRI.nii
%
basepth = '/Volumes/Chris5TB/Universe/master';
if ~exist('basepth','var')
    basepth = pwd;
end

tmpdir = fullfile(pwd, 'tmp');
if ~exist(tmpdir, 'dir')
   mkdir(tmpdir) 
end
subjs = dir(fullfile(basepth, 'M*'));
fnms={subjs.name};
[~,idx]=sort(fnms);
subjs=subjs(idx);
for s = 1: numel(subjs)
    if ~subjs(s).isdir, continue; end
    if ~isempty(strfind(subjs(s).name,'_')), continue; end
    subjpth = fullfile(basepth, subjs(s).name);
    %fprintf('%d/%d %s\n', s, numel(subjs), subjs(s).name);
    %simages = [];
    visits = dir(fullfile(subjpth, '*'));
    for v = 1: numel(visits)
        if ~visits(v).isdir, continue; end
        if visits(v).name(1) == '.', continue; end
        %fprintf(' %s\n', visits(v).name);
        visitpth = fullfile(subjpth, visits(v).name);
        img = representative_img(visitpth);
    end
end
%end

function img = representative_img(pth)
img = [];
imgs = dir(fullfile(pth, [ '*']));
for i = 1 : numel(imgs)
    if imgs(i).name(1) == '.', continue; end
    if imgs(i).bytes < 1
        fprintf('>>>\t%s\t%s\n', pth, imgs(i).name);
        return;
    end
end

