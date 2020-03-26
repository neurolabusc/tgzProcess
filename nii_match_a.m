function nii_match_a (basepth)
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
if ~exist('basepth','var')
    basepth = pwd;
end
cd(basepth);
subjs = dir('M*');
tmpdir = fullfile(pwd, 'tmp');
if ~exist(tmpdir, 'dir')
   mkdir(tmpdir) 
end
for s = 1: numel(subjs)
    if ~subjs(s).isdir, continue; end
    if ~isempty(strfind(subjs(s).name,'_')), continue; end
    fprintf('---');
    subjpth = fullfile(basepth, subjs(s).name);
    fprintf('%s\n', subjs(s).name);
    cd(subjpth);
    outdir = fullfile(tmpdir, subjs(s).name);
    mkdir(outdir);
    %simages = [];
    visits = dir('*');
    modalities = {'T1','T2','DTI','fMRI'};
    for v = 1: numel(visits)
        if ~visits(v).isdir, continue; end
        if visits(v).name(1) == '.', continue; end
        fprintf(' %s\n', visits(v).name);
        visitpth = fullfile(subjpth, visits(v).name);
        cd(visitpth);
        img = representative_img();
        if isempty(img)
            fprintf(' ??? No scans for %s\n', visits(v).name);
            continue; 
        end
        %fprintf('  %s\n', image);
        [~,n,x] = fileparts(img);
        innm = fullfile(visitpth, [n, x]);
        outnm = fullfile(outdir, [n, x]);
        fprintf('  %s -> %s\n', innm, outnm);
        copyfile(innm, outnm);
        if strcmpi(x, '.gz')
            gunzip(outnm);
            delete(outnm);
        end
        %simages{end+1} = fullfile(visitpth, image);
    end
    %if isempty(simages), continue; end
    %for i = 1 : numel(simages);
    %    [a,b,c] = fileparts(simages{i});
    %    outnm = fullfile(outdir, [b,c]);
    %    copyfile(simages{1}, outnm);
    %end
end
cd(basepth);
%end

function img = representative_img()
modalities = {'T1','T2','DTI','fMRI'};
for i = 1 : numel(modalities)
    images = dir([modalities{i}, '*.nii.gz']);
    if ~isempty(images)
       img = images(1).name;
       return
    end
    images = dir([modalities{i}, '*.nii']);
    if ~isempty(images)
       img = images(1).name;
       return
    end
end
