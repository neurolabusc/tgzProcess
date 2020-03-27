function nii_match_b (basepth)
%Normalize all scans from an individual and display
%
%Images are basepth/Participant/Session/*.nii
%For example:
% M2002/
%       T1_201801202.nii
%       fMRI.nii
% M2012/
%    T1_201801202.nii
%    fMRI_201901202.nii
%
if ~exist('basepth','var')
    basepth = pwd;
end
basepth = '/Users/chris/Downloads/harvestZip/tgzProcess/tmp'
%cd(basepth);
%subjs = dir('M*');
subjs = dir(fullfile(basepth, 'M*'));
fnms={subjs.name};
[~,idx]=sort(fnms);
subjs=subjs(idx);
for s = 1: numel(subjs)
    if ~subjs(s).isdir, continue; end
    if ~isempty(strfind(subjs(s).name,'_')), continue; end
    subjpth = fullfile(basepth, subjs(s).name);
    fprintf('%d/%d\t%s\n', s, numel(subjs), subjs(s).name);
    %cd(subjpth);
    %simages = [];
    %visits = dir('*.nii');
    visits = dir(fullfile(subjpth, '*.nii'));
    for v = 1: numel(visits)
        if visits(v).isdir, continue; end
        if visits(v).name(1) == '.', continue; end
        if visits(v).name(1) == 'w', continue; end %already warped
        visitname = fullfile(subjpth, visits(v).name);
        warpname = fullfile(subjpth, ['w', visits(v).name]);
        if exist(warpname, 'file'), continue; end
        %fprintf(' %s\n', visits(v).name);
        affine_norm(visitname);
    end
    visits = dir(fullfile(subjpth, 'w*.nii'));
    imgs = [];
    for v = 1: numel(visits)
        if visits(v).isdir, continue; end
        warpname = fullfile(subjpth, visits(v).name);
        imgs = [imgs, cellstr(warpname)];
    end
    nii_show(imgs,subjs(s).name);
    
end
%cd(basepth);
%end

function affine_norm(fnm)
[~,nm] = fileparts(fnm);
if startsWith(nm,'w','IgnoreCase',true)
   return; 
end
template = fullfile(spm('Dir'),'toolbox','OldNorm','EPI.nii,1');
if startsWith(nm,'T1_','IgnoreCase',true)
    template = fullfile(spm('Dir'),'toolbox','OldNorm','T1.nii,1');
end
if startsWith(nm,'T2_','IgnoreCase',true)
    template = fullfile(spm('Dir'),'toolbox','OldNorm','T2.nii,1');
end
spm('defaults', 'FMRI');
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source = cellstr([fnm,',1']);
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = cellstr([fnm,',1']);
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = cellstr(template);
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 0;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg = 1000;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [-78 -112 -70; 78 76 85];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = [2 2 2];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
spm_jobman('run',matlabbatch);
%affine_norm()
