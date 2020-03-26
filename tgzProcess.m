function tgzProcess(fnms, outDir)
%requires installation of latest dcm2niix: https://github.com/rordenlab/dcm2niix
%requires installation of nii_tool.m: part of https://github.com/xiangruili/dicm2nii
%Convert tgz file(s) to BIDs
% fnms: tgz file(s) of DICOMs to convert
%Versions
% 20200326 : includes support for MUSC .zip files
%Examples
% tgzProcess('20080601_1421.tgz');
% cd '/Volumes/TbStick/DICOM_JuliusSynolgy'
% d = dir('*.tgz')
% fnms = cellstr(strcat(pth,char(d.name)))
% tgzProcess(fnms);

tmpDir = '/Users/rorden/tst/org'; %fast cache drive
dcm2niix ='/usr/local/bin/dcm2niix';

if ~exist('outDir', 'var')
    fprintf('Select output master directory\n');
    outDir = uigetdir(pwd,'Select output master directory');
end
if ~exist(outDir, 'dir')
    outDir = pwd;
    error('outDir not specified, using %s\n', outDir);
end
%fnms = '/Users/rorden/tst/org/tgz/20070309_143935.tgz'; %20080521_143146.tgz
if ~exist('tmpDir', 'var') || ~exist(tmpDir, 'dir')
    tmpDir = fullfile(pwd,'tmp');
    if ~exist(tmpDir, 'dir') 
        mkdir(tmpDir)
    end
    warning('tmpDir not specified (this should be your fastest drive), using %s\n', tmpDir);
end
if ~exist(outDir, 'dir')
    error('Unable to output folder "%s"\n"', outDir);
end
dcm2niixCheckSub(dcm2niix);
if ~exist('fnms','var')
    fprintf('Select tgz file(s), or press "Cancel" to select DICOM folders\n');
    [files,pth] = uigetfile({'*.tgz';'*.*'},'Select the DICOM tgz(s)', 'MultiSelect', 'on');
    if ~isempty(files) && ~isa(files,'double')
        fnms = cellstr(strcat(pth,char(files)));
    else
        fprintf('Select folder that contains uncompressed DICOM files\n');
        fnms = uigetdir(pwd);
    end
end
if isempty(fnms), return; end
fnms= cellstr(fnms);
%get ids
ids = {};
exps = {};
for i = 1 : numel(fnms)
    fnm = deblank(fnms{i});
    [~,n,x] = fsl_filepartsSub(fnm);
    if strcmpi(x,'.zip') %MUSC uses .zip files
        prompt={sprintf('Enter ID (e.g. "M4002" for %s', n),'Enter the study (e.g. "POLAR"):'};    
    else
        prompt={sprintf('Enter ID (e.g. "M2002" for %s', n),'Enter the study (e.g. "LIME"):'};
    end
    name=sprintf('Details for %s', n);
    numlines=1;
    defaultanswer={'','POLAR'};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    if isempty(answer) || isempty(answer{1}) || isempty(answer{2})
        error('You must answer these questions');
    end
    ids{i} = answer{1};
    exps{i} = answer{2};
    
    %answer{1} = ['M' last4digitsof filname not including extension]
    %answer{2} = [LARC or POLAR depending on filename]
    
    
end
%process each file
for i = 1 : numel(fnms)
    fnm = deblank(fnms{i});
    [~,n,x] = fsl_filepartsSub(fnm);
    if n(1) == '.', continue; end %hidden file, parent folder, etc.
    fprintf('>> %d/%d >> %s\n', i,numel(fnms), n);
    outDirN = fullfile(outDir,ids{i});
    if ~exist(outDirN,'dir'), mkdir(outDirN); end
    outDirN = fullfile(outDirN,[exps{i},'_', n,'_',ids{i}]);
    if ~exist(outDirN,'dir'), mkdir(outDirN); end
    if isfolder(fnm)
        cmd = sprintf('%s -8 -z y -f %%t_%%2s_%%p -l y -i y -b y -t y -o %s %s', dcm2niix, outDirN, fnm);
        fprintf('Running \n %s\n', cmd);
        system(cmd,'-echo');
    elseif strcmpi(x,'.tgz') || strcmpi(x,'.tar.gz') || strcmpi(x,'.zip')
        tmpDirN = fullfile(tmpDir, ['temp',n]);
        if ~exist(tmpDirN,'dir'), mkdir(tmpDirN); end
        %dcms = untar(fnm, outDir);
        if strcmpi(x,'.zip') 
            unzip(fnm, tmpDirN);
        else
            untar(fnm, tmpDirN);
        end
        %convert DICOM -> Json
        cmd = sprintf('%s -8 -z y -f %%t_%%2s_%%p -l y -i y -b y -t y -o %s %s', dcm2niix, outDirN, tmpDirN);
        fprintf('Running \n %s\n', cmd);
        system(cmd,'-echo');
        rmdir(tmpDirN, 's');
    else
        error('Unknown file type %s\n');
    end
    %append prefixes, e.g. 'T1_', 'fMRI_', etc.
    nii_set_modality_prefix(outDirN,true);
    %combine all DTI scans of the same dimension/phaseEncoding direction
    nii_merge_dti_all(outDirN);
    %anonymize
    nii_deface_all(outDirN)
end
%end tgzProcess()

function [pth, nam, ext] = fsl_filepartsSub(fileName)
% a.nii.gz has the extension ".nii.gz" not ".nii"
[pth, nam, ext] = fileparts(fileName);
if (length(ext)==3)  && min((ext=='.gz')==1)
	[pth, nam, ext2] = fileparts( fullfile(pth, nam)); %remove .nii for .nii.gz
    ext = [ext2 ext];
end
%end fsl_filepartsSub()

function dcm2niixCheckSub(dcm2niix)
if ~exist(dcm2niix, 'file')
   error('Unable to find dcm2niix "%s"', dcm2niix);
end
minVer = 20190902;
cmd = sprintf('%s -f %%s', dcm2niix);
[~,cmdout] = system(cmd);
key = 'v1.0.';
pos = strfind(cmdout,key);
v = 19770703;
if ~isempty(pos)
	pos = pos(1)+numel(key);
	v = str2num(cmdout(pos:pos+7));
end
if v < minVer
	error('Update dcm2niix (from %d to at least %d) %s\n', v, minVer, dcm2niix);
end
%end dcm2niixCheckSub()
