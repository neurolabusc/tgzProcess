function names = nii_deface_crop(fnms, isSaveBitmap, isReduceVox)
% Face strip images, removes neck fat reduces volume to reasonable bounding box
% FORMAT names = spm_deface(fnms)
% fnms   - cell array of NIfTI file names
% saveBmp - if true image saved for subsequent review
%
% names        - cell array of de-faced images
%
% Extends spm_deface by John Ashburner Copyright (C) 2013-2014 Wellcome Trust Centre for Neuroimaging
%         $Id: spm_deface.m 6086 2014-07-03 16:08:44Z guillaume $


if ~exist('fnms','var')
    [P, sts] = spm_select(Inf,'image','Select images to strip the face from');
    if ~sts, return; end
    %P = 'T1_CT184.nii.gz';z
elseif isstruct(fnms)
    P = fnms.images;
else
    P = fnms;
end
if ~exist('isSaveBitmap','var')
    isSaveBitmap = true;
end
if ~exist('isReduceVox','var')
    isReduceVox = true;
end

P = cellstr(P);
names = cell(size(P));
tpm = spm_load_priors8(spm_vol(fullfile(spm('Dir'),'tpm','TPM.nii')));
for i=1:numel(P)
    names{i} = deface(P{i},tpm, isReduceVox, isSaveBitmap);
end
%end spm_deface_crop()

function outname = deface(fnm,tpm, isReduceVox, isSaveBitmap)
isOverwrite = true; 
isGz = false;
if endsWith(upper(fnm),'.GZ')
    isGz = true;
    nm = gunzip(fnm);
    fnm = nm{:};
end
nul     = [0 -1.1 0.98 120; ... %eyeballs
    0 1.0 0.8 132; ... %posterior neck
    0 0 1 80; ... %bottom slices
    0 0 -1 100; ... %top slices
   0 1 0 130; ... %posterior slices
    0 -1 0 100; ... %anterior slices
    1 0 0 90; ... %left slices
    -1 0 0 90]; %right slices
%spm_maff8()
Nii     = nifti(fnm);
if startsWith(Nii.descrip,'^') || (size(Nii.dat,4) > 1)
    if (size(Nii.dat,4) > 1)
        fprintf('Skipping file: 4D image %s\n', fnm);
    else
        fprintf('Skipping file: already cropped %s\n', fnm);
    end
   if isGz
    delete(fnm);
   end
   outname = [];
   return;
end
M       = spm_maff8(fnm,4,20,tpm,[],'mni');
d       = [size(Nii.dat) 1];
[i,j,k] = ndgrid(1:d(1),1:d(2),1:d(3));
nul1    = nul(1,:)*M*Nii.mat;
msk     = nul1(1)*i + nul1(2)*j + nul1(3)*k + nul1(4) < 0;
for o = 2: 8%size(nul,1)
    nul1    = nul(o,:)*M*Nii.mat;
    mskj     = nul1(1)*i + nul1(2)*j + nul1(3)*k + nul1(4) < 0;
    msk = msk | mskj ;
end
inDim = Nii.dat.dim;
if ~isReduceVox
    mnD = [1,1,1];
    mxD = inDim(1:3);
else
    nx = inDim(1);
    ny = inDim(2);
    nz = inDim(3);
    x = sum(sum(msk,3)'); %X
    x = (x ~= (ny*nz));
    mnD(1) = find(x, 1, 'first');
    mxD(1) = find(x, 1, 'last');
    y = sum(sum(msk,3)); %y
    y = (y ~= (nx*nz));
    mnD(2) = find(y, 1, 'first');
    mxD(2) = find(y, 1, 'last');

    z = sum(sum(msk,2)); %Z
    z = (z ~= (nx*ny));
    mnD(3) = find(z, 1, 'first');
    mxD(3) = find(z, 1, 'last');
    vxD =mxD-mnD+1;
    for j = 1 : 3
       if mod(vxD(j),2)
        if (mnD(j) > 1)
            mnD(j) = mnD(j) -1;
        elseif (mxD(j) < inDim(j))
            mxD(j) = mxD(j) +1;
        end
       end
    end
end
vxD =mxD-mnD+1;

vx = prod(vxD);
if vx <= 1, error('image not coregistered'); end;
pct = 100* vx/(inDim(1)*inDim(2)*inDim(3));
[~,n] = fileparts(Nii.dat.fname);
if (pct < 100)
	fprintf('%s cropping %dx%dx%d -> %dx%dx%d (%g%%) %s\n', mfilename, inDim(1), inDim(2), inDim(3), mxD(1)-mnD(1)+1,mxD(2)-mnD(2)+1, mxD(3)-mnD(3)+1, pct, n);
end
v2m = Nii.mat; %voxel2mm transform
m2v=inv(v2m); %mm2voxel transform

origin= (mnD-1)*v2m(1:3,1:3)' + v2m(1:3,4)';

outname   = spm_file(Nii.dat.fname,'prefix','crop_');
Noo2     = Nii;
Noo2.descrip = ['^', Noo2.descrip];
Noo2.mat(1:3,4) = origin;

Noo2.dat.fname = outname;
Noo2.dat.dim = vxD;
create(Noo2);
for k=1:size(Nii.dat,6),
    for j=1:size(Nii.dat,5),
        for i=1:size(Nii.dat,4),
            F       = Nii.dat(:,:,:,i,j,k);
            F(msk)  = NaN;
            Noo2.dat(:,:,:) = F(mnD(1):mxD(1), mnD(2):mxD(2),mnD(3):mxD(3));
        end
    end
end
if isSaveBitmap
    %plotOrthoSub(Noo2.dat(:,:,:), outname);
    nii_show(outname, outname)
end
if isOverwrite
    movefile(outname,fnm);
    outname = fnm;
end
if ~isOverwrite && isGz %delete decompressed duplicate
    delete(fnm);
end
if isGz
   nm = gzip(outname);
   delete(outname);
   outname = nm{:};
end

%end defaceSub()

% function plotOrthoSub(img, Caption)
% xhair = [size(img,1)/2, size(img,2)/2, size(img,3)/2 ];
% xhair = round([size(img,1)/2, size(img,2)/2, size(img,3)/2 ]);
% if min(xhair(:)) < 1, return; end
% set(gcf,'color','w');
% threshLo = min(img(:));
% threshHi = mean(img(:));
% if  numel(img) > 1000
%     imgS = sort( img(isfinite(img(:))));
%     imgS = imgS(imgS ~= 0);
%     pct = round(numel(imgS) * 0.98); %brightest 2%
%     threshHi = imgS(pct);
%     pct = round(numel(imgS) * 0.02); %darkest 2%
%     threshLo = imgS(pct);
% end
% if threshLo == threshHi %e.g. binary image map
%     threshLo = min(img(:));
%     threshHi = max(img(:));
% end
% img(~isfinite(img)) = 0;
% img = img - threshLo; %translate so darkest voxel is 0
% thresh = threshHi - threshLo;
% img(img > thresh) = thresh;
% img = 63 * (img / thresh); %matlab color scheme have 64 indices "size(bone)"
% sz = size(img);
% colormap(gray); %colormap(bone)
% ax = img(:,:,xhair(3));
% ax(xhair(1),:) = 128;
% ax(:,xhair(2)) = 128;
% cor = squeeze(img(:,xhair(2),:));
% cor(xhair(1),:) = 128;
% cor(:,xhair(3)) = 128;
% sag = squeeze(img(xhair(1),:,:));
% sag(xhair(2),:) = 128;
% sag(:,xhair(3)) = 128;
% im = zeros( sz(3)+sz(2), sz(1)+sz(2));
% im(1:sz(1),1:sz(2)) = ax;
% im(1:sz(1),(1:sz(3))+sz(2)) = cor;
% im((1:sz(2))+sz(1),(1:sz(3))+sz(2)) = sag;
% %scale output
% plotImgSub ( flipud(im'), Caption,0, 0, 1, 1);
% %print('-dpsc', '-append', [Caption, '.ps']);
% print([Caption, '.png'],'-dpng')
% %end plotOrthoSub
% 
% function plotImgSub ( Img, Caption, X, Y, wid, ht)
% subplot('Position',[X Y wid ht]); %width, height
% image((Img));
% set(gca,'XTickLabel', [],'XTick',[],'YTick',[]);
% axis image
% w = wid/2;
% h = ht/2;
% annotation('textbox',[X+w Y w h],'String',Caption,'FontSize',18,'fontn','Arial', 'color','red', 'BackgroundColor', 'white', 'LineStyle','none')
% %end plotImgSub()

