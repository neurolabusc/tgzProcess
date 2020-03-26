function nii_show(imgs, png)
%Show a mosaic of different NIfTI images
% imgs : filenames of NIfTI images
% png  : optional, create png file with this filename
% nii_show('T1.nii.gz')
% nii_show({'T1.nii.gz', 'T2.nii.gz'})

if ~exist('imgs', 'var')
   [files,pth] = uigetfile({'*.nii;*.gz;';'*.*'},'Select the Image[s]', 'MultiSelect', 'on'); 
   imgs = strcat(pth, files);
end
imgs = cellstr(imgs);
if numel(imgs) < 1, return; end
h = findobj('type','figure','name','mat2ortho'); %re-use if available
if isempty(h), h = figure('Name','mat2ortho','NumberTitle','off'); end; %make sure we do not use SPM graphics
figure(h); %make current
clf;
XYZmm = [0;0;0];
for i = 1: numel(imgs)
    fnm = imgs{i};
    nii = nii_tool('load', fnm);
    nii = nii_reorientSub(nii, true);
    nii.hdr.dim = nii.hdr.dim(1:3);
    nii.img = nii.img(:,:,:,1);
    [~,nam,~] = fileparts(fnm);
    plotOrthoSub(nii.hdr, nii.img, XYZmm, i, numel(imgs), nam);
end
if ~exist('png', 'var'), return; end
[~,~,x] = fileparts(png);
if ~strcmpi(x,'.png')
    png = [png, '.png'];
end
print(png,'-dpng')

function plotOrthoSub(hdr, img, XYZmm, Slot, NumSlots, Caption)
xhair = mm2voxSub(hdr, XYZmm);
set(gcf,'color','w');
%img = img - min(img(:)); %set minimum to zero
%img(img < 0) = 0;
img = double(img);
threshLo = min(img(:));
threshHi = mean(img(:));
if  numel(img) > 1000 
    imgS = sort( img(isfinite(img(:))));
    imgS = imgS(imgS ~= 0);
    pct = round(numel(imgS) * 0.98); %brightest 2%
    threshHi = imgS(pct);
    pct = round(numel(imgS) * 0.02); %darkest 2%
    threshLo = imgS(pct);
end
if threshLo == threshHi %e.g. binary image map
    threshLo = min(img(:));
    threshHi = max(img(:));   
end
img(~isfinite(img)) = 0;
img = img - threshLo; %translate so darkest voxel is 0
thresh = threshHi - threshLo; 
img(img > thresh) = thresh;
img = 63 * (img / thresh); %matlab color scheme have 64 indices "size(bone)"
sz = size(img);
colormap(gray); %colormap(bone)
ax = img(:,:,xhair(3));
ax(xhair(1),:) = 128;
ax(:,xhair(2)) = 128;
cor = squeeze(img(:,xhair(2),:));
cor(xhair(1),:) = 128;
cor(:,xhair(3)) = 128;
sag = squeeze(img(xhair(1),:,:));
sag(xhair(2),:) = 128;
sag(:,xhair(3)) = 128;
im = zeros( sz(3)+sz(2), sz(1)+sz(2));
im(1:sz(1),1:sz(2)) = ax;
im(1:sz(1),(1:sz(3))+sz(2)) = cor;
im((1:sz(2))+sz(1),(1:sz(3))+sz(2)) = sag;
%scale output
rows = ceil(sqrt(NumSlots)); %e.g. 2..4 items shown in 2x2 mosaic, 5..9 in 3x3
scale = 1/rows;
col = mod(Slot-1,rows) * scale;
row = floor((NumSlots - Slot)/rows) * scale; %top->bottom, for bottom->top: row = floor((Slot-1)/rows) * scale;
plotImgSub ( flipud(im'), Caption, col, row, scale, scale);
%plotOrthoSub

function xhair = mm2voxSub(hdr, XYZmm)
mat = [hdr.srow_x; hdr.srow_y; hdr.srow_z; 0 0 0 1];
mInv = inv(mat);
xhair = mInv * [XYZmm; 1]; %convert from voxels to millimeters

xhair = round(xhair(1:3))';
%xhair(xhair < 1) = 1;
%xhair = min(xhair, hdr.dim);
%end mm2voxSub

function plotImgSub ( Img, Caption, X, Y, wid, ht)
subplot('Position',[X Y wid ht]); %width, height
image((Img));
set(gca,'XTickLabel', [],'XTick',[],'YTick',[]);
axis image
w = wid/2;
h = ht/2;
annotation('textbox',[X+w Y w h],'String',Caption,'FontSize',18,'fontn','Arial', 'color','red', 'BackgroundColor', 'white', 'LineStyle','none')
%end plotImgSub()

function XYZmm = getCenterOfIntensitySub(hdr, img)
XYZmm = ones(3,1);
img = img - min(img(:));
img(isnan(img)) = 0;
%find center of mass in each dimension (total mass divided by weighted location of mass
% img = [1 2 1; 3 4 3];
sumTotal = sum(img(:));
coivox = ones(4,1);
coivox(1) = sum(sum(sum(img,3),2)'.*(1:size(img,1)))/sumTotal; %dimension 1
coivox(2) = sum(sum(sum(img,3),1).*(1:size(img,2)))/sumTotal; %dimension 2
coivox(3) = sum(squeeze(sum(sum(img,2),1))'.*(1:size(img,3)))/sumTotal; %dimension 3
XYZmm = hdr.mat * coivox; %convert from voxels to millimeters
XYZmm = XYZmm(1:3);
%end setCenterOfIntensitySub()

%following functiosn from Xiangrui Li nii_viewer.m
% . retain his BSD 2-Clause license
%https://github.com/xiangruili/dicm2nii
function [nii, perm, flp] = nii_reorientSub(nii, leftHand)
% reorient nii to diagnal major
[R, frm] = nii_xform_mat(nii.hdr);
dim = nii.hdr.dim(2:4);
pixdim = nii.hdr.pixdim(2:4);
[R, perm, flp] = reorient(R, dim, leftHand);
fps = bitand(nii.hdr.dim_info, [3 12 48]) ./ [1 4 16];
if ~isequal(perm, 1:3)
    nii.hdr.dim(2:4) = dim(perm);
    nii.hdr.pixdim(2:4) = pixdim(perm);
    nii.hdr.dim_info = [1 4 16] * fps(perm)' + bitand(nii.hdr.dim_info, 192);
    nii.img = permute(nii.img, [perm 4:8]);
end
sc = nii.hdr.slice_code;
if sc>0 && flp(fps==3)
    nii.hdr.slice_code = sc+mod(sc,2)*2-1; % 1<->2, 3<->4, 5<->6
end
if isequal(perm, 1:3) && ~any(flp), return; end
if frm(1) == nii.hdr.sform_code % only update matching form
    nii.hdr.srow_x = R(1,:);
    nii.hdr.srow_y = R(2,:);
    nii.hdr.srow_z = R(3,:);
end
if frm(1) == nii.hdr.qform_code
    nii.hdr.qoffset_x = R(1,4);
    nii.hdr.qoffset_y = R(2,4);
    nii.hdr.qoffset_z = R(3,4);
    R0 = normc(R(1:3, 1:3));
    dcm2quat = dicm2nii('', 'dcm2quat', 'func_handle');
    [q, nii.hdr.pixdim(1)] = dcm2quat(R0);
    nii.hdr.quatern_b = q(2);
    nii.hdr.quatern_c = q(3);
    nii.hdr.quatern_d = q(4);
end
for i = find(flp), nii.img = flip(nii.img, i); end

% Return xform mat and form_code: form_code may have two if not to ask_code
function [R, frm] = nii_xform_mat(hdr, ask_code)
% [R, form] = nii_xform_mat(hdr, asked_code);
% Return the transformation matrix from a NIfTI hdr. By default, this returns
% the sform if available. If the optional second input, required form code, is
% provided, this will try to return matrix for that form code. The second
% optional output is the form code of the actually returned matrix.
fs = [hdr.sform_code hdr.qform_code]; % sform preferred
if fs(1)==fs(2), fs = fs(1); end % sform if both are the same
f = fs(fs>=1 & fs<=4); % 1/2/3/4 only
if isempty(f) || ~strncmp(hdr.magic, 'n', 1) % treat it as Analyze
    frm = 0;
    try % try spm style Analyze
        [pth, nam, ext] = fileparts(hdr.file_name);
        if strcmpi(ext, '.gz'), [~, nam] = fileparts(nam); end
        R = load(fullfile(pth, [nam '.mat']));
        R = R.M;
    catch % make up R for Analyze: suppose xyz order with left storage 
        R = [diag(hdr.pixdim(2:4)) -(hdr.dim(2:4).*hdr.pixdim(2:4)/2)'; 0 0 0 1];
        R(1,:) = -R(1,:); % use left handed
    end
    return;
end
if numel(f)==1 || nargin<2 || isempty(ask_code) % only 1 avail or no ask_code
    frm = f;
else % numel(f) is 2, numel(ask_code) can be 1 or 2
    frm = f(f == ask_code(1));
    if isempty(frm) && numel(ask_code)>1, frm = f(f == ask_code(2)); end
    if isempty(frm) && any(f==2), frm = 2; end % use confusing code 2
    if isempty(frm), frm = f(1); end % no match to ask_code, use sform
end
if frm(1) == fs(1) % match sform_code or no match
    R = [hdr.srow_x; hdr.srow_y; hdr.srow_z; 0 0 0 1];
else % match qform_code
    R = quat2R(hdr);
end

% Reorient 4x4 R
function [R, perm, flp] = reorient(R, dim, leftHand)
% [R, perm, flip] = reorient(R, dim, leftHand)
% Re-orient transformation matrix R (4x4), so it will be diagonal major and
% positive at diagonal, unless the optional third input is true, which requires
% left-handed matrix, where R(1,1) will be negative. 
% The second input is the img space dimension (1x3). 
% The perm output, like [1 2 3] or a permutation of it, indicates if input R was
% permuted for 3 axis. The third output, flip (1x3 logical), indicates an axis 
% (AFTER perm) is flipped if true.
a = abs(R(1:3,1:3));
[~, ixyz] = max(a);
if ixyz(2) == ixyz(1), a(ixyz(2),2) = 0; [~, ixyz(2)] = max(a(:,2)); end
if any(ixyz(3) == ixyz(1:2)), ixyz(3) = setdiff(1:3, ixyz(1:2)); end
[~, perm] = sort(ixyz);
R(:,1:3) = R(:,perm);
flp = R([1 6 11]) < 0; % diag(R(1:3, 1:3))
if nargin>2 && leftHand, flp(1) = ~flp(1); end
rotM = diag([1-flp*2 1]);
rotM(1:3, 4) = (dim(perm)-1) .* flp; % 0 or dim-1
R = R / rotM; % xform matrix after flip

% normalize columns
function v = normc(M)
v = bsxfun(@rdivide, M, sqrt(sum(M .* M)));
% v = M ./ sqrt(sum(M .* M)); % since 2016b