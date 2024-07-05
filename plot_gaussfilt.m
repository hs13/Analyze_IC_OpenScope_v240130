blankICwcfg0 = double(imread('/Users/hyeyoung/Documents/CODE/Display_IC/visOpenScope/visICtxiwcfg0/100000.tif'));
halfhorz = blankICwcfg0(round(size(blankICwcfg0,1)/2),:);
cumhalfhorz = cumsum(halfhorz);
npix16deg = nnz(cumhalfhorz==cumhalfhorz(round(length(cumhalfhorz)/2)))+1;
figure; plot(cumhalfhorz)

fs=14;
impath = '/Users/hyeyoung/Documents/CODE/Display_IC/visOpenScope/visICtxiwcfg1/';
imlist = {'110106', '110111', '110506', '110511'};
figure
for ii = 1:numel(imlist)
A = imread([impath imlist{ii} '.tif']);
%B = imgaussfilt(A,npix16deg/2);
B = imgaussfilt(A, 'FilterSize', round(npix16deg/2)*2+1);
subplot(2,2,ii);
imagesc(B)
title(imlist{ii})
axis equal
end

figure
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', 'gauss filter std 4 deg', 'edgecolor', 'none', 'FontSize', fs)
for ii = 1:numel(imlist)
A = imread([impath imlist{ii} '.tif']);
B = imgaussfilt(A,npix16deg/4);
subplot(2,2,ii);
imagesc(B)
title(imlist{ii}, 'FontSize', fs)
axis equal
axis off
colormap gray
end

figure
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', 'gauss filter std 5 deg', 'edgecolor', 'none', 'FontSize', fs)
for ii = 1:numel(imlist)
A = imread([impath imlist{ii} '.tif']);
B = imgaussfilt(A,(20/16)*npix16deg/4);
subplot(2,2,ii);
imagesc(B)
title(imlist{ii}, 'FontSize', fs)
axis equal
axis off
colormap gray
end

figure
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', 'gauss filter std 2.5 deg', 'edgecolor', 'none', 'FontSize', fs)
for ii = 1:numel(imlist)
A = imread([impath imlist{ii} '.tif']);
B = imgaussfilt(A,(10/16)*npix16deg/4);
subplot(2,2,ii);
imagesc(B)
title(imlist{ii}, 'FontSize', fs)
axis equal
axis off
colormap gray
end

% Neuropixels monitor spanned 120° × 95° of visual space 
% Neuropixels n=18 exclusively center-responsive neurons out of 1,804 % 0.0100
% 2p data n=310 exclusively center-responsive neurons out of 18,576 V1 layer 2/3 neurons % 0.0167

((16-10)/(90-10))^2 % 0.0056
(pi*(16-10)^2)/((95-10)*(120-10)) % 0.0121
