% Add helper functions to path
addpath('/Users/harveymannering/Downloads/connectomes-data/Toolboxes/2016_01_16_BCT/');
addpath('/Users/harveymannering/Downloads/connectomes-data/Toolboxes/NIfTI_20140122/');
addpath('/Users/harveymannering/Downloads/connectomes-data/Toolboxes/covshrink-kpm/');

% Load nii files
filename = '/Users/harveymannering/Downloads/connectomes-data/Task1Data/tractor/diffusion/dti_FA.nii.gz';
result = load_nii(filename);
disp(result)

% Get the fa map from the results object
fa_images = result.img;
fa_images(isnan(fa_images)) = 0;

% Plot original image
slice_num = 15;
subplot(1,9,1);
imshow(fa_images(:,:,slice_num))
title('Original')

for i = 1:8
    % Threshold fa map
    threshold = 0.1 * i;
    fa_images_threshold = fa_images;
    fa_images_threshold(fa_images_threshold < threshold) = 0;
    
    % Plot
    subplot(1,9,1+i);
    imshow(fa_images_threshold(:,:,slice_num));
    title(threshold);
end

% To improve visualization we could:
%  - Rewindow
%  - Colourize