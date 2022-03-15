% Add helper functions to path
addpath('connectomes-data/Toolboxes/2016_01_16_BCT/');
addpath('connectomes-data/Toolboxes/NIfTI_20140122/');
addpath('connectomes-data/Toolboxes/covshrink-kpm/');

% Load nii files
filename = 'connectomes-data/Task1Data/tractor/diffusion/dti_FA.nii.gz';
result = load_nii(filename);
disp(result)

% Get the fa map from the results object
fa_images = result.img;
% fa_images(isnan(fa_images)) = 0;
% we would figure out it later in question 3

% Plot original image
slice_num = 15;
subplot(3,3,1);
imshow(fa_images(:,:,slice_num))
title('Original')

for i = 1:8
    % Threshold fa map
    threshold = 0.1 * i;
    fa_images_threshold = fa_images;
    fa_images_threshold(fa_images_threshold < threshold) = 0;
    
    % Plot
    subplot(3,3,1+i);
    imshow(fa_images_threshold(:,:,slice_num));
    title(threshold);
end

% To improve visualization we could:
%  - Rewindow
%  - Colourize


connectome_file = {'FA.1_graph.csv'; 'FA.2_graph.csv'; 'FA.3_graph.csv';...
            'FA.4_graph.csv'; 'FA.5_graph.csv'; 'FA.6_graph.csv'; ...
            'FA.7_graph.csv'; 'FA.8_graph.csv'};
for i = 1:8
    fa_file_path = ['connectomes-data/Task1Data/connectomes/', connectome_file{i}];
    % skip the first three lines 
    struct_connect{i} = csvread(fa_file_path,3,0);

    % this function return density, number of vertices, and number of edges in 3
    [density{i}(1), density{i}(2), density{i}(3)] = density_und(struct_connect{i}); 

    % mean shortest path
    % distance matrices 
    distance_mat{i} = distance_bin(struct_connect{i}); 
    % this function return network characteristic path length and network global efficiency
    % Characteristic path length is defined here as the mean shortest path length between all pairs of nodes
    % About missing edges and disconnected, it will be infinite value, it
    % would be zero as third argument setting
    [char_path{i}(1), char_path{i}(2)]  = charpath(distance_mat{i},0,0); 

    % efficiency
    efficiency{i} = 1 / efficiency_bin(struct_connect{i}); 
    
    % mean clustering coefficient
    mean_cluster_coeff{i} = mean(clustering_coef_bu(struct_connect{i})); 
end

