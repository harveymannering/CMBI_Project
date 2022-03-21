close all;
clear all;

% Add helper functions to path
addpath('connectomes-data/Toolboxes/2016_01_16_BCT/');
addpath('connectomes-data/Toolboxes/NIfTI_20140122/');

FA_threshold = [0.1:0.1:0.8];
% Load nii files
filename = 'connectomes-data/Task1Data/tractor/diffusion/dti_FA.nii.gz';
result = load_nii(filename);

% Get the fa map from the results object
fa_images = result.img;
fa_images(isnan(fa_images)) = 0;

% Plot original image
slice_num = 15;
subplot(3,3,1);
imshow(fa_images(:,:,slice_num))
title('Original')

%% question 2
count_subvoxel = zeros([1,8]);
for i = 1:8
    % Threshold fa map
    threshold = 0.1 * i;
    fa_images_threshold = fa_images;
    fa_images_threshold(fa_images_threshold < threshold) = 0;
    
    % Plot
    subplot(3,3,1+i);
    imshow(fa_images_threshold(:,:,slice_num));

    slice_img = fa_images_threshold(:,:,slice_num);
    count_subvoxel(i) = length(slice_img(slice_img>0));
    title(threshold);
end

figure
hold on 
bar(FA_threshold,count_subvoxel)
xlabel('FA threshold')
ylabel('The number of voxels thresholded out')

%% question 3
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

display_graphs(density,char_path,efficiency,mean_cluster_coeff,'FA Threshold');

%% question 4
parcellation_func = load_nii('connectomes-data/Task1Data/tractor/functional/parcellation.nii.gz');
data_func = load_nii('connectomes-data/Task1Data/tractor/functional/data.nii.gz');

% Load cortical region information
file_id = fopen('connectomes-data/Task1Data/tractor/functional/parcellation.lut');
region = textscan(file_id, '%d %s %s %s %s %s %s', 'Headerlines', 3); % start from foutrh line
region_ids = region{1}; 
[num_regions,~] = size(region_ids);

% now we are going to compute the average time series for each cortical region
% the size of rs-fMRI is 64×64×30×15, and the fourth dimension is time 
avg_time_series = zeros(num_regions, 15);
for i = 1:num_regions
    id = region_ids(i);

    % use mask to select voxels 
    cortical_region = (parcellation_func.img == id); 
    voxels = data_func.img(cortical_region(:,:,:,ones(1,15)));
    [size_voxels,~] = size(voxels); 
    region_voxels = reshape(voxels, [size_voxels/15, 15]); 
    
    % average time series of this cortical region's all voxel
    avg_time_series(i,:) = mean(region_voxels, 1);
end

% corshrink
for i = 1:8
    corshrink_matrix{i} = corshrink(avg_time_series',FA_threshold(i));
end

% resulting matrix at a correlation of 0.1, binarise it
for t = 1:8
    corshrink_new = corshrink_matrix{t};
    % binarise it
    corshrink_new(corshrink_new < 0.1) = 0;
    corshrink_new(corshrink_new >= 0.1) = 1;
    binary_corshrink{t} = corshrink_new;
end

%  graph metrics 
for i = 1:8
    % function like q3
    [cor_density{i}(1), cor_density{i}(2), cor_density{i}(3)] = density_und(binary_corshrink{i}); 
    cor_dist_mat{i} = distance_bin(binary_corshrink{i}); 
    [cor_char_path{i}(1), cor_char_path{i}(2)]  = charpath(cor_dist_mat{i},0,0); 
    cor_efficiency{i} = 1 / efficiency_bin(binary_corshrink{i}); 
    cor_mean_cluster_coeff{i} = mean(clustering_coef_bu(binary_corshrink{i})); 
end

display_graphs(cor_density,cor_char_path,cor_efficiency,cor_mean_cluster_coeff,'Lambda');

% discarding negative correlations with absolute value greater than 0.1
for t = 1:8
    corshrink_new = corshrink_matrix{t};
    corshrink_new = (corshrink_new > -0.1);
    discard_neg_corshrink{t} = corshrink_new;
end

for i = 1:8
    [discard_neg_density{i}(1), discard_neg_density{i}(2), discard_neg_density{i}(3)] = density_und(discard_neg_corshrink{i});
    discard_neg_dist_mat{i} = distance_bin(discard_neg_corshrink{i}); 
    [discard_neg_char_path{i}(1), discard_neg_char_path{i}(2)]  = charpath(discard_neg_dist_mat{i},0,0); 
    discard_neg_efficiency{i} = 1 / efficiency_bin(discard_neg_corshrink{i}); 
    discard_neg_mean_cluster_coeff{i} = mean(clustering_coef_bu(discard_neg_corshrink{i})); 
end

display_graphs(discard_neg_density,discard_neg_char_path,discard_neg_efficiency,discard_neg_mean_cluster_coeff,'Lambda');

% function definition
% mannually setting lambda
function [Rhat] = corshrink(x,input_lambda)
    % Eqn on p4 of Schafer and Strimmer 2005
    [~, p] = size(x);
    sx = makeMeanZero(x); sx = makeStdOne(sx); % convert S to R
    [r, ~] = varcov(sx);
    lambda = input_lambda;
    lambda = min(lambda, 1); lambda = max(lambda, 0);
    Rhat = (1-lambda)*r;
    Rhat(logical(eye(p))) = 1;
end

function [S, VS] = varcov(x)
    % s(i,j) = cov X(i,j)
    % vs(i,j) = est var s(i,j)
    [n,p] = size(x);
    xc = makeMeanZero(x); 
    S = cov(xc);
    XC1 = repmat(reshape(xc', [p 1 n]), [1 p 1]); % size p*p*n !
    XC2 = repmat(reshape(xc', [1 p n]),  [p 1 1]); % size p*p*n !
    VS = var(XC1 .* XC2, 0,  3) * n/((n-1)^2);
end

function xc = makeMeanZero(x)
% make column means zero
    [n,~] = size(x);
    m = mean(x);
    xc = x - ones(n, 1)*m; 
end

function xc = makeStdOne(x)
    % make column  variances one
    [n,~] = size(x);
    sd = ones(n, 1)*std(x);
    xc = x ./ sd; 
end

function display_graphs(density,char_path,efficiency,mean_cluster_coeff,x_label)
    FA_threshold = [0.1:0.1:0.8];
    density_data = zeros([1,8]);
    mean_shortest_path = zeros([1,8]);
    efficiency_data = zeros([1,8]);
    mean_cluster_coeff_data = zeros([1,8]);
    for i = 1:8
        density_data(i) = density{i}(1);
        mean_shortest_path(i) = char_path{i}(1);
        efficiency_data(i) = efficiency{i};
        mean_cluster_coeff_data(i) = mean_cluster_coeff{i};
    end

    figure
    hold on 
    plot(FA_threshold,mean_shortest_path)
    hold on 
    plot(FA_threshold,efficiency_data)
    xlabel(x_label)
    ylabel('Value')
    legend('mean shortest path','efficiency')

    figure
    plot(FA_threshold,density_data)
    hold on 
    plot(FA_threshold,mean_cluster_coeff_data)
    xlabel(x_label)
    ylabel('Value')
    legend('edge density','mean clustering coefficient')
end