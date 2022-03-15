% Add helper functions to path
addpath('/Users/harveymannering/Downloads/connectomes-data/Toolboxes/2016_01_16_BCT/');
addpath('/Users/harveymannering/Downloads/connectomes-data/Toolboxes/NIfTI_20140122/');
addpath('/Users/harveymannering/Downloads/connectomes-data/Toolboxes/covshrink-kpm/');

% Load csv files
s = zeros(68, 68, 19);
f = zeros(68, 68, 19);
for i = 1:19
    s(:,:,i) = csvread(strcat('/Users/harveymannering/Downloads/connectomes-data/Task2Data/', int2str(31+i), '_WFA_68.csv'), 1, 0);
    f(:,:,i) = csvread(strcat('/Users/harveymannering/Downloads/connectomes-data/Task2Data/', int2str(31+i), '_rsfMRI_68.csv') , 1, 0);
end
% Calucluate indirect structural connectivity matrix
t = zeros(size(s));

% Iterate through every patient p
for p = 1:size(s,3)
    % Iterate through (and calculate) every element of t
    for t_i = 1:size(s, 1)
        for t_j = 1:size(s, 2)
    
            % Find the greatest minimum weight in all available two-step chains
            greatest_min_weight = -realmax;
            for k = 1:size(s, 1)
                min_weight = min(s(t_i, k, p), s(k, t_j, p));
                if (min_weight > greatest_min_weight)
                    greatest_min_weight = min_weight;
                end
            end
    
            % Set matrix element of the indirect structural connectivity matrix
            t(t_i, t_j, p) = greatest_min_weight;
        end 
    end 
end 

%% Fit Model 1: f = alpha + beta * s
alpha = zeros(size(s,1),size(s,2)); % 68 x 68
beta = zeros(size(s)); % 68 x 68

% Iterate over every element of the alpha/beta matrix
for i = 1:size(s, 1)
    for j = 1:size(s, 2)
        design_matrix = [ones(1,19) ; reshape(s(i,j,:),1,19)];

        % Solve for the alpha and beta coefficients
        coeff = reshape(f(i,j,:),1,19) \ design_matrix;
        alpha(i,j) = coeff(1);
        beta(i,j) = coeff(1);
    end
end