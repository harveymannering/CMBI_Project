%% Load data
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

%% Calculate t matrix
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
                if (min_weight > greatest_min_weight && min_weight ~= 0)
                    greatest_min_weight = min_weight;
                end
            end
    
            % Set matrix element of the indirect structural connectivity matrix
            t(t_i, t_j, p) = greatest_min_weight;
        end 
    end 
end 

%% Fit Model 1: f = alpha + beta * s
alpha1 = zeros(size(s,1),size(s,2)); % 68 x 68
beta1 = zeros(size(s,1),size(s,2)); % 68 x 68

% Iterate over every element of the alpha/beta matrix
for i = 1:size(s, 1)
    for j = 1:size(s, 2)
        % Put variables into matrix form: Y = coefficients * X
        X = [ones(19, 1) reshape(s(i,j,:), 19, 1)];
        Y = reshape(f(i,j,:), 19, 1);

        % Solve for the alpha and beta coefficients
        coeff = pinv(X' * X) * X' * Y;
        alpha1(i,j) = coeff(1);
        beta1(i,j) = coeff(2);

    end
end

%% Fit Model 2: f = alpha + beta * s + y * s^2
alpha2 = zeros(size(s,1),size(s,2)); % 68 x 68
beta2 = zeros(size(s,1),size(s,2)); % 68 x 68
y2 = zeros(size(s,1),size(s,2)); % 68 x 68

% Iterate over every element of the alpha/beta/y matrix
for i = 1:size(s, 1)
    for j = 1:size(s, 2)
        % Put variables into matrix form: Y = coefficients * X
        s_slice = reshape(s(i,j,:), 19, 1);
        X = [ones(19, 1) s_slice (s_slice.^2)];
        Y = reshape(f(i,j,:), 19, 1);

        % Solve for the alpha and beta coefficients
        coeff = pinv(X' * X) * X' * Y;
        alpha2(i,j) = coeff(1);
        beta2(i,j) = coeff(2);
        y2(i,j) = coeff(3);
    end
end


%% Fit Model 3: f = alpha + beta * t
alpha3 = zeros(size(s,1),size(s,2)); % 68 x 68
beta3 = zeros(size(s,1),size(s,2)); % 68 x 68

% Iterate over every element of the alpha/beta matrix
for i = 1:size(s, 1)
    for j = 1:size(s, 2)
        % Put variables into matrix form: Y = coefficients * X
        X = [ones(19, 1) reshape(t(i,j,:), 19, 1)];
        Y = reshape(f(i,j,:), 19, 1);

        % Solve for the alpha and beta coefficients
        coeff = pinv(X' * X) * X' * Y;
        alpha3(i,j) = coeff(1);
        beta3(i,j) = coeff(2);
    end
end

%% Fit Model 4: f = alpha + beta * t + y * t^2
alpha4 = zeros(size(s,1),size(s,2)); % 68 x 68
beta4 = zeros(size(s,1),size(s,2)); % 68 x 68
y4 = zeros(size(s,1),size(s,2)); % 68 x 68

% Iterate over every element of the alpha/beta/y matrix
for i = 1:size(s, 1)
    for j = 1:size(s, 2)
        % Put variables into matrix form: Y = coefficients * X
        t_slice = reshape(t(i,j,:), 19, 1);
        X = [ones(19, 1) t_slice (t_slice.^2)];
        Y = reshape(f(i,j,:), 19, 1);

        % Solve for the alpha and beta coefficients
        coeff = pinv(X' * X) * X' * Y;
        alpha4(i,j) = coeff(1);
        beta4(i,j) = coeff(2);
        y4(i,j) = coeff(3);
    end
end

%% Fit Model 5: f = alpha + beta * s + y * t
alpha5 = zeros(size(s,1),size(s,2)); % 68 x 68
beta5 = zeros(size(s,1),size(s,2)); % 68 x 68
y5 = zeros(size(s,1),size(s,2)); % 68 x 68

% Iterate over every element of the alpha/beta/y matrix
for i = 1:size(s, 1)
    for j = 1:size(s, 2)
        % Put variables into matrix form: Y = coefficients * X
        X = [ones(19, 1) reshape(s(i,j,:), 19, 1) reshape(t(i,j,:), 19, 1)];
        Y = reshape(f(i,j,:), 19, 1);

        % Solve for the alpha and beta coefficients
        coeff = pinv(X' * X) * X' * Y;
        alpha5(i,j) = coeff(1);
        beta5(i,j) = coeff(2);
        y5(i,j) = coeff(3);
    end
end