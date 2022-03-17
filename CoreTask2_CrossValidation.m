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

% Cross validation loop
for c = 1:19
    % Split data into training/validation set
    training_set = s(:,:,1:size(s,3)~=c);
    training_set_targets = f(:,:,1:size(f,3)~=c);
    validation_set = s(:,:,c);
    validation_set_targets = f(:,:,c);


    alpha1 = zeros(size(s,1), size(s,2)); % 68 x 68
    beta1 = zeros(size(s,1), size(s,2)); % 68 x 68
    
    % Iterate over every element of the alpha/beta matrix
    for i = 1:size(s, 1)
        for j = 1:size(s, 2)
            % Put variables into matrix form: Y = coefficients * X
            X = [ones(18, 1) reshape(training_set(i,j,:), 18, 1)];
            Y = reshape(training_set_targets(i,j,:), 18, 1);
    
            % Solve for the alpha and beta coefficients
            coeff = pinv(X' * X) * X' * Y;
            alpha1(i,j) = coeff(1);
            beta1(i,j) = coeff(2);
    
        end
    end

    % Calculate sum of square errors
    model1_sum = zeros(68);
    res1 = validation_set_targets(:,:)-(alpha1+beta1.*validation_set(:,:));
    model1_sum = model1_sum + res1.^2;
    SSE = sum(sum(model1_sum));
    disp(SSE);
    
    % Calculate AIC/BIC
    aic1 = 2*3 + 19.*log(model1_sum.*1/19);
    bic1 = 3*log(19) + 19.*log(model1_sum.*1/19);

end

%% Fit Model 2: f = alpha + beta * s + y * s^2
% Cross validation loop
for c = 1:19
    % Split data into training/validation set
    training_set = s(:,:,1:size(s,3)~=c);
    training_set_targets = f(:,:,1:size(f,3)~=c);
    validation_set = s(:,:,c);
    validation_set_targets = f(:,:,c);


    alpha2 = zeros(size(s,1), size(s,2)); % 68 x 68
    beta2 = zeros(size(s,1), size(s,2)); % 68 x 68
    y2 = zeros(size(s,1), size(s,2)); % 68 x 68
    
    % Iterate over every element of the alpha/beta/y matrix
    for i = 1:size(s, 1)
        for j = 1:size(s, 2)
            % Put variables into matrix form: Y = coefficients * X
            s_slice = reshape(training_set(i,j,:), 18, 1);
            X = [ones(18, 1) s_slice (s_slice.^2)];
            Y = reshape(training_set_targets(i,j,:), 18, 1);
    
            % Solve for the alpha and beta coefficients
            coeff = pinv(X' * X) * X' * Y;
            alpha2(i,j) = coeff(1);
            beta2(i,j) = coeff(2);
            y2(i,j) = coeff(3);
        end
    end
    
    model2_sum = zeros(68);
    res2 = validation_set_targets(:,:) - (alpha2 + beta2.*validation_set(:,:) + y2.*validation_set(:,:).^2);
    model2_sum = model2_sum + res2.^2;
    SSE = sum(sum(model2_sum));
    %disp(SSE);

    aic2 = 2*4 + 19.*log(model2_sum./19);
    bic2 = 4*log(19) + 19.*log(model2_sum./19);

end

%% Fit Model 3: f = alpha + beta * t
% Cross validation loop
for c = 1:19
    % Split data into training/validation set
    training_set = t(:,:,1:size(t,3)~=c);
    training_set_targets = f(:,:,1:size(f,3)~=c);
    validation_set = t(:,:,c);
    validation_set_targets = f(:,:,c);

    alpha3 = zeros(size(s,1),size(s,2)); % 68 x 68
    beta3 = zeros(size(s,1),size(s,2)); % 68 x 68
    
    % Iterate over every element of the alpha/beta matrix
    for i = 1:size(s, 1)
        for j = 1:size(s, 2)
            % Put variables into matrix form: Y = coefficients * X
            X = [ones(18, 1) reshape(training_set(i,j,:), 18, 1)];
            Y = reshape(training_set_targets(i,j,:), 18, 1);
    
            % Solve for the alpha and beta coefficients
            coeff = pinv(X' * X) * X' * Y;
            alpha3(i,j) = coeff(1);
            beta3(i,j) = coeff(2);
        end
    end
    
    % There seems to be some NaN values in the alpha3 and beta3 variables,
    % so these need to be set to zero
    alpha3(isnan(alpha3)) = 0;
    beta3(isnan(beta3)) = 0;

    % Calculate sum of squares
    model3_sum = zeros(68);
    res3 = validation_set_targets(:,:) - (alpha3 + beta3.*validation_set(:,:));
    model3_sum = model3_sum + res3.^2;
    model3_sum(isinf(model3_sum)) = 0; % Set inifity values to zeros 
    SSE = sum(sum(model3_sum));
    %disp(SSE);
    
    aic3 = 2*3 + 19.* log(model3_sum./19);
    bic3 = 3*log(19) + 19.*log(model3_sum./19);
end

%% Fit Model 4: f = alpha + beta * t + y * t^2
% Cross validation loop
for c = 1:19
    % Split data into training/validation set
    training_set = t(:,:,1:size(t,3)~=c);
    training_set_targets = f(:,:,1:size(f,3)~=c);
    validation_set = t(:,:,c);
    validation_set_targets = f(:,:,c);

    alpha4 = zeros(size(s,1),size(s,2)); % 68 x 68
    beta4 = zeros(size(s,1),size(s,2)); % 68 x 68
    y4 = zeros(size(s,1),size(s,2)); % 68 x 68
    
    % Iterate over every element of the alpha/beta/y matrix
    for i = 1:size(s, 1)
        for j = 1:size(s, 2)
            % Put variables into matrix form: Y = coefficients * X
            t_slice = reshape(training_set(i,j,:), 18, 1);
            X = [ones(18, 1) t_slice (t_slice.^2)];
            Y = reshape(training_set_targets(i,j,:), 18, 1);
    
            % Solve for the alpha and beta coefficients
            coeff = pinv(X' * X) * X' * Y;
            alpha4(i,j) = coeff(1);
            beta4(i,j) = coeff(2);
            y4(i,j) = coeff(3);
        end
    end

    % There seems to be some NaN values in the alpha3 and beta3 variables,
    % so these need to be set to zero
    alpha4(isnan(alpha4)) = 0;
    beta4(isnan(beta4)) = 0;
    y4(isnan(y4)) = 0;
    
    % Calculate sum of square errors
    model4_sum = zeros(68);
    res4 = validation_set_targets(:,:) - (alpha4 + beta4.*validation_set(:,:) + y4.*validation_set(:,:).^2);
    model4_sum = model4_sum + res4.^2;
    model4_sum(isnan(model4_sum)) = 0; % Set NaN values to zeros again 
    model4_sum(isinf(model4_sum)) = 0; % Set inifity values to zeros 
    SSE = sum(sum(model4_sum));
    %disp(SSE);

    % AIC/BIC calculation
    aic4 = 2*4 + 19.* log(model4_sum./19);
    bic4 = 4*log(19) + 19.*log(model4_sum./19);
   
end


%% Fit Model 5: f = alpha + beta * s + y * t
% Cross validation loop
for c = 1:19
    % Split data into training/validation set
    training_set_s = s(:,:,1:size(t,3)~=c);
    training_set_t = t(:,:,1:size(t,3)~=c);
    training_set_targets = f(:,:,1:size(f,3)~=c);
    validation_set_s = s(:,:,c);
    validation_set_t = t(:,:,c);
    validation_set_targets = f(:,:,c);

    alpha5 = zeros(size(s,1),size(s,2)); % 68 x 68
    beta5 = zeros(size(s,1),size(s,2)); % 68 x 68
    y5 = zeros(size(s,1),size(s,2)); % 68 x 68
    
    % Iterate over every element of the alpha/beta/y matrix
    for i = 1:size(s, 1)
        for j = 1:size(s, 2)
            % Put variables into matrix form: Y = coefficients * X
            X = [ones(18, 1) reshape(training_set_s(i,j,:), 18, 1) reshape(training_set_t(i,j,:), 18, 1)];
            Y = reshape(training_set_targets(i,j,:), 18, 1);
    
            % Solve for the alpha and beta coefficients
            coeff = pinv(X' * X) * X' * Y;
            alpha5(i,j) = coeff(1);
            beta5(i,j) = coeff(2);
            y5(i,j) = coeff(3);
        end
    end
        
    % There seems to be some NaN values in the alpha3 and beta3 variables,
    % so these need to be set to zero
    alpha5(isnan(alpha5)) = 0;
    beta5(isnan(beta5)) = 0;
    y5(isnan(y5)) = 0;
    
    % Calculate sum of square errors
    model5_sum = zeros(68);
    res5 = validation_set_targets(:,:) - (alpha5 + beta5.*validation_set_s(:,:) + y5.*validation_set_t(:,:));
    model5_sum = model5_sum + res5.^2;
    model5_sum(isnan(model5_sum)) = 0; % Set NaN values to zeros again 
    model5_sum(isinf(model5_sum)) = 0; % Set inifity values to zeros 
    SSE = sum(sum(model5_sum));
    disp(SSE);
    
    
    aic5 = 2*4 + 19.* log(model5_sum./19);
    bic5 = 4*log(19) + 19.*log(model5_sum./19);
end
