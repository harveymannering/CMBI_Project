clear all;
cd D:\Download\connectomes-data\Task2Data

% Based on Harvey's previous code

% Load csv files
s = zeros(68, 68, 19);
f = zeros(68, 68, 19);
for i = 1:19
    s(:,:,i) = csvread(strcat('D:\Download\connectomes-data\Task2Data\', int2str(31+i), '_WFA_68.csv'), 1, 0);
    f(:,:,i) = csvread(strcat('D:\Download\connectomes-data\Task2Data\', int2str(31+i), '_rsfMRI_68.csv') , 1, 0);
end

% Calucluate indirect structural connectivity matrix
t = zeros(size(s));
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

%% model 1
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

model1_sum = zeros(68);
for k = 1:19
    res1 = f(:,:,k)-(alpha1+beta1.*s(:,:,k));
    model1_sum = model1_sum + res1.^2;
end

aic1 = 2*3 + 19.*log(model1_sum.*1/19);
bic1 = 3*log(19) + 19.*log(model1_sum.*1/19);

cv_sum1 = 0;
for cv = 1:19
    f_test = f(:,:,cv);
    f_train = f(:,:,[1:cv-1,cv+1:19]);
    s_test = s(:,:,cv);
    s_train = s(:,:,[1:cv-1,cv+1:19]);
    
    cv_alpha1 = zeros(size(s,1),size(s,2));
    cv_beta1 = zeros(size(s,1),size(s,2));
    
    for m = 1:1:size(s, 1)
        for n = 1:size(s, 2)
            cv_X1 = [ones(18, 1) reshape(s_train(m,n,:), 18, 1)];
            cv_Y1 = reshape(f_train(m,n,:), 18, 1);
            cv_coeff1 = pinv(cv_X1' * cv_X1) * cv_X1' * cv_Y1;
            cv_alpha1(i,j) = cv_coeff1(1);
            cv_beta1(i,j) = cv_coeff1(2);
        end
    end
    cv_res1 = f_test - (cv_alpha1 + cv_beta1.*s_test);
    cv_sum1 = cv_sum1 + cv_res1.^2;
end


%% model 2
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

model2_sum = zeros(68);
for k = 1:19
    res2 = f(:,:,k) - (alpha2 + beta2.*s(:,:,k) + y2.*s(:,:,k).^2);
    model2_sum = model2_sum + res2.^2;
end

aic2 = 2*4 + 19.*log(model2_sum./19);
bic2 = 4*log(19) + 19.*log(model2_sum./19);

cv_sum2 = 0;
for cv = 1:19
    f_test = f(:,:,cv);
    f_train = f(:,:,[1:cv-1,cv+1:19]);
    s_test = s(:,:,cv);
    s_train = s(:,:,[1:cv-1,cv+1:19]);
    
    cv_alpha2 = zeros(size(s,1),size(s,2)); 
    cv_beta2 = zeros(size(s,1),size(s,2)); 
    cv_y2 = zeros(size(s,1),size(s,2)); 
    
    for m = 1:1:size(s, 1)
        for n = 1:size(s, 2)
            cv_s_slice2 = reshape(s_train(m,n,:), 18, 1);
            cv_X2 = [ones(18, 1) cv_s_slice2 (cv_s_slice2.^2)];
            cv_Y2 = reshape(f_train(m,n,:), 18, 1);
            cv_coeff2 = pinv(cv_X2' * cv_X2) * cv_X2' * cv_Y2;
            cv_alpha2(m,n) = cv_coeff2(1);
            cv_beta2(m,n) = cv_coeff2(2);
            cv_y2(m,n) = cv_coeff2(3);
        end
    end
    cv_res2 = f_test - (cv_alpha2 + cv_beta2.*s_test + cv_y2.*s_test.^2);
    cv_sum2 = cv_sum2 + cv_res2.^2;
end

%% model 3
alpha3 = zeros(size(s,1),size(s,2)); % 68 x 68
beta3 = zeros(size(s,1),size(s,2)); % 68 x 68

% Iterate over every element of the alpha/beta matrix
for i = 1:size(s, 1)
    for j = 1:size(s, 2)
        % Put variables into matrix form: Y = coefficients * X
        X = [ones(19, 1) reshape(t(i,j,:), 19, 1)];
        Y = reshape(f(i,j,:), 19, 1);

        % Solve for the alpha and beta coefficients
        coeff = inv(X' * X) * X' * Y;
        alpha3(i,j) = coeff(1);
        beta3(i,j) = coeff(2);
    end
end

model3_sum = zeros(68);
for k = 1:19
    res3 = f(:,:,k) - (alpha3 + beta3.*t(:,:,k));
    model3_sum = model3_sum + res3.^2;
end

aic3 = 2*3 + 19.* log(model3_sum./19);
bic3 = 3*log(19) + 19.*log(model3_sum./19);

%% model 4
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
        coeff = inv(X' * X) * X' * Y;
        alpha4(i,j) = coeff(1);
        beta4(i,j) = coeff(2);
        y4(i,j) = coeff(3);
    end
end

model4_sum = zeros(68);
for k = 1:19
    res4 = f(:,:,k) - (alpha4 + beta4.*t(:,:,k) + y4.*t(:,:,k).^2);
    model4_sum = model4_sum + res4.^2;
end

aic4 = 2*4 + 19.* log(model4_sum./19);
bic4 = 4*log(19) + 19.*log(model4_sum./19);

%% model 5
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
        coeff = inv(X' * X) * X' * Y;
        alpha5(i,j) = coeff(1);
        beta5(i,j) = coeff(2);
        y5(i,j) = coeff(3);
    end
end

model5_sum = zeros(68);
for k = 1:19
    res5 = f(:,:,k) - (alpha5 + beta5.*s(:,:,k) + y5.*t(:,:,k));
    model5_sum = model5_sum + res5.^2;
end

aic5 = 2*4 + 19.* log(model5_sum./19);
bic5 = 4*log(19) + 19.*log(model5_sum./19);
