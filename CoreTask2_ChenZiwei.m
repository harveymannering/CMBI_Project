cd D:\Download\connectomes-data\Task2Data
% Based on Harvey's previous code
% Load csv files
s = zeros(68, 68, 19);
f = zeros(68, 68, 19);
for i = 1:19
    s(:,:,i) = csvread(strcat('D:\Download\connectomes-data\Task2Data\', int2str(31+i), '_WFA_68.csv'), 1, 0);
    f(:,:,i) = csvread(strcat('D:\Download\connectomes-data\Task2Data\', int2str(31+i), '_rsfMRI_68.csv') , 1, 0);
end
% model 1
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
for t = 1:19
    res1 = f(:,:,t)-(alpha1+beta1.*s(:,:,t));
    model1_sum = model1_sum + res1.^2;
end

aic1 = 2*3 + 19.*log(model1_sum.*1/19);
bic1 = 3*log(19) + 19.*log(model1_sum.*1/19);


% The following is my previous code
% %  indirect structural connectivity matrix
% str_data_dir = dir('*_WFA_68.csv');
% f_data_dir = dir('*_rsfMRI_68.csv');
% for t = 1:size(f_data_dir,1) 
%     f_data{t} = readmatrix(f_data_dir(t).name);
% end
% 
% for t = 1:size(str_data_dir,1) 
%     s_matrix = readmatrix(str_data_dir(t).name);
%     str_data{t} = s_matrix;
%     len = length(s_matrix);
%     t_matrix = zeros(len);
%     
%     for i = 1:len
%         for j = 1:len 
%             
%             min_list = [];
%             for k = 1:len
%                 if s_matrix(i,k)~=0 && s_matrix(k,j)~=0
%                     min_list = [min_list, min(s_matrix(i,k),s_matrix(k,j))];
%                 end
%             end
%             
%             if length(min_list)~=0
%                 t_matrix(i,j) = max(min_list);  
%             else
%                 t_matrix(i,j) = 0;
%             end
%         end
%     end
%     ind_str_matrix{t} = t_matrix;  
% end
% 
% % Fit the model
% % model 1
% alpha{1} = zeros(len);
% beta{1} = zeros(len);
% num_measurements{1} = zeros(len);
% 
% for i = 1:len
%     for j = 1:len
%              
%         func = [];
%         structural = [];
%         
%         for t = 1:size(f_data,2)
%             if str_data{t}(i,j) ~= 0
%                 if f_data{t}(i,j) ~=0
%                     func(end+1) = f_data{t}(i,j);
%                     structural(end+1) = str_data{t}(i,j);
%                 end
%             end
%         end
%         
%         if length(func) ~=0
%             veclen = length(structural);
%             design_mat = [ones(veclen,1),structural'];
%             X = pinv(design_mat'*design_mat)*design_mat'*func';
%             alpha{1}(i,j) = X(1);
%             beta{1}(i,j) = X(2);
%             num_measurements{1}(i,j) = veclen;
% 
%         end
%         
%     end
% end
% 
% mod1_sum = zeros(68);
% for t = 1:19
%     f_store = f_data{t};
%     alpha_val = alpha{1};
%     beta_val = beta{1};
% 
%     res = f_store - (alpha_val+beta_val.*str_data{t});
%     res = (res.^2);
%     mod1_sum = mod1_sum+res;
% end
% 
% aic{1} = 2*3 + num_measurements{1} .* log ((1/num_measurements{1}).*mod1_sum);

