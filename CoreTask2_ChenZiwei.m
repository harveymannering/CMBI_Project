
cd D:\Download\connectomes-data\Task2Data
%  indirect structural connectivity matrix
str_data_dir = dir('*_WFA_68.csv');
f_data_dir = dir('*_rsfMRI_68.csv');
for t = 1:size(f_data_dir,1) 
    f_data{t} = readmatrix(f_data_dir(t).name);
end

for t = 1:size(str_data_dir,1) 
    s_matrix = readmatrix(str_data_dir(t).name);
    str_data{t} = s_matrix;
    len = length(s_matrix);
    t_matrix = zeros(len);
    
    for i = 1:len
        for j = 1:len 
            
            min_list = [];
            for k = 1:len
                if s_matrix(i,k)~=0 && s_matrix(k,j)~=0
                    min_list = [min_list, min(s_matrix(i,k),s_matrix(k,j))];
                end
            end
            
            if length(min_list)~=0
                t_matrix(i,j) = max(min_list);  
            else
                t_matrix(i,j) = 0;
            end
        end
    end
    ind_str_matrix{t} = t_matrix;  
end

% Fit the model
% model 1
alpha{1} = zeros(len);
beta{1} = zeros(len);
num_measurements{1} = zeros(len);

for i = 1:len
    for j = 1:len
             
        func = [];
        structural = [];
        
        for t = 1:size(f_data,2)
            if str_data{t}(i,j) ~= 0
                if f_data{t}(i,j) ~=0
                    func(end+1) = f_data{t}(i,j);
                    structural(end+1) = str_data{t}(i,j);
                end
            end
        end
        
        if length(func) ~=0
            veclen = length(structural);
            design_mat = [ones(veclen,1),structural'];
            X = pinv(design_mat'*design_mat)*design_mat'*func';
            alpha{1}(i,j) = X(1);
            beta{1}(i,j) = X(2);
            num_measurements{1}(i,j) = veclen;

        end
        
    end
end

mod1_sum = zeros(68);
for t = 1:19
    f_store = f_data{t};
    alpha_val = alpha{1};
    beta_val = beta{1};

    res = f_store - (alpha_val+beta_val.*str_data{t});
    res = (res.^2);
    mod1_sum = mod1_sum+res;
end

aic{1} = 2*3 - num_measurements{1} .* log (num_measurements{1}.*mod1_sum);

