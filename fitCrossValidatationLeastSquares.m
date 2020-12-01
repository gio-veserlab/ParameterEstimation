function [parameters, error, residual, weight, mean_obj, ...
    params] = fitCrossValidatationLeastSquares(training_data, delta_h, fold, ...
    mw, Ts, rhos, params_initial)
%%fold generation

% stratified sampling
groups = training_data(:,1);

% syntax: crossvalind('methodname',length of the training data or stratification,number of folds)
% random sampliing: indices = crossvalind('Kfold',100,fold);

indices = crossvalind('Kfold',groups,fold);

% Sometimes folds created might be empty if there is a small number of datapoints
% (e.g. 10 datapoints with 5 folds).This loop makes sure that the total number 
% of folds match the non-empty number of created folds. 
while length(unique(indices)) < fold
    indices = crossvalind('Kfold',groups,fold);
end

%% initialization

% each fold generates a set of kinetic parameters and stored in "params" matrix 
params = zeros(length(params_initial),fold);

% mean objective function of each fold on its respective test data
mean_obj = zeros(1,fold);

%% training

% For each fold parameter estimation is performed using nonlinear least
% squares method. After the parameter set is derived, error is calculated on
% the hold-out fold.

for i=1:fold
    % separate experimental data into test and training set for each fold.
    test = find(indices == i); % test set
    train = find(indices ~= i); % training set
    
    % non-linear least squares solution for training set
    options = optimoptions('lsqnonlin','Display','iter', ...
        'MaxFunctionEvaluations',3000,'OptimalityTolerance',1e-16, ...
        'FunctionTolerance',1e-16,'StepTolerance',1e-6);
    
    para = lsqnonlin(@(params) ... 
        modelError(params,training_data, mw, rhos, Ts, train, delta_h),...
        params_initial,[0 0 0],[Inf Inf Inf],options);
    
    % calculation of mean objective function on test data for each fold
    params(:,i) = para(:);
    % sum of the square error is used for weighting
    mean_obj(i) = mean(modelError(params(:,i),training_data, mw, rhos, Ts, ...
        test, delta_h).^2);
end

%% averaging all folds and finalization of kinetic parameters

% initialization of the final parameter set and weights for each hold. 
params_normalized = zeros(length(params_initial),1);
weight = zeros(fold,1);

for i=1:fold
    % objective function for each fold inversely proportional to the weight
    weight(i) = 1/mean_obj(i)^2;
    % final parameters are the weighted sum of the parameters estimated for
    % each fold.
    params_normalized = params_normalized + params(:,i) * weight(i);
end
% normalization
parameters = params_normalized / sum(weight);


residual = modelError(parameters,training_data, mw, rhos, Ts,1:length(training_data), delta_h);
error = sum(residual.^2);

end