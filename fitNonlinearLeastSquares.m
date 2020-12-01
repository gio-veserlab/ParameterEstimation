function [parameters, error, residual] = ...
    fitNonlinearLeastSquares(training_data, delta_h, mw, Ts, rhos, ...
    params_initial)

% indices is used in modelError function.(Indices change only for CV
% method. All other methods use all available data)

indices = 1:length(training_data);

% options for lsqnonlin function
% options = optimoptions('lsqnonlin','Display','iter', ... 
%    'MaxFunctionEvaluations',3000,'OptimalityTolerance',1e-16, ...
%    'FunctionTolerance',1e-16,'StepTolerance',1e-6);

% lsqnonlin used to estimate parameters 
[para,resnorm,residual,~,~] = lsqnonlin(@(params) ...
    modelError(params,training_data, mw, rhos, Ts,indices, ...
    delta_h),params_initial,[0 0 0],[Inf Inf Inf]);

% the return values are prepared
parameters = para; 
error = resnorm; 
end