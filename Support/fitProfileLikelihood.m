function error = fitProfileLikelihood(training_data, para, ...
    idx, delta_h,  mw, Ts, rhos, params_initial)
%% This function is an alternative form of nonlinear least squares ...
% Instead of estimating all parameters, all but one parameters are refitted.  

% indices are used in modelError function. All estimation methods except CV 
% uses all data
indices = 1:length(training_data);

options = optimoptions('lsqnonlin','Display','iter');
% depending on the index of specified parameter a different nonlinear least
% squares. 

if idx==1
    params_init = params_initial(2:3);
    [param,~,~,~,~] = lsqnonlin(@(params) ...
        modelError([para,params(1),params(2)],training_data, mw, rhos, ...
        Ts,indices, delta_h),params_init,[0 0],[Inf Inf],options);
elseif idx==2
    params_init = [params_initial(1) params_initial(3)];
    [param,~,~,~,~] = lsqnonlin(@(params) ...
        modelError([params(1),para,params(2)],training_data, mw, rhos, ...
        Ts,indices, delta_h),params_init,[0 0],[Inf Inf],options);
else
    params_init = params_initial(1:2);
    [param,~,~,~,~] = lsqnonlin(@(params) ...
        modelError([params(1),params(2),para],training_data, mw, rhos, ...
        Ts,indices, delta_h),params_init,[0 0],[Inf Inf],options);
end

% After estimating parameters, a full set of three parameters is created by
% adding the fixed parameter
if idx == 1
    parax = [para, param(1), param(2)];
elseif idx==2
    parax = [param(1), para, param(2)];
else
    parax = [param(1), param(2), para];
end

% the objective function is calculated with the full final parameter set. 
errors = modelError(parax, training_data, mw, rhos, Ts, indices, delta_h);
error = sum(errors.^2);
end