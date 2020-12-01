function [error, variance] = modelError(params, data, mw, rhos, Ts, indices, delta_h)
% This function defines the model error used in all parameter estimation
% methods

%% Kinetic parameters

k0_forward = params(1);
k0_backward = params(2);
Ea_forward = params(3);
Ea_backward = Ea_forward + delta_h;

%% Experimental data parsing

selected_data = data(indices,:);

% T, WHSV, H2SDMS, weight, conversion and selectivity data are collected
% from excel file. 

T_list = selected_data(:,1); % inlet temperatures
ratio_list = selected_data(:,2); % inlet feed ratio
ghsv_list = selected_data(:,3); % ghsv
c_experimental_list = selected_data(:,4:7); % concentration data 

%% Model error 

% for each data point model is evaluated and the error is
% lsqnonlin function requires array of errors as input. 

% Returns simulated exit concentration values
for i=1:length(T_list)
    [cA(i,1), cB(i,1), cC(i,1), cD(i,1)] = runModel(T_list(i),ghsv_list(i), ...
        ratio_list(i), k0_forward, k0_backward, Ea_forward, Ea_backward, ...
        mw, rhos, Ts); 
end

% Error initialization
error = zeros(length(selected_data) * 4, 1);

% Error calculation 1x400 array
for i=1:length(T_list)
    error((i-1)*4 + 1) =  (cA(i) - c_experimental_list(i,1));
    error((i-1)*4 + 2) =  (cB(i) - c_experimental_list(i,2));
    error((i-1)*4 + 3) =  (cC(i) - c_experimental_list(i,3));
    error(i*4) =  (cD(i) - c_experimental_list(i,4));
end


%% Please ignore for least squares & cross-val least squares!! 
% (Variance calculation is only needed for the genetic algorithm.)
all_values = [cA(:,1); cB(:,1); cC(:,1); cD(:,1)];
variance = var(all_values);

end

