function [objective_functions] = calculateProfileLikelihood(parameters, ...
    training_data, prefix, number_of_points, mw, Ts, rhos, ...
    params_initial, delta_h)
%% This function generates likelihood profiles of each parameter for a given range

%% Initialization

% the number of plots is set to number of parameters
number_of_plots = length(parameters);

% initialization of the objective functions
objective_functions = zeros(number_of_points,number_of_plots);

%% PL Evaluation

for i=1:number_of_plots
    % first two parameters: pre-exponential factors). The third parameter 
    % is the activation energy 
    if i<3
        minimum_parameter = parameters(i) * 1e-2; 
        maximum_parameter = parameters(i) * 1e2;
    else
        minimum_parameter = parameters(i) * 10^-0.5; 
        maximum_parameter = parameters(i) * 10^0.5;
    end

    % A number of parameters are created within the specified range
    
    params_to_run = logspace(log10(minimum_parameter), ...
        log10(maximum_parameter), number_of_points);
    
    for j=1:length(params_to_run)
        params = parameters;
        params(i) = params_to_run(j);
        % This function solves the parameter estimation problem using
        % specified parameter using nonlinear least squares method. The
        % difference between original nonlinear least square application
        % and this application is the specification of one parameter. As a
        % result, this function solves the same optimization problem
        % without traning for specified parameter. 
        error = fitProfileLikelihood(training_data, ...
            params(i), i, delta_h,  mw, Ts, rhos, params_initial);
        objective_functions(j,i) = error;
    end
    
    % The minimum objective functions (profile likelihood) are plotted and 
    % saved for each parameter. The prefix term gives information about the
    % identifiers for the data. Therefore, the profile likelihood plots
    % have described specific naming. 
    figure
    plot(log10(params_to_run), objective_functions(:,i));
    set(gca,'FontSize',16)
    xlabel(strcat('log_{10} \theta_{',num2str(i),'}'));
    
    set(gca,'FontSize',16)

    savefig(strcat(prefix, '_pl_', num2str(i),'.fig'));
    save(strcat(prefix, '_pl_', num2str(i)));
    
    set(gca,'FontSize',16)

    
end

end
