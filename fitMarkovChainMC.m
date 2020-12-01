function [parameters, error, residual, accept_ratio, errors, best_params, ...
    minimum_obj, all_params] = fitMarkovChainMC(training_data, ...
    max_iter, term_iter, final_term_iter, param_logs, param_low_limits, ...
    param_high_limits, scale_factors,T,mw,Ts,rhos, delta_h, params_initial)

% indices is necessary for use of modelError function that is optimized at
% for this exercise, we will use all the data
indices = 1:length(training_data);

% before starting the iterations the initial parameters set to loop 
% parameters and initial objective function is calculated. 
current_params = params_initial;
best_params = params_initial; 
new_params = params_initial;
all_params = zeros(max_iter,3);
[model_error, variance] = modelError(current_params,training_data, mw, ...
    rhos, Ts, indices, delta_h);
current_obj = sum(model_error.^2);
% a couple of counters created to get the total number of acceptance and
% to check the termination criteria. 
iter = 1; 
true_accept = 0;
rand_accept = 0;
unchange_iter = 0;
unchange_iter_term = 0;

norm_factor = variance; % normalization factor is set to the variance of 
% the values.
minimum_obj = current_obj;
errors = zeros(max_iter, 1);

while iter < max_iter
    
    % the random walk stage. If the random walk is out of the parameter 
    % ranges, the random walk process restarts.
    while 1 
        for i=1:length(params_initial)
            if param_logs(i) == 1
                new_params(i) = exp(log(current_params(i)) + randn * ...
                    scale_factors(i) * ((max_iter - iter) / max_iter));
            else
                new_params(i) = current_params(i) + randn * ...
                    scale_factors(i) * ((max_iter - iter) / max_iter);
            end
        end
        if sum(new_params<=param_low_limits) == 0 & sum(new_params>= ...
                param_high_limits) == 0
            break
        end
    end
    % after random walk step, the objective function is calculated. 
    new_obj = sum(modelError(new_params,training_data, mw, rhos, Ts, ...
        indices, delta_h).^2);
    iter = iter + 1;
    % depending on the new objective function the random walk is accepted
    % or rejected. 
    if new_obj < current_obj
        % if the new objective function is smaller than the current
        % objective function, new parameter is accepted. 
        current_params = new_params; 
        current_obj = new_obj; 
        true_accept = true_accept + 1;
        unchange_iter = 0;
    else
        % if the new objective function is greater than the current
        % objective function, depending on the Metropolis criterion random
        % walk could be accepted or rejected. The change in the objective
        % function is normalized with the variance of the output terms,
        % temperature term (which adjusts the acceptance rate) and the
        % iteration count which will slowly reduce the acceptance rate with
        % iterations
        if exp(-1/2*(new_obj - current_obj)/variance /T/ ((max_iter - iter) / max_iter)) > rand
            current_params = new_params; 
            current_obj = new_obj; 
            rand_accept = rand_accept + 1;
            unchange_iter = 0;
        else
            % if rejected no parameter changes but the iteration count
            % increases. 
            unchange_iter = unchange_iter + 1;
        end
    end
    % if the new objective function is lower than the best objective
    % function, update best objective function and the best parameter.
    if new_obj < minimum_obj
        unchange_iter_term = 0;
        best_params = new_params; 
        minimum_obj = new_obj; 
    else
        unchange_iter_term = unchange_iter_term + 1;
    end
    % if the best objective function did not get better in a certain
    % amounts of iterations, it will get back to best parameter. 
    if unchange_iter_term > term_iter
        current_params = best_params;
        unchange_iter_term = 0;
    end
    % if the best objective function did not get better in a certain number
    % of iterations, the algorithm is terminated. 
    if unchange_iter_term > final_term_iter
        break
    end
    
    % after accepting or rejecting, the acceptance ratio is calculated. 
    accept_ratio = (true_accept + rand_accept) / iter; 
    if mod(iter, 100) == 0
        % the variance term is updated at every 100 iterations. 
        [model_error, variance] = modelError(best_params,training_data, ...
            mw, rhos, Ts, indices, delta_h);
    end
end

% the return values are prepared at the end of function. 
parameters = current_params; 
error = current_obj;
residual = modelError(current_params,training_data, mw, rhos, Ts, ...
    indices, delta_h);
accept_ratio = (true_accept + rand_accept) / iter; 

end