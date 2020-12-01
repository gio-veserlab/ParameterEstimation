function [parameters, error, residual, generation] = fitGeneticAlgorithm(training_data, ...
    max_generation, term_generation, elitism, ...
    mutation_rate, population_size, param_low_limits, ...
    param_high_limits, param_logs, scale_factors, mw, Ts, rhos, delta_h, ...
    params_initial)

% indices is necessary for use of modelError function that is optimized at
% for this exercise, we will use all the data
indices = 1:length(training_data);

% initial population is generated using Latin hypercube sampling. 
% initialization of the population
population = zeros(population_size, length(params_initial));
% the lhsdesign function generates numbers between zero and one
X = lhsdesign(population_size,length(params_initial));
% pre-exponential factors are determined using the logspace and activation
% energy is calculated using linear scaling. 
for idx=1:population_size
    population(idx,1:2) = param_low_limits(1:2) .* ...
        (10.^(log10(param_high_limits(1:2)./param_low_limits(1:2)) ...
        .* X(idx,1:2)));
    population(idx,3) = param_low_limits(3) + (param_high_limits(3) - ...
        param_low_limits(3)) * X(idx,3);
end

% total number of parameter sets which selected as elite group.  
elite_size = floor(population_size * elitism);

% the initial generation is set to zero and the objective function array is
% created. 
generation = 0;
objective_functions = zeros(population_size,1);

% Initially, the best objective function is set to infinity and the best
% parameter set to zero. gene_term is the counter that denotes the last
% time the objective function improved. 
best_objective = Inf;
best_parameter = params_initial;
gene_term = 0; 

while generation < max_generation
    
    % if the best objective function has not been changed for a certain
    % number of generations, iteration is terminated, otherwi
    if gene_term == term_generation
        break
    end
    
    generation = generation + 1;
    fprintf('The current generation is %d\n', generation)
    % The objective function of individual parameter sets is calculated.  
    for i=1:length(population)
        objective_functions(i,1) = mean(modelError(population(i,:), ...
            training_data, mw, rhos, Ts, indices, delta_h).^2);
    end
    
    fprintf('Sorting objective functions\n')
    % objective functions are sorted in ascending order and the top elite
    % group of the population selected as elite population. 
    [~,sorted_indices] = sort(objective_functions,'ascend');
    elite_population = population(sorted_indices(1:elite_size),:);
    
    % the new population that is created as offsprings of elite population
    % is initialized. 
    generated_population = zeros(population_size-elite_size,length(params_initial));
    
    % the best objective function and corresponding parameter set is
    % selected. 
    new_best_objective = objective_functions(sorted_indices(1));
    new_best_parameter = population(sorted_indices(1),:);
    
    % two parents are randomly selected from the elite population and the
    % crossover point is also detected. The parameters from the first
    % parameter set until the crossover point and the parameters from the
    % second parameter set after the crossover point are concatenated and
    % the final parameter is produced. The same concatenation is performed
    % on the remaining parameters. The loop guarantees that all elite
    % parameter sets used at least once in the generation of new
    % population. 
    fprintf('Crossing elites for new population\n')
    while 1
        use_rate = zeros(1,elite_size);
        for i=1:length(generated_population)/2
           first = randi(elite_size);
           second = randi(elite_size); 
           while second == first
                second = randi(elite_size); 
           end
           use_rate(first) = use_rate(first) + 1;
           use_rate(second) = use_rate(second) + 1;
           %cross = randi(length(params_initial)-1);
           r1 = rand; 
           r2 = rand; 
           generated_population(2*i-1,1:2) = exp(log(elite_population(first,1:2)) * r1 +  log(elite_population(second,1:2)) * (1-r1));
           generated_population(2*i-1,3) = elite_population(first,3) * r1 + elite_population(second,3) * (1-r1);
           generated_population(2*i,1:2) = exp(log(elite_population(first,1:2)) * r2 +  log(elite_population(second,1:2)) * (1-r2));
           generated_population(2*i,3) = elite_population(first,3) * r2 + elite_population(second,3) * (1-r2);
        end
        if sum(use_rate == 0) > elite_size
            continue
        else
            break
        end
    end
    
    % the generated population and the elite population are concatenated. 
    premutated_population = [elite_population; generated_population];
    new_population = premutated_population;
    
    % mutation rows and columns are generated randomly and mutation
    % operations are applied to the data. 
    total_mut_number = floor(population_size * length(params_initial) * mutation_rate);
    mutation_rows = randi(population_size,total_mut_number,1);
    mutation_cols = randi(length(params_initial),total_mut_number,1);
    
    % Mutation is applied in logscale if the parameter is a preexponential
    % factor and applied in linear scale if the parameter is activation
    % energy. 
    fprintf('Mutating population\n')
    for i=1:length(mutation_rows)
        % check to see if the parameter scales in log or linear scale
       if ismember(mutation_cols(i), param_logs)
           new_population(mutation_rows(i),mutation_cols(i)) = ...
               exp(log(new_population(mutation_rows(i),mutation_cols(i))) + ...
               randn * scale_factors(mutation_cols(i)));
       else
           new_population(mutation_rows(i),mutation_cols(i)) = ...
               new_population(mutation_rows(i),mutation_cols(i)) + ...
               randn * scale_factors(mutation_cols(i));
       end
    end
    
    population = new_population;
    
    % if the new best objective function better than the previous best
    % objective function, the best objective function is updated. The 
    % parameter that tracks the last time best objective function improved
    % is set to zero. Otherwise, the term incremented by 1. 
    if new_best_objective < best_objective
        best_objective = new_best_objective;
        best_parameter = new_best_parameter; 
        gene_term = 0;
    else
        gene_term = gene_term + 1;
    end
end 

% output parameters are prepared for return. 
parameters = best_parameter; 
error = best_objective; 
residual = modelError(parameters,training_data, mw, rhos, Ts, indices, delta_h);

end
