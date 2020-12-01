%% This script run the full parameter estimation and model evaluation pipeline for the reaction
% system defined as A + B <--> C + D 

% Model assumes Arrhenius type kinetics, only mass balance is considered (ambient P, isothermal rxn)

% The parameter estimation workflow needs an experimental concentration data that has
% three input conditions for each data point: temperature, inlet molar
% ratio of A/B and the gas hourly space velocity (residence time) 
%%
% DO NOT RUN THIS SCRIPT AS WHOLE. RUN EACH SECTION INDIVIDUALLY %%%

addpath('Support/');
%%  Property definitions 

% molecular weights of each component
mw = [28;18;44;2]; %g/mol A, B, C, D

% density values of each component for a given temperature range (normal pressure)
physical_prop = xlsread('physicalproperties.xlsx');

% temperature values are given in first column of the xls file
Ts = physical_prop(:,1); % temperatures (C)

% density values are given from the third to sixth column (4 components)
rhos = physical_prop(:,3:6); % vapor densities (kg/m3)   

% Reaction enthalpy

delta_h = 40000; % J/mol

%% Read experimental data

% This function loads the data provided in an excel spreadsheet. 
training_data = xlsread('data.xlsx');

%% train with non-linear least squares (LSQNONLIN)

% initial guess - must be specified for each parameter estimation problem! 
params_initial = [3e5 2.5e8 50000]; % k0_fwd k0_bwd Ea_fwd

% parameter estimation by nonlinear least squares
[lsq_parameters, lsq_training_error, lsq_training_resi] = ...
    fitNonlinearLeastSquares(training_data, delta_h, mw, Ts, rhos, ...
    params_initial);

%% train with cross-validation nonlinear least squares

% Number of folds - must be defined by the user
fold=10;

params_initial = [3e5 2.5e8 50000]; % k0_fwd k0_bwd Ea_fwd

% parameter estimation by cross-validation nonlinear least squares
[cvlsq_parameters, cvlsq_error, cvresiduals, ...
    cv_weights, cvlsq_errors, cv_ind_params] = ...
    fitCrossValidatationLeastSquares(training_data, delta_h, fold, mw, Ts, ...
    rhos, params_initial);

%% train with Markov chain Monte Carlo (MCMC)

% Hyperparameters (algorithm specific parameters) are defined for MCMC 

max_iter = 25000; % maximum number of iterations before termination
term_iter = 200; % maximum number of iterations that will return walker back to
% best solution
final_term_iter = 200; % maximum number of iterations that will lead to termination
% if the best objective function is not improved
param_logs = [1,1,0]; % The logscale for random walks defined here. For 
% preexponential factors a logscale is used, for the activation energy a 
% linear scale is used for random walks. 
scale_factors = [1,1,5000]; % for linear scaling, a factor needs to be defined. 
% The higher and lower limits for each parameters need to be specified. 
param_low_limits = [100,1000,10000];
param_high_limits = [1e12, 1e12, 100000];
% temperature is set here to adjust the acceptance ratios of random walks. 
T = 5;
% initial parameter needs to be defined. 
params_initial = [3e5 2.5e8 50000];

% parameter estimation by MCMC
[mcmc_parameters, mcmc_error, mcmcresidual, accept_ratio] = ...
    fitMarkovChainMC(training_data, max_iter, term_iter, ...
    final_term_iter, param_logs, param_low_limits, param_high_limits, ...
    scale_factors, T, mw, Ts, rhos, delta_h, params_initial);

%%  train with genetic algorithm (GA)

% Hyperparameters (algorithm specific parameters) are defined for GA 

max_generation = 200; % total number of generations before the termination
% of the genetic algorithm 
term_generation = 20; % if the best offspring did not improve for this 
% many generation algorithm will terminate.
elitism = 0.1; % elitism parameter represents the portion of the offsprings
% which generate the next generations (survives)
mutation_rate = 0.1; % after population is generated a portion of 
% individual parameters are changed using a random walk, this simulates the
% mutation in real population
population_size = 400; % total number of parameter sets in each generation
param_logs = [1,1,0]; % The logscale for random walks defined here. For 
% preexponential factors a logscale is used, for the activation energy a 
% linear scale is used for random walks. 
scale_factors = [1,1,5000]; % for linear scaling, a factor needs to be defined. 
% The higher and lower limits for each parameters need to be specified.

% The higher and lower limits for each parameters need to be specified. 
param_low_limits = [100, 1000,10000];
param_high_limits = [1e9, 1e12, 100000];

% parameter estimation by GA
[ga_parameters, ga_error, garesidual, generation] = fitGeneticAlgorithm(training_data, ...
    max_generation, term_generation, elitism, ...
    mutation_rate, population_size, param_low_limits, ...
    param_high_limits, param_logs, scale_factors, mw, Ts, rhos, delta_h, ...
    params_initial);
  
%% PROFILE LIKELIHOOD ESTIMATION 

% the parameters for profile likelihood calculations are defined: 
% The profile likelihood algorithm will create a discrete set of parameters
% between minumum and maximum value. The size of this discrete set needs 
% to be specified. 

number_of_points = 11;

% this function calculates the profile likelihood of each parameter in a
% specified range and creates plots for analysis. 

% IMPORTANT! PL CURRENTLY USES FITTED PARAMETERS OF NONLINEAR LSQ MODEL.
% change from lsq_parameters to cvlsq_parameters if needed.

objective_functions = calculateProfileLikelihood(lsq_parameters, training_data, ...
    'data', number_of_points, mw, Ts, rhos, ...
    params_initial, delta_h);

%% RELATIVE SENSITIVITY

% Relative sensitivity is calculated at specific input conditions.
% Must be determined based on the operating variable range

T = 250; % temperature (C)
ghsv = 3000; % specified gas hour space velocity (hr-1)
ratio = 1; % specified input ratio 

dr_h = 0.0127; % diameter of the tube (m)
Lr = 0.3048; % length of the tube (m)
rho_bulk = 800; % bulk density of the catalyst (kg/m3)

% this function calculates the relative sensitivity of the output variables
% to parameters at given condition. The output returns a raw relative
% sensitivity matrix and a percent relative sensitivity matrix. 

% rows=parameters, columns=species
[rel_sens_matrix, percent_rel_sens] = calculateRelativeSensitivity(T, ...
    ghsv, ratio, lsq_parameters, delta_h, mw, Ts, rhos, dr_h, Lr, rho_bulk);