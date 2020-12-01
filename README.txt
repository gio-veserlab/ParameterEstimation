This code bundle performs:

1. Kinetic parameter estimation of a simple model reaction (A+B<->C+D) from 4 different regression / optimization algorithms (least squares, cross-val least sq, Markov chain Monte Carlo, genetic algorithm)

2. Model evaluation via local sensitivity and identifiability (profile likelihood) analysis

Experimental concentration data (please format based on current data.xls), molecular weights and density values of each species in the reaction, reactor geometry info are required to run the scripts.

Brief descriptions of each script is provided below:

mainParameterEstimation: Main script that calls and runs all parameter estimation methods and model evaluation functions

fitXX: fitting routines of each algorithm

calculateProfileLikelihood: Generates likelihood profiles of each kinetic parameter. Must be run after parameters are estimated. 

calculateRelativeSensitivity: Calculates relative sensitivities of each species to kinetic parameters. Must be run after parameters are estimated. 

Support folder:

steadyReactor: Defines the material balance equations

runModel: Runs the steady-state reactor and calculates exit concentrations

modelError: Defines error functions to be minimized during kinetic fittings. Used by all methods.

fitProfileLikelihood: Parameter re-fitting script based on nonlinear least squares. Used only for identifiability analysis. Difference from regular least squares is that one of the kinetic parameters is always fixed during refit.

Sample outputs are provided in the demo.