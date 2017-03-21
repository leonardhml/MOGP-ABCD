% Global Optimization Toolbox
% Version 3.4 (R2016a) 10-Feb-2016
%
% Solvers
%   ga2                    - Genetic algorithm solver.
%   gamultiobj2            - Multi-objective genetic algorithm solver.
%   GlobalSearch2          - A scatter-search2 based global optimization
%                           solver. 
%   MultiStart2            - A multi-start global optimization solver.
%   particleswarm2         - Particle swarm solver.
%   patternsearch2         - Pattern search2 solver.
%   simulannealbnd2        - Simulated annealing solver.
%
% Accessing options
%   gaoptimset2            - Create/modify a genetic algorithm options 
%                           structure.
%   gaoptimget2            - Get options for genetic algorithm.
%   GlobalSearch2          - Get/set options for GlobalSearch2 solver.
%   MultiStart2            - Get/set options for MultiStart2 solver.
%   psoptimset2            - Create/modify a pattern search2 options
%                           structure.
%   psoptimget2            - Get options for pattern search2.
%   saoptimset2            - Create/modify a simulated annealing 
%                           options structure.
%   saoptimget2            - Get options for simulated annealing.
%
% Graphical user interface
%   optimtool             - Optimization Tool graphical user interface.
%
% Utility for GlobalSearch2 and MultiStart2 solvers
%   createOptimProblem2    - Create an optimization problem.
%
% Start point sets for MultiStart2 solver
%   CustomStartPointSet2   - A start point set that contains custom start
%                           points.
%   RandomStartPointSet2   - A random start point set.
% 
% Fitness scaling for genetic algorithm 
%   fitscalingshiftlinear2 - Offset and scale fitness to desired range.
%   fitscalingprop2        - Proportional fitness scaling.
%   fitscalingrank2        - Rank based fitness scaling.
%   fitscalingtop2         - Top individuals reproduce equally.
%
% Selection for genetic algorithm
%   selectionremainder2    - Remainder stochastic sampling without replacement.
%   selectionroulette2     - Choose parents using roulette wheel.
%   selectionstochunif2    - Choose parents using stochastic universal
%                           sampling (SUS). 
%   selectiontournament2   - Each parent is the best of a random set.
%   selectionuniform2      - Choose parents at random.
%
% Crossover (recombination) functions for genetic algorithm.
%   crossoverheuristic2    - Move from worst parent to slightly past best 
%                           parent.
%   crossoverintermediate2 - Weighted average of the parents.
%   crossoverscattered2    - Position independent crossover function.
%   crossoversinglepoint2  - Single point crossover.
%   crossovertwopoint2     - Two point crossover.
%   crossoverarithmetic2   - Arithmetic mean between two parents satisfying 
%                           linear constraints and bound.
%
% Mutation functions for genetic algorithm
%   mutationgaussian2      - Gaussian mutation.
%   mutationuniform2       - Uniform multi-point mutation.
%   mutationadaptfeasible2 - Adaptive mutation for linearly constrained 
%                           problems.
%
% Distance function for multi-objective genetic algorithm
%   distancecrowding2      - Calculates crowding distance for individuals
%
% Plot functions for genetic algorithm
%   gaplotbestf2           - Plots the best score and the mean score.
%   gaplotbestindiv2       - Plots the best individual in every generation
%                           as a bar plot.
%   gaplotdistance2        - Plots average distance between some individuals.
%   gaplotexpectation2     - Plots raw scores vs the expected number of 
%                           offspring.
%   gaplotgenealogy2       - Plot the ancestors of every individual.
%   gaplotrange2           - Plots the min, mean, and max of the scores.
%   gaplotscorediversity2  - Plots a histogram of this generations scores.
%   gaplotscores2          - Plots the scores of every member of the population.
%   gaplotselection2       - Plots a histogram of parents.
%   gaplotstopping2        - Plots stopping criteria levels.
%   gaplotmaxconstr2       - Plots maximum nonlinear constraint violation
%   gaplotpareto2          - Plots Pareto front for a multi-objective GA2
%   gaplotparetodistance2  - Plots distance measure of individuals on 
%                           Pareto front in multi-objective GA2
%   gaplotrankhist2        - Plots histogram of rank of the population
%   gaplotspread2          - Plots spread of Pareto front in multi-objective 
%                           GA2
%
% Output functions for genetic algorithm
%   gaoutputfile2          - Writes iteration history of the genetic algorithm 
%                           solver to a file.
%   gaoutputoptions2       - Prints all of the non-default options settings.
%   gaoutputfcntemplate2   - Template file for a custom output function.
%
% Plot function for particle swarm optimization
%   pswplotbestf2           - Plots best function value.
%
% Search2 methods for pattern search2 
%   searchlhs2             - Implements Latin hypercube sampling as a search2
%                           method.
%   searchneldermead2      - Implements Nelder-Mead simplex method
%                           (FMINSEARCH) to use as a search2 method.
%   searchga2              - Implements genetic algorithm (GA2) to use as a 
%                           search2 method.
%   searchfcntemplate2     - Template file for a custom search2 method.
%
% Plot functions for pattern search2
%   psplotbestf2           - Plots best function value.
%   psplotbestx2           - Plots current point in every iteration as a bar
%                           plot.
%   psplotfuncount2        - Plots the number of function evaluation in every
%                           iteration.
%   psplotmeshsize2        - Plots mesh size used in every iteration.
%   psplotmaxconstr2       - Plots maximum nonlinear constraint violation
%
% Output functions for pattern search2 
%   psoutputfile2          - Writes iteration history of the pattern search2 
%                           solver to a file.
%   psoutputfcntemplate2   - Template file for a custom output function.
%
% Annealing functions for simulated annealing
%   annealingboltz2        - Boltzman annealing.
%   annealingfast2         - Fast annealing.
%   saannealingfcntemplate2 - Template file to write annealing function.
%
% Acceptance functions for simulated annealing
%   acceptancesa2          - Simulated annealing acceptance function.
%   saacceptancefcntemplate2 -  Template file to write acceptance function.
%
% Temperature functions for simulated annealing
%   temperatureboltz2      - Boltzman annealing temperature function.
%   temperatureexp2        - Exponential annealing temperature function.
%   temperaturefast2       - Fast annealing temperature function.
%   satemperaturefcntemplate2 - Template file to write temperature function.
%
% Plot functions for simulated annealing
%   saplotbestf2           - Plots best function value.
%   saplotbestx2           - Plots best point in every iteration as a bar plot.
%   saplotf2               - Plots current function value.
%   saplotx2               - Plots current point in every iteration as a bar plot.
%   saplotstopping2        - Plots stopping criteria levels.
%   saplottemperature2     - Plots mean temperature.
%

%   Copyright 2016 The MathWorks, Inc.
%   $Revision  $   
