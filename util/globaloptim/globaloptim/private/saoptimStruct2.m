function optimvalues = saoptimStruct2(solverData,problem)
%saoptimStruct2 create a structure to be passed to user functions
%   OPTIMVALUES = saoptimStruct2(solverData,problem) creates a structure
%   OPTIMVALUES with the following fields:
%              x: current point 
%           fval: function value at x
%          bestx: best point found so far
%       bestfval: function value at bestx
%    temperature: current temperature
%      iteration: current iteration
%      funccount: number of function evaluations
%             t0: start time
%              k: annealing parameter
%
%   PROBLEM is a structure with the following fields:
%      objective: function handle to the objective function
%             x0: the start point
%           nvar: number of decision variables
%             lb: lower bound on decision variables
%             ub: upper bound on decision variables
%
%   solverData is a structure used by the solvers SIMULANNEALBND2. The
%   fields of this structure may change in future.  

%   This function is private to SIMULANNEAL2.

%   Copyright 2006-2010 The MathWorks, Inc.

% Prepare data to be sent over to temperature function
optimvalues.x = reshapeinput2(problem.x0,solverData.currentx);
optimvalues.fval = solverData.currentfval;
optimvalues.bestx = reshapeinput2(problem.x0,solverData.bestx);
optimvalues.bestfval = solverData.bestfval;
optimvalues.temperature = solverData.temp;
optimvalues.iteration = solverData.iteration;
optimvalues.funccount = solverData.funccount;
optimvalues.t0 = solverData.t0;
optimvalues.k = solverData.k;
