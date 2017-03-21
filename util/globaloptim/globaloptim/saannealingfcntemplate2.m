function newx = saannealingfcntemplate2(optimValues,problem)
%SAANNEALINGFCNTEMPLATE2 Template to write annealing function
%   NEWX = SAANNEALINGFCNTEMPLATE2(optimValues,problem) generate a point
%   based on the current point and the current temperature
%
%   OPTIMVALUES is a structure containing the following information:
%              x: current point 
%           fval: function value at x
%          bestx: best point found so far
%       bestfval: function value at bestx
%    temperature: current temperature
%      iteration: current iteration 
%             t0: start time
%              k: annealing parameter 'k'
%
%   PROBLEM is a structure containing the following information:
%      objective: function handle to the objective function
%             x0: the start point
%           nvar: number of decision variables
%             lb: lower bound on decision variables
%             ub: upper bound on decision variables
%

%   Copyright 2006-2010 The MathWorks, Inc.

% A simple annealing function would be generate new point using a uniform
% distribution

currentx = optimValues.x;
nvar = numel(currentx);
newx = currentx;
newx(:) = currentx(:) + sqrt(optimValues.temperature).*rand(nvar,1);
% New point must be in bounds
newx = sahonorbounds2(newx,optimValues,problem);