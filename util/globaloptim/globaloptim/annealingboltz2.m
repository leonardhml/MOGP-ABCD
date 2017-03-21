function newx = annealingboltz2(optimValues,problem)
%ANNEALINGBOLTZ2 Generates a point using multivariate normal distribution.
%   NEWX = ANNEALINGBOLTZ2(optimValues,problem) generates a point based
%   on the current point and the current temperature using multivariate
%   normal distribution.
%
%   OPTIMVALUES is a structure containing the following information:
%              x: current point 
%           fval: function value at x
%          bestx: best point found so far
%       bestfval: function value at bestx
%    temperature: current temperature
%      iteration: current iteration 
%             t0: start time
%              k: annealing parameter
%
%   PROBLEM is a structure containing the following information:
%      objective: function handle to the objective function
%             x0: the start point
%           nvar: number of decision variables
%             lb: lower bound on decision variables
%             ub: upper bound on decision variables
%
%   Example:
%    Create an options structure using ANNEALINGBOLTZ2 as the annealing
%    function
%    options = optimoptions2('simulannealbnd2','AnnealingFcn',@annealingboltz2);

%   Copyright 2006-2015 The MathWorks, Inc.

currentx = optimValues.x;
nvar = numel(currentx);
newx = currentx;
y = randn(nvar,1);
y = y./norm(y);
newx(:) = currentx(:) + sqrt(optimValues.temperature).*y;

newx = sahonorbounds2(newx,optimValues,problem);
