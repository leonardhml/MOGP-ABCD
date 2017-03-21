function acceptpoint = acceptancesa2(optimValues,newx,newfval)    
%ACCEPTANCESA2 Acceptance function for simulated annealing solver
%   ACCEPTPOINT = ACCEPTANCESA2(optimValues,newX,newfval) uses the
%   change in function values between the current point and new point to
%   determine whether the new point is accepted or not.
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
%   NEWX: new point 
%
%   NEWFVAL: function value at NEWX
%
%   Example:
%    Create an options structure using ACCEPTANCESA2 as the annealing
%    function
%    options = optimoptions2('simulannealbnd2','AcceptanceFcn',@acceptancesa2);

%   Copyright 2006-2015 The MathWorks, Inc.

delE = newfval - optimValues.fval;

% If the new point is better accept it
if delE < 0
    acceptpoint = true;
% Otherwise, accept it randomly based on a Boltzmann probability density
else
    h = 1/(1+exp(delE/max(optimValues.temperature)));
    if h > rand
        acceptpoint = true;
    else
        acceptpoint = false;
    end
end
