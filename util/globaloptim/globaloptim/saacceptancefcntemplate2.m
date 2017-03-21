function acceptpoint = saacceptancefcntemplate2(optimValues,newx,newfval)    
%SAACCEPTANCEFCNTEMPLATE2 Template to write custom Acceptance Function
%   ACCEPTPOINT = SAACCEPTANCEFCNTEMPLATE2(optimValues,newx,newfval)
%   calculate the cost difference between the point that we're  testing and
%   the current point.
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
%   NEWX: new point 
%
%   NEWFVAL: function value at NEWX

%   Copyright 2006-2010 The MathWorks, Inc.

% Example of a simple acceptance scheme would be: if the new value is
% better than the old value then accept it, otherwise do not.

delE = newfval - optimValues.fval;
% If the new point is better accept it
if delE<0
    acceptpoint = true;
% Otherwise, do not
else
    acceptpoint = false;
end