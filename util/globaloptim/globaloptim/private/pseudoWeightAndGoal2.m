function [pweight, pgoal] = pseudoWeightAndGoal2(fval,scores)
%pseudoWeightAndGoal2 Calculate pseudo weight and goal.
%
%   This function is private to GAMULTIOBJ2

%   Copyright 2007 The MathWorks, Inc.

fmin = min(scores,[],1);
fmax = max(scores,[],1);
pgoal = fmin;

totalwt = sum((fmax - fval)./(1 + fmax - fmin));
pweight = (fmax - fval)./((1 + fmax - fmin)*totalwt);
