function expectation = fitscalingprop2(scores,nParents)
%FITSCALINGPROP2 Proportional fitness scaling (single objective only).
%   EXPECTATION = FITSCALINGPROP2(SCORES,NPARENTS) calculates the EXPECTATION
%   using the SCORES and number of parents NPARENTS. The least sophisticated 
%   fitness scaling strategy, expectation is simply proportional to the raw 
%   fitness score. This strategy has weaknesses when raw scores are not in a
%   "good range".
%
%   Example:
%   Create an options structure using FITSCALINGPROP2 as the fitness scaling
%   function
%     options = optimoptions2('ga2','FitnessScalingFcn',@fitscalingprop2);
%

%   Copyright 2003-2015 The MathWorks, Inc.

scores = scores(:);
% Rotate scores around their mean because we are minimizing so lower scores
% must yield a higher expectation.
scores = 2 * mean(scores) - scores;

% Negative expectations are simply not allowed so here we make sure that
% doesn't happen by sliding the whole vector up until everything is
% non-negative.
m = min(scores);
if(m < 0)
    scores = scores - m;
end

% Normalize: expectation should sum to one.
expectation = nParents * scores ./ sum(scores);
