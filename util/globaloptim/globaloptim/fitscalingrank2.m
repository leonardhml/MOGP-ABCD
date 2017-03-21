function expectation = fitscalingrank2(scores,nParents)
% FITSCALINGRANK2 Rank based fitness scaling (single objective only).
%   EXPECTATION = FITSCALINGRANK2(SCORES,NPARENTS) calculates the
%   EXPECTATION using the SCORES and number of parents NPARENTS.
%   This relationship can be linear or nonlinear.
%
%   Example:
%   Create an options structure using FITSCALINGRANK2 as 
%   the fitness scaling function
%     options = optimoptions2('ga2','FitnessScalingFcn',@fitscalingrank2);

%   Copyright 2003-2015 The MathWorks, Inc.

scores = scores(:);
[~,i] = sort(scores);

expectation = zeros(size(scores));
expectation(i) = 1 ./ ((1:length(scores))  .^ 0.5);

expectation = nParents * expectation ./ sum(expectation);
