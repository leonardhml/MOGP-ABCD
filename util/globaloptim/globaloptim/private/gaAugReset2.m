function  innerOptions = gaAugReset2(innerPopulation,options,state,gaAugParam)
%gaAugReset2 Create options structure for sub-problem in ALGA.
%   [exitFlag,innerState,options,step] = gaAugReset2(Iterate, ...
%     state,options,step,currentTolFun)

%   Copyright 2005-2015 The MathWorks, Inc.

% Initialize options structure used by sub-problem solver
innerOptions = options;
% a variety of data used in various places
innerOptions.TimeLimit = options.TimeLimit - toc(state.StartTime);
innerOptions.FitnessLimit = -Inf;
innerOptions.TolFun = gaAugParam.currentTolFun;  % Use the current value of tolerance
% Other settings
% Modify display if a sub-problem is solved
innerOptions.Verbosity = 0;
innerOptions.Display = 'off';
innerOptions.HybridFcn = [];

% InitialPopulation can be set for all generations > 1
if ~isempty(innerPopulation)
    innerOptions.InitialPopulation = innerPopulation;
end
% Users will not be able to provide initial scores because ALGA will 
% optimize a sub-problem formulation not fitness function. 
innerOptions.InitialScores = [];

