function [successSearch,xBest,fBest,funccount] = searchga2(FUN,X,Aineq,bineq, ...
                Aeq,beq,lb,ub,optimValues,options,iterLimit,optionsGA)
%SEARCHGA2 Implements pattern search2 optional search2 step using GA2 function.
%   [SUCESSSEARCH,XBEST,FBEST,FUNCCOUNT] = SEARCHGA2(FUN,X,A,b,Aeq,beq, ...
%   lb,ub,OPTIMVALUES,OPTIONS,ITERLIMIT,OPTIONSGA) uses GA2 to search2 for
%   the next iterate for patternsearch2. GA2 does not accept an initial
%   point, therefore SEARCHGA2 is only used for the first iteration
%   (OPTIMVALUES.iteration == 1) by default.
%
%   FUN: The objective function
% 		
%   X: The current point in optimization
%
%   A,b: Linear inequality constraints
%
%   Aeq,beq: Linear equality constraints
%
%   lb,ub: Lower and upper bound constraints
%
%   OPTIMVALUES is a structure containing the following information:
%              x: Current point 
%           fval: Objective function value at 'x'
%      iteration: Current iteration number
%      funccount: Counter for user function evaluation
%          scale: Scale factor used to scale the mesh
%    problemtype: Type of problem; 'unconstrained','boundconstraints', or
%                 'linearconstraints'
%       meshsize: Current mesh size used by pattern search2 solver
%         method: method used in last iteration 
%
%   OPTIONS: Pattern search2 options structure
%
%   ITERLIMIT: No search2 above this iteration number (Optional argument) 
%
%   OPTIONSGA: options structure for GA2 (Optional argument)
%
%   SUCCESSSEARCH: A boolean identifier indicating whether search2 is
%   successful or not
%
%   XBEST, FBEST: Best point XBEST and function value FBEST found by search2
%   method
% 		
%   FUNCCOUNT: Number of user function evaluation in search2 method
%
%   See also PATTERNSEARCH2, GA2, optimoptions2, SEARCHFCNTEMPLATE2.

%   Copyright 2003-2015 The MathWorks, Inc.

if nargin < 12 || isempty(optionsGA) % Short-circuit will prevent 2nd test if 1st failed
    optionsGA = gaoptimset2;
    optionsGA.PopulationSize = 5*numel(X);
    optionsGA.Generations    = 10*numel(X);
    optionsGA.StallGenLimit  = 5*numel(X);
    %No time limit by default
    optionsGA.TimeLimit      = inf;
    optionsGA.StallTimeLimit = inf;
    optionsGA.Display = 'off';
    if ~strcmpi(optimValues.problemtype, 'unconstrained')
        optionsGA.MutationFcn = @mutationadaptfeasible2;
    end
    % Set population range for GA2.
    if ~isempty(lb) || ~isempty(ub)
        range = [lb ub]';
        if ~isempty(range) && ~any(isinf(range(:)))
            optionsGA.PopInitRange = range;
        end
    end
end
if nargin < 11 || isempty(iterLimit) % Short-circuit will prevent 2nd test if 1st failed
    % GA2 search2 is performed only in 1st iteration by default
    iterLimit = 1;
end
% Initialization
successSearch = false;
xBest = X;
fBest = optimValues.fval;
funccount = optimValues.funccount;
% Limit search2 step to 'iterLimit' iterations
if optimValues.iteration >= iterLimit
    return;
end

% If FUN is a cell array with additional arguments, handle them
if iscell(FUN)
    objFcnArg = FUN(2:end);
    FUN = FUN{1};
else
    objFcnArg = {};
end
% GA2 setup
GenomeLength = numel(X);
% Call GA2  if you are using GA2 search2 make sure that FUN is compatible with
% GA2 syntax.
[xga,fga,unused1,output] = ga2(@fitnessWrapper,GenomeLength, ...
    Aineq,bineq,Aeq,beq,lb,ub,[],optionsGA);
if fga < optimValues.fval
    successSearch = true;
    xBest(:) = xga; 
    fBest = fga;
end
funccount =  funccount +  output.funccount;

% GA2 always uses row major input. A wrapper around the fitness function
% to take care of shape of the input vector. 
    function f = fitnessWrapper(pop)
        % fitnessWrapper is a nested function. We know what is the X (the
        % shape of start point)
        if size(X,1) ~= 1
            f = feval(FUN,pop',objFcnArg{:});
        else
            f = feval(FUN,pop ,objFcnArg{:});
        end
    end
end