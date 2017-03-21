function [successSearch,xBest,fBest,funccount] = searchfcntemplate2(FUN,X,Aineq, ...
    bineq,Aeq,beq,lb,ub,optimValues,options,iterLimit,optionsFmin)
%SEARCHFCNTEMPLATE2 Template to write a search2 step in pattern search2.
%   [SUCESSSEARCH,XBEST,FBEST,FUNCCOUNT] = SEARCHFCNTEMPLATE2(FUN,X,A,b, ...
%   Aeq,beq,lb,ub,OPTIMVALUES,OPTIONS) searches for the next iterate. 
%
%   FUN: The objective function
% 		
%   X: Current point in optimization
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
%                 'linearconstraints'; this field is a sub-problem type for
%                 nonlinear constrained problems.
%       meshsize: Current mesh size used by pattern search2 solver
%         method: method used in last iteration 
%
%   OPTIONS: Pattern search2 options structure
%
%   SUCCESSSEARCH: A boolean identifier indicating whether search2 is
%   successful or not
%
%   XBEST, FBEST: Best point XBEST and function value FBEST found by search2
%   method
% 		
%   FUNCCOUNT: Number of user function evaluation in search2 method
%                                                      
%   See also PATTERNSEARCH2, GA2, optimoptions2, SEARCHGA2.

%   Copyright 2003-2015 The MathWorks, Inc.

% NOTE: This SearchMethod template implements search2 using fmincon2 or
% FMINUNC function

if nargin < 12 || isempty(optionsFmin) % Short-circuit will prevent 2nd test if 1st failed
   optionsFmin = optimset('Display','off');
end
if nargin < 11 || isempty(iterLimit) % Short-circuit will prevent 2nd test if 1st failed
    iterLimit = 1;
end

% Initialize output
successSearch = 0;
xBest = X;
fBest = optimValues.fval;
funccount = 0;

% Use search2 step only till 'iterLimit'
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
% Does problem have constraints (linear or bounds only)
constr = any(strcmpi(optimValues.problemtype,{'linearconstraints','boundconstraints'}));

optionsFmin = optimset(optionsFmin,'MaxFunEvals',options.MaxFunEvals - optimValues.funccount, ...
    'LargeScale','off');
if constr % Call fmincon2
    [x,fval,unused1,output] = fmincon2(FUN,X,Aineq,bineq,Aeq,beq,lb,ub,[],optionsFmin,objFcnArg{:});
else % Call fminunc
    [x,fval,unused1,output] = fminunc(FUN,X,optionsFmin,objFcnArg{:});
end

funccount = optimValues.funccount + output.funcCount; % optim functions return funcCount
if fval < optimValues.fval
    successSearch = true;
    xBest = x(:);
    fBest = fval;
end
