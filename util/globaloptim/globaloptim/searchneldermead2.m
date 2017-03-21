function [successSearch,xBest,fBest,funccount] = searchneldermead2(FUN,X,Aineq,bineq, ...
                Aeq,beq,lb,ub,optimValues,options,iterLimit,optionsNM)
%SEARCHNELDERMEAD2 PATTERNSEARCH2 optional search2 step using FMINSEARCH.
%   [SUCCESSSEARCH,NEXTITERATE,OPTIMSTATE] = SEARCHNELDERMEAD2(FUN,X,A,b,Aeq, ...
%   beq,lb,ub,OPTIMVALUES,OPTIONS,ITERLIMIT,OPTIONSNM) searches for a lower
%   function value uses the Nelder-Mead method.
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
%   OPTIONSNM: Options structure for FMINSEARCH (Optional argument)
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
 
if nargin < 12 || isempty(optionsNM) % Short-circuit will prevent 2nd test if 1st failed
   optionsNM = optimset('Display','off');
end
if nargin < 11 || isempty(iterLimit) % Short-circuit will prevent 2nd test if 1st failed
    iterLimit = 1;
end
% Initialize output
successSearch = 0;
xBest = X;
fBest = optimValues.fval;
funccount = optimValues.funccount;

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
% Call FMINSEARCH 
optionsNM = optimset(optionsNM,'MaxFunEvals',options.MaxFunEvals - optimValues.funccount);
[x,fval,unused1,output] = fminsearch(FUN,X,optionsNM,objFcnArg{:});
funccount = funccount + output.funcCount; % optim functions return funcCount

if fval < optimValues.fval
    successSearch = true;
    xBest = x(:);
    fBest = fval;
end
