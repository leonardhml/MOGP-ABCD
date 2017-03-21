function [successSearch,xBest,fBest,funccount] = searchlhs2(FUN,X,Aineq,bineq, ...
                Aeq,beq,lb,ub,optimValues,options,iterLimit,factors)
%SEARCHLHS2 PATTERNSEARCH2 optional search2 step using LHS.
%   [SUCESSSEARCH,XBEST,FBEST,FUNCCOUNT] = SEARCHLHS2(FUN,X,A,b,Aeq,beq, ...
%   lb,ub,OPTIMVALUES,OPTIONS,ITERLIMIT,FACTORS) searches for the next
%   iterate in a Latin hypercube design space. SEARCHLHS2 is only used for
%   the first iteration (OPTIMVALUES.iteration == 1) by default.
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
%   FACTORS: The design level for LHS (Optional argument)
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

successSearch = 0;
xBest = X;
fBest = optimValues.fval;
funccount = optimValues.funccount;
numberofVariables  = numel(X);
range = [];

if nargin < 12 || isempty(factors) % Short-circuit will prevent 2nd test if 1st failed
    factors = 15*numberofVariables; 
end
if nargin < 11 || isempty(iterLimit) % Short-circuit will prevent 2nd test if 1st failed
    iterLimit = 1;  
end
   

% Use search2 step only till iterLimit.
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
constr = any(strcmpi(optimValues.problemtype, ...
    {'linearconstraints','nonlinearconstr','boundconstraints'}));

% LHS design calculations for 'maximin' criterion
Points = lhspoint2(factors,numberofVariables);
if ~isempty(lb) || isempty(ub)
    range = [lb ub]';
end

if ~isempty(range) && ~any(isinf(range(:)))
    limit = range(2,:) - range(1,:);
    Points = repmat(range(1,:)',1,size(Points,2)) +  repmat(limit',1,size(Points,2)) .* Points;
    span = size(Points,2);
    sites = struct('x',cell(span,1,1),'f',cell(span,1,1));
    for k = 1:span
        sites(k).x = Points(:,k);
    end
else 
    Points = -1 + Points*2;
    span = size(Points,2);
    sites = struct('x',cell(span,1,1),'f',cell(span,1,1));
    for k = 1:span
        sites(k).x = X(:) + optimValues.meshsize*optimValues.scale.*Points(:,k); 
    end
end

% Create structure 'Iterate' to be used by local function 'psnextfeasible2'
Iterate.x = X(:); Iterate.f = optimValues.fval;

% Find an iterate with lower objective using psnextfeasible2 utility. This
% utility works for Poll2 options so, override the 'CompletePoll' option by
% 'CompleteSearch' option. Also, override 'NotVectorizedPoll' option by
% 'NotVectorizedSearch'
options.CompletePoll = options.CompleteSearch;
options.NotVectorizedPoll = options.NotVectorizedSearch;
[successSearch,Iterate,~,funccount] = psnextfeasible2(FUN,X,sites, ...
    Iterate,Aineq,bineq,Aeq,beq,lb,ub,options.TolBind,constr,objFcnArg,optimValues.funccount,options);

if successSearch
    xBest(:) = Iterate.x;
    fBest = Iterate.f;
end
