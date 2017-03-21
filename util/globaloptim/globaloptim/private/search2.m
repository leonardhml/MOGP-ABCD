function [successSearch,nextIterate,optimState] = search2(FUN,Xin,Iterate,MeshSize,Aineq,bineq, ...
    Aeq,beq,NullBasisAeq,lb,ub,problemtype,objFcnArg,optimState,options)
%SEARCH2 Implements a generic search2 step as described in GPS.
%        FUN: The objective function on which SEARCH2 step is implemented.
%
%        ITERATE: Incumbent point around which polling will be done. Iterate Stores
%        the current point 'x' and function value 'f' at this point.
%
%        MESHSIZE: Current mesh size used in SEARCH2 step.
%
%        SCALE: Scale factor used to scale the design points.
%
%        PROBLEMTYPE: This flag is passed to the SEARCH2 routines, indicating that the
%        problem is 'unconstrained', 'boundconstraints', 'linearconstraints'
%
%        NEXTITERATE: Successful iterate after polling is done. If POLL2 is NOT
%        successful, NEXTITERATE is same as ITERATE.
%
%        SUCCESSSEARCH: A boolean identifier indicating whether SEARCH2 is
%         successful or not.

%   Copyright 2003-2009 The MathWorks, Inc.


nextIterate=Iterate;
successSearch =false;
searchtype = options.SearchMethod;
constr = any(strcmpi(problemtype,{'linearconstraints','nonlinearconstr','boundconstraints'}));
feasible = true;

% If [], then do nothing
if isempty(searchtype)
    return;
end


% If not a function_handle
if ischar(searchtype)
    % Override 'PollMethod' option by 'SearchMethod' option
    options.PollMethod = searchtype;
    % Override the 'CompletePoll' option by 'CompleteSearch' option
    options.CompletePoll = options.CompleteSearch;
    % Override 'NotVectorizedPoll' option by 'NotVectorizedSearch'
    options.NotVectorizedPoll = options.NotVectorizedSearch;
    [successSearch,nextIterate,optimState] = poll2(FUN,Xin,Iterate,MeshSize,Aineq,bineq, ...
        Aeq,beq,NullBasisAeq,lb,ub,problemtype,objFcnArg,optimState,options);
    return;
elseif isa(searchtype,'function_handle') % A function handle
    prepareSearchFcnInputArg(); % Prepare search2 function input
    try
       [successSearch,nextIterate.x(:),nextIterate.f,optimState.FunEval] = feval(searchtype,FUN,Xin,Aineq,bineq, ...
            Aeq,beq,lb,ub,optimValues,options,options.SearchMethodArg{:});
     %  optimState.FunEval = optimState.FunEval + FunEval;
    catch userFcn_ME
        try
            % Try with old syntax
            [successSearch,nextIterate,optimState.FunEval] = callSearchFunction;
            % The old syntax worked but it will be obsolete in future version
            if optimState.Iter < 1
                warning(message('globaloptim:search2:oldSyntaxSearchMethod'));
            end
        catch unused_ME
            % Assume that the syntax is new and error happened otherwise
            gads_ME = MException('globaloptim:search2:InvalidSearchType', ... 
                'Failure in SearchMethod evaluation.');
            userFcn_ME = addCause(userFcn_ME,gads_ME);
            rethrow(userFcn_ME)
        end

    end
end

% This should not happen; make sure that when a custom
% search2 is used nextIterate is better than Iterate.
if constr
    feasible =   isTrialFeasible2(nextIterate.x,Aineq,bineq,Aeq,beq,lb,ub,options.TolBind);
end
if nextIterate.f > Iterate.f || ~feasible
    nextIterate = Iterate;
    successSearch =0;
end

%-------helper function to prepare input arguments to search2 function
    function prepareSearchFcnInputArg()
       if  ~isempty(objFcnArg) % Must be a cell array syntax
           FUN = {FUN,objFcnArg{:}};
       end
       optimValues.problemtype = problemtype;
       optimValues.x = Iterate.x;
       optimValues.fval = Iterate.f;
       optimValues.iteration = optimState.Iter;
       optimValues.scale = optimState.scale;
       optimValues.meshsize = MeshSize;
       optimValues.method = optimState.how;
       optimValues.funccount = optimState.FunEval;
       options.MaxFunEvals = max(0,options.MaxFunEvals - optimState.FunEval);
       Xin(:) = Iterate.x;

    end
%-------helper function to call search2 function with old syntax----------
    function [successSearch,nextIterate,FunEval] = callSearchFunction()
        % This function will call search2 functions using old syntax
        % Create A, L, U form for linear/bound constraints (old syntax)
    [unused1,A,L,U,unused2,unused3,unused4,IndIneqcstr,IndEqcstr] = ...
    aluform2(Iterate.x,Aineq,bineq,Aeq,beq,lb,ub,length(Iterate.x),problemtype,0,[],Xin,[]);

        searchOptions = struct('completesearch',options.CompleteSearch, ...
            'meshsize',MeshSize, ...
            'indineqcstr',IndIneqcstr, ...
            'indeqcstr',IndEqcstr, ...
            'problemtype',problemtype, ...
            'iteration',optimState.Iter, ...
            'scale',optimState.scale, ...
            'notvectorized',options.NotVectorizedSearch, ...
            'cache',options.Cache, ...
            'cachetol',options.CacheTol, ...
            'cachelimit',options.CacheSize);

        [successSearch,nextIterate,FunEval] = feval(searchtype,FUN,Xin,Iterate,options.TolBind, ...
            A,L,U,optimState.FunEval,options.MaxFunEvals,searchOptions,objFcnArg,options.SearchMethodArg{:});

    end %----------End of callSearchFunction-----------
end
