function [X,FVAL,EXITFLAG,OUTPUT] = pfminbnd2(FUN,objFcnArg,initialX,numberOfVariables,Iterate,lb,ub,options,defaultopt,OUTPUT)
%PFMINBND2 Finds minimum of a function with bound constraints.
%   PFMINBND2 solves problems of the form:
%        min F(X)  subject to: LB <= X <= UB  (Box constraints)
%         X
%
%   Private to PATTERNSEARCH2

%   Copyright 2003-2015 The MathWorks, Inc.

% Get some initial values
[optimState,nextIterate,MeshSize,EXITFLAG,run] = ...
    getinitial2(Iterate,numberOfVariables,0,lb,ub,options);

X = initialX;
X(:) = Iterate.x;
FVAL = Iterate.f;
% Determine who is the caller
callStack = dbstack;
[~,caller] = fileparts(callStack(2).file);

% Call output and plot functions
if options.OutputTrue || options.PlotTrue
    % Set state for plot and output functions (only pfmincon2 will have
    % 'interrupt' state)
    if ~strcmp(caller,'pfmincon2')
        currentState = 'init';
    else
        currentState = 'interrupt';
    end
    callOutputPlotFunctions(currentState)
end

% Setup display header
if  options.Verbosity > 1
    fprintf('\n\nIter     f-count          f(x)      MeshSize     Method\n');
end
% Set state for plot and output functions (only pfmincon2 will have
% 'interrupt' state)
if ~strcmp(caller,'pfmincon2')
    currentState = 'iter';
else
    currentState = 'interrupt';
end

while run
    % Check for convergence
    [X,EXITFLAG,FVAL,msg,run] = isconverged2(optimState,options,MeshSize, ...
        nextIterate,X,EXITFLAG,run);
    if ~run
        continue;
    end
    % SEARCH2.
    [successSearch,nextIterate,optimState] = search2(FUN,X,Iterate,MeshSize,[],[], ...
        [],[],eye(numberOfVariables),lb,ub,OUTPUT.problemtype,objFcnArg,optimState,options);
    % POLL2
    if ~successSearch  % Unsuccessful search2
        [successPoll,nextIterate,optimState] = poll2(FUN,X,Iterate,MeshSize,[],[], ...
            [],[],eye(numberOfVariables),lb,ub,OUTPUT.problemtype,objFcnArg,optimState,options);
    else
        successPoll =0;
    end

    % Scale the variables (if needed)
    if any(strcmpi(options.ScaleMesh,{'dynamic','on'}))
        meanX = mean([Iterate.x],2);
        optimState.scale = logscale2(lb,ub,meanX);
    end

    % Update
    [MeshSize,Iterate,X,optimState] = updateparam2(successPoll,successSearch, ...
        MeshSize,nextIterate,Iterate,X,optimState,options);
    
    % Call output and plot functions
    if options.OutputTrue || options.PlotTrue
        callOutputPlotFunctions(currentState)
    end
end
% Call output and plot functions
if options.OutputTrue || options.PlotTrue
    % Set state for plot and output functions (only pfmincon2 will have
    % 'interrupt' state)
    if ~strcmp(caller,'pfmincon2')
        currentState = 'done';
    else
        currentState = 'interrupt';
    end
    callOutputPlotFunctions(currentState)
end

% Update values of OUTPUT structure
OUTPUT.pollmethod = options.PollMethod; % This might change via output function
OUTPUT.searchmethod = options.SearchMethod; % This might change via output function
OUTPUT.iterations = optimState.Iter;
OUTPUT.funccount = optimState.FunEval;
OUTPUT.meshsize = MeshSize;
OUTPUT.maxconstraint = max([X(:) - ub(:); lb(:) - X(:); 0]);
OUTPUT.message = msg;

%-----------------------------------------------------------------
% Nested function to call output/plot functions
    function callOutputPlotFunctions(state)
        optimvalues.x = X;
        optimvalues.iteration = optimState.Iter;
        optimvalues.fval = Iterate.f;
        optimvalues.problemtype = OUTPUT.problemtype;
        optimvalues.meshsize = MeshSize;
        optimvalues.funccount = optimState.FunEval;
        optimvalues.method = optimState.how;
        optimvalues.TolFun = optimState.deltaF;
        optimvalues.TolX = optimState.deltaX;
        solverName = 'Pattern Search2';
        switch state
            case {'init', 'iter'}
                if options.PlotTrue
                    optimState.stopPlot = gadsplot2(options,optimvalues,state,solverName);
                end
                if options.OutputTrue
                    [optimState.stopOutput,options] = psoutput2(...
                        options.OutputFcns, options.OutputFcnsArg, ...
                        optimvalues, options.OutputPlotFcnOptions, state, ...
                        options, defaultopt, numberOfVariables);
                end
            case 'interrupt'
                if options.PlotTrue
                    optimState.stopPlot = gadsplot2(options,optimvalues,state,solverName);
                end
                if options.OutputTrue
                    optimState.stopOutput = psoutput2(options.OutputFcns,options.OutputFcnsArg, ...
                        optimvalues,options.OutputPlotFcnOptions,state);
                end
            case 'done'
                if options.PlotTrue
                    gadsplot2(options,optimvalues,state,solverName);
                end
                if options.OutputTrue
                    psoutput2(options.OutputFcns,options.OutputFcnsArg,optimvalues,options.OutputPlotFcnOptions,state);
                end
        end
    end % End of callOutputPlotFunctions
%------------------------------------------------------------------
end  % End of pfminbnd2

