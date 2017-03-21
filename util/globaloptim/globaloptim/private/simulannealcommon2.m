function [x,fval,exitflag,output,solverData,problem,options] = ...
    simulannealcommon2(FUN,x0,Aineq,bineq,Aeq,beq,lb,ub,options,defaultopt)
%simulannealcommon2 is called by simulanneal2 to perform common
%initialization tasks

% This function is private to simulanneal2.

%   Copyright 2006-2015 The MathWorks, Inc.

% Only function_handle are allowed
if isempty(FUN) ||  ~isa(FUN,'function_handle')
    error(message('globaloptim:simulannealcommon2:needFunctionHandle'));
end

% Use default options if empty
if ~isempty(options) && ~isa(options,'struct')
    error(message('globaloptim:simulannealcommon2:fifthInputNotStruct'));
elseif isempty(options)
    options = defaultopt;
end
% All inputs should be double
try
    dataType = superiorfloat(x0,lb,ub);
    if ~isequal('double', dataType)
        error(message('globaloptim:simulannealcommon2:dataType'))
    end
catch
    error(message('globaloptim:simulannealcommon2:dataType'))
end

% Keep a copy of the user options structure
user_options = options;

% Get all default options and merge with non-default options
options = saoptimset2(defaultopt,options);

if ~isempty(x0)
    x = x0;
    initialX = x0(:);
    numberOfVariables = length(initialX);
else
    error(message('globaloptim:simulannealcommon2:initialPoint'));
end

% Remember the random number states used
output.iterations = 0;
output.funccount   = 0;
output.message   = '';
dflt = RandStream.getGlobalStream;
output.rngstate = struct('state',{dflt.State}, 'type',{dflt.Type});

% Initialize data structures for problem and also for solver state
problem = struct('objective', FUN, 'x0', x0, 'nvar', numberOfVariables);
solverData  = struct('t0',cputime,'currentx',initialX);

problem.bounded = any(isfinite(ub)) || any(isfinite(lb));

% Call savalidate2 to ensure that all the fields of options are valid
options = savalidate2(options,problem);

% Make sure that the bound-constraints are valid
[lb, ub, msg, exitflag] = checkbound2(lb, ub, numberOfVariables);
if exitflag < 0
    fval = [];
    x(:) = solverData.currentx;
    output.message = msg;
    if options.Verbosity > 0
        fprintf('%s\n',msg)
    end
    return;
end
problem.lb = lb; problem.ub = ub;
% Make sure the problem is not unbounded
if problem.bounded
    if strcmpi(options.DataType,'custom')
        problem.bounded = false;
        if options.Verbosity > 0
            warning(message('globaloptim:simulannealcommon2:boundsWithCustom'));
        end
    end
end
% Determine problemtype
if problem.bounded
    output.problemtype = 'boundconstraints';
else
    output.problemtype = 'unconstrained';
end

% For calling an OutputFcn or PlotFcn, make an options object (to be
% updated later) that can be passed in
if ~isempty(options.OutputFcns) || ~isempty(options.PlotFcns)
    options.OutputPlotFcnOptions = copyForOutputAndPlotFcn(...
        optim.options.SimulannealbndOptions2, options);
else
    options.OutputPlotFcnOptions = [];
end

% Aineq, bineq, Aeq, beq are not part of problem because linear
% constraints are not implemented yet
% Find initial feasible point
[solverData.currentx,~,~,~,~,~,~,msg,exitflag] = ...
    preProcessLinearConstr2(solverData.currentx,Aineq,bineq,Aeq,beq,lb,ub,0, ...
                    numberOfVariables,output.problemtype,options.Verbosity);
if exitflag < 0
    fval =    [];
    x(:) =  solverData.currentx;
    output.message = msg;
    if options.Verbosity > 0
        fprintf('%s\n',msg)
    end
    return
end
% Print diagnostic information if asked
if options.Verbosity > 2
    sadiagnose2(user_options,problem);
end
% Additional fields for solver parameters need initialization
solverData = samakedata2(solverData,problem,options);
% Initialize fval
fval = [];

