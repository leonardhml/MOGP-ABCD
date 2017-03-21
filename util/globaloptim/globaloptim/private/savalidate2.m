function options  = savalidate2(options,problem)
% SAVALIDATE2 is used to ensure that the options structure for a SIMULANNEAL2
%   problem is valid and each field in it has an acceptable type.
%
%   This function is private to SIMULANNEAL2.

%   Copyright 2006-2012 The MathWorks, Inc.

% Sanity check for the options structure
options = saoptimset2(options);

% Determine the verbosity
switch  options.Display
    case {'off','none'}
        options.Verbosity = 0;
    case 'final'
        options.Verbosity = 1;
    case 'iter'
        options.Verbosity = 2;
    case 'diagnose2'
        options.Verbosity = 3;
    otherwise
        options.Verbosity = 1;
end

% MaxFunEvals default value is string
if strcmpi(options.MaxFunEvals,'3000*numberofvariables')
    options.MaxFunEvals = 3000*problem.nvar;
end
% StallIterLimit default value is string
if strcmpi(options.StallIterLimit,'500*numberofvariables')
    options.StallIterLimit = 500*problem.nvar;
end

if strcmpi(options.HybridInterval, 'end')
    options.HybridInterval = [];
elseif strcmpi(options.HybridInterval, 'never')
    options.HybridFcn = [];
elseif ~isempty(options.HybridInterval)
    nonNegInteger2('HybridInterval', options.HybridInterval)
end

positiveInteger2('MaxIter', options.MaxIter)
positiveInteger2('MaxFunEvals', options.MaxFunEvals);
positiveInteger2('StallIterLimit', options.StallIterLimit)
positiveScalarArray2('InitialTemperature',options.InitialTemperature);

positiveScalar2('TimeLimit', options.TimeLimit)
realScalar2('ObjectiveLimit', options.ObjectiveLimit)

nonNegInteger2('DisplayInterval', options.DisplayInterval)
nonNegInteger2('ReannealInterval', options.ReannealInterval)
nonNegInteger2('PlotInterval', options.PlotInterval)

nonNegScalar2('TolFun',options.TolFun);

% These functions not only savalidate2, they separate fcns from args for
% speed.
[options.AnnealingFcn,options.AnnealingFcnArgs] = functionHandleOrCell2('AnnealingFcn',options.AnnealingFcn);
[options.TemperatureFcn,options.TemperatureFcnArgs] = functionHandleOrCell2('TemperatureFcn',options.TemperatureFcn);
[options.AcceptanceFcn,options.AcceptanceFcnArgs] = functionHandleOrCell2('AcceptanceFcn',options.AcceptanceFcn);

if ~isempty(options.HybridFcn)
    [options.HybridFcn,options.HybridFcnArgs] = functionHandleOrCell2('HybridFcn',options.HybridFcn);
    hybridFcnName = func2str(options.HybridFcn);
    
    unconstrainedHybridFcns = {'fminsearch','fminunc','patternsearch2'};
    constrainedHybridFcns = {'fmincon2','patternsearch2'};
    allHybridFcns = union(unconstrainedHybridFcns,constrainedHybridFcns);
    % Check for a valid hybrid function
    % NOTE: that a hybrid function that is outside of the set above would
    % be caught by one of the checks below. However, the error message
    % would then be slightly misleading. With a dedicated check here, the
    % message is the most clear.
    if ~any(strcmpi(hybridFcnName, allHybridFcns))
        error(message('globaloptim:savalidate2:InvalidHybridFcn', ...
            strjoin(upper(allHybridFcns), ', ') ));
    end
    
    % All HybridFcns only apply to double data types
    if ~strcmpi(options.DataType,'double')
        msg = getString(message('globaloptim:validate2:NonDoubleHybridFcn', ...
            'DataType','double'));
        error('globaloptim:savalidate2:InvalidHybrid',msg);
    end    
    
    % Check for a valid hybrid function for constrained problems
    if problem.bounded && ~any(strcmpi(hybridFcnName,constrainedHybridFcns))
        msg = getString(message('globaloptim:validate2:NotConstrainedHybridFcn', ...
                upper(hybridFcnName),strjoin(upper(constrainedHybridFcns),', ')));
        error('globaloptim:savalidate2:NotConstrainedHybridFcn',msg);
    elseif ~problem.bounded && ~any(strcmpi(hybridFcnName,unconstrainedHybridFcns))
        msg = getString(message('globaloptim:validate2:NotUnconstrainedHybridFcn', ...
                upper(hybridFcnName),strjoin(upper(unconstrainedHybridFcns),', ')));
        error('globaloptim:savalidate2:NotUnconstrainedHybridFcn',msg);
    end
    
    % If the user has set a hybrid function, they can specify options for the
    % hybrid function. If a user has passed a SolverOptions2 object for these
    % options, convert the options object to a structure. Note that we will not
    % warn here if a user specifies a solver with a different solver's options.
    if ~isempty(options.HybridFcnArgs) ...
            && isa(options.HybridFcnArgs{1}, 'optim.options.SolverOptions2')
        % It is possible for a user to pass in a vector of options to the
        % solver. Silently use the first element in this array.
        options.HybridFcnArgs{1} = options.HybridFcnArgs{1}(1);
        
        % Extract the options structure
        options.HybridFcnArgs{1} = extractOptionsStructure(options.HybridFcnArgs{1});
    end
    
end

% this special case takes an array of function cells
[options.PlotFcns,options.PlotFcnsArgs] = functionHandleOrCellArray2('PlotFcns',options.PlotFcns);
[options.OutputFcns,options.OutputFcnsArgs] = functionHandleOrCellArray2('OutputFcns',options.OutputFcns);

% Make the length of the options.InitialTemperature to nvar
len_temp = length(options.InitialTemperature);
if len_temp > problem.nvar
    warning(message('globaloptim:savalidate2:lengthOfInitialTemperature'));
    options.InitialTemperature(len_temp:end) = [];
elseif len_temp < problem.nvar
    temperature = zeros(problem.nvar,1);
    temperature(1:len_temp) = options.InitialTemperature;
    temperature(len_temp+1 : end)= options.InitialTemperature(end);
    options.InitialTemperature = temperature;
end
% We want initialtemperature to be a column vector
options.InitialTemperature = options.InitialTemperature(:);

