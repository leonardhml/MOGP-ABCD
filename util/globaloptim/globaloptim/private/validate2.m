function [options,gLength,fitness,nonlcon] = validate2(options,type,gLength,fitness,nonlcon,user_options)
%VALIDATE2 validates the contents2 of the fitness function, genome length and
%   options struct. 
%   [OUT,nvars,fitness,constr] = VALIDATE2(GenomeLength,FitnessFcn,Nonlcon,IN,type) 
%   validates the FitnessFcn, GenomeLength and the structure IN. OUT is a
%   structure which have all the fields in IN and it gets other fields
%   like FitnessFcn, GenomeLength, etc. The output 'nvars' is the number of
%   variables, 'fitness' is the function handle to the fitness function,
%   and 'nonlcon' is the function handle to the nonlinear constraint
%   function. 
%
%   This function is private to GA2 and GAMULTIOBJ2.

%   Copyright 2003-2015 The MathWorks, Inc.

if nargin < 6
    user_options = options;
end
% Make sure user_options is consistent with gaoptimset2
user_options = gaoptimset2(user_options);

% range check each field
stringSet2('PopulationType',options.PopulationType,{'doubleVector','custom','bitString'});
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
end

validNumberofVariables(gLength)
% PopulationSize validation
if ischar(options.PopulationSize)
    if strcmpi(options.PopulationSize,'15*numberofvariables')
        options.PopulationSize = 15*gLength;
        options.PopulationSize  = floor(options.PopulationSize);
    elseif strcmpi(options.PopulationSize,'50 when numberOfVariables <= 5, else 200')
        if gLength <= 5
            options.PopulationSize = 50;
        else
            options.PopulationSize = 200;
        end
    end
end
positiveIntegerArray2('PopulationSize',options.PopulationSize);
% If population size is a matrix then we want to get the row vector expansion
options.PopulationSize = options.PopulationSize(:)';

% Generations validation
if ischar(options.Generations) 
    if strcmpi(options.Generations,'200*numberofvariables')
        options.Generations = 200*gLength;
        options.Generations = floor(options.Generations);
    elseif strcmpi(options.Generations,'100*numberofvariables')
        options.Generations = 100*gLength;
        options.Generations = floor(options.Generations);
    end
end
positiveInteger2('Generations',options.Generations);

% These options does not apply to gamultiobj2
if ~isempty(options.FitnessLimit) 
    realScalar2('FitnessLimit',options.FitnessLimit);
end
if ~isempty(options.StallTimeLimit)
    positiveScalar2('StallTimeLimit',options.StallTimeLimit);
end
if ~isempty(options.FitnessScalingFcn) 
    [options.FitnessScalingFcn,options.FitnessScalingFcnArgs] = functionHandleOrCell2('FitnessScalingFcn',options.FitnessScalingFcn);
end
% Elite count validation
if ischar(options.EliteCount) && strcmpi(options.EliteCount,'0.05*PopulationSize')
    options.EliteCount = ceil(0.05*mean(options.PopulationSize));
elseif ~isempty(options.EliteCount)
    % Protect against EliteCount greater than PopulationSize.
    if options.EliteCount >= sum(options.PopulationSize)
        error(message('globaloptim:validate2:EliteCountGTPop'));
    end
end

% These fields does not apply to GA2
if ~isempty(options.ParetoFraction)
    realUnitScalar2('ParetoFraction',options.ParetoFraction);
end
if ~isempty(options.DistanceMeasureFcn)
     [options.DistanceMeasureFcn,options.DistanceMeasureFcnArgs] = functionHandleOrCell2('DistanceMeasureFcn',options.DistanceMeasureFcn);
end

stringSet2('Vectorized',options.Vectorized,{'on','off'});
realUnitScalar2('CrossoverFraction',options.CrossoverFraction);
positiveInteger2('MigrationInterval',options.MigrationInterval);
realUnitScalar2('MigrationFraction',options.MigrationFraction);
stringSet2('MigrationDirection',options.MigrationDirection,{'both','forward'});
nonNegScalar2('TolFun',options.TolFun);
nonNegScalar2('TolCon',options.TolCon);
positiveScalar2('TimeLimit',options.TimeLimit);
positiveInteger2('StallGenLimit',options.StallGenLimit);

positiveInteger2('PlotInterval',options.PlotInterval);

if ~isempty(options.UseParallel)
    options.SerialUserFcn = ~validateopts_UseParallel(options.UseParallel,true,true);
else
    options.SerialUserFcn = true;
end

% Creation function for constrained GA2 has different default
% Creation function for linearly constrained GA2 has different default

% Creation function for constrained GA2 has different default
options.CreationFcnArgs = {};
if isempty(user_options.CreationFcn) 
    if ~isempty(nonlcon) && strcmpi(options.NonlinConAlgorithm,'penalty')
        % Make gacreationnonlinearfeasible2 the default for nonlinearly
        % constrained problems using the penalty algorithm.
        options.CreationFcn = @gacreationnonlinearfeasible2;
    elseif strcmp(type,'linearconstraints')
        options.CreationFcn = @gacreationlinearfeasible2;
    end
else
    [options.CreationFcn,options.CreationFcnArgs] = functionHandleOrCell2('CreationFcn',options.CreationFcn);
end

% Crossover function for linearly constrained GA2 has different default
if isempty(user_options.CrossoverFcn) && strcmp(type,'linearconstraints')
    options.CrossoverFcn = @crossoverintermediate2;
    options.CrossoverFcnArgs = {};
else
    [options.CrossoverFcn,options.CrossoverFcnArgs] = functionHandleOrCell2('CrossoverFcn',options.CrossoverFcn);
end

[options.SelectionFcn,options.SelectionFcnArgs] = functionHandleOrCell2('SelectionFcn',options.SelectionFcn);
if options.MultiObjective && ~isempty(user_options.SelectionFcn)
    % Check that the selection function is tournament or custom
    selectionFcnName = func2str(options.SelectionFcn);
    if any(strcmpi(selectionFcnName,optim.options.GamultiobjOptions2.InvalidSelectionFcns))
        error(message('globaloptim:validate2:InvalidMultiObjSelectionFcn'));
    end
end


% Mutation function validation
if isempty(user_options.MutationFcn) && ~strcmp(type,'unconstrained')
    % Mutation function for constrained GA2 has different default
    options.MutationFcn = @mutationadaptfeasible2;
    options.MutationFcnArgs = {};
else
    [options.MutationFcn,options.MutationFcnArgs] = functionHandleOrCell2('MutationFcn',options.MutationFcn);
end

if ~isempty(options.HybridFcn)
    [options.HybridFcn,options.HybridFcnArgs] = functionHandleOrCell2('HybridFcn',options.HybridFcn);
    hybridFcnName = func2str(options.HybridFcn);
    
    unconstrainedHybridFcns = {'fminsearch','fminunc','patternsearch2'};
    constrainedHybridFcns = {'fmincon2','patternsearch2'};
    multiObjHybridFcns = {'fgoalattain'};
    allHybridFcns = [union(unconstrainedHybridFcns,constrainedHybridFcns) multiObjHybridFcns];
   
    % Check for a valid hybrid function
    % NOTE: that a hybrid function that is outside of the set above would
    % be caught by one of the checks below. However, the error message
    % would then be slightly misleading. With a dedicated check here, the
    % message is the most clear.
    stringSet2('HybridFcn',hybridFcnName,allHybridFcns);
    
    % All HybridFcns only apply to double population types
    if ~strcmpi(options.PopulationType,'doubleVector')
        error(message('globaloptim:validate2:NonDoubleHybridFcn', ...
            'PopulationType','doubleVector'));
    end    
    
    % Check for invalid combination of problem-type and Hybrid function
    % e.g. HybridFcn = fminunc, for a constrained problem
    if options.MultiObjective && ~any(strcmpi(hybridFcnName,multiObjHybridFcns))
        % gamultiobj2 hybrid function must be fgoalattain
        error(message('globaloptim:validate2:NotMultiObjHybridFcn', ...
                      upper(hybridFcnName),upper(multiObjHybridFcns{:})));
    elseif ~options.MultiObjective % ga2 hybrid function
        % Check for a valid hybrid function for constrained problems
        % NOTE: we have to check for nonlinear constrained problems with
        % nonlcon since "type" is changed to "subtype" in gacommon2 before
        % calling here.
        hasNonlinearConstraints = ~isempty(nonlcon);
        if (hasNonlinearConstraints || any(strcmp(type,{'linearconstraints','boundconstraints'}))) && ...
            ~any(strcmpi(hybridFcnName,constrainedHybridFcns))
            error(message('globaloptim:validate2:NotConstrainedHybridFcn', ...
                upper(hybridFcnName),strjoin(upper(constrainedHybridFcns),', ')));
        elseif ~hasNonlinearConstraints && strcmp(type,'unconstrained') && ...
               ~any(strcmpi(hybridFcnName,unconstrainedHybridFcns))
            error(message('globaloptim:validate2:NotUnconstrainedHybridFcn', ...
                upper(hybridFcnName),strjoin(upper(unconstrainedHybridFcns),', ')));
        end
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

[options.PlotFcns,options.PlotFcnsArgs] = functionHandleOrCellArray2('PlotFcns',options.PlotFcns);
[options.OutputFcns,options.OutputFcnsArgs] = functionHandleOrCellArray2('OutputFcns',options.OutputFcns);

options.FitnessFcnArgs = {};
options.NonconFcnArgs = {};
if ~isempty(fitness)
    [fitness,FitnessFcnArgs] = functionHandleOrCell2('FitnessFcn',fitness);
    fitness = createAnonymousFcn2(fitness,FitnessFcnArgs);
else
    fitness = [];
end
if ~isempty(nonlcon)
    [nonlcon,NonconFcnArgs] = functionHandleOrCell2('NonconFcn',nonlcon);
    nonlcon = createAnonymousFcn2(nonlcon,NonconFcnArgs);
else
    nonlcon = [];
end

% gacreationnonlinearfeasible2 needs a handle to the nonlinear constraints
if strcmpi(func2str(options.CreationFcn),'gacreationnonlinearfeasible2')
    options.CreationFcnArgs = [{nonlcon} options.CreationFcnArgs];
    
    % Extra validation for gacreationnonlinearfeasible2
    % - This function only works with double data
    if ~strcmpi(options.PopulationType,'doubleVector')
        error(message('globaloptim:validate2:CreateNonlinFeasBadPopType'));
    end
    % - This function is only used for nonlinearly constrained problems
    if isempty(nonlcon)
        if strcmp(type,'linearconstraints')
            options.CreationFcn = @gacreationlinearfeasible2;
        else
            options.CreationFcn = @gacreationuniform2;
        end
        options.CreationFcnArgs = {};
    end
end

% Adjust PopInitRange, if necessary
options = rangeCorrection(gLength,options);

% Additional checks for 'bitString' population type
if strcmpi(options.PopulationType,'bitString') 
    % Verify that if population type is 'bitString' then initial population is
    % not 'logical' data type (common mistake in input)    
    if islogical(options.InitialPopulation)
        options.InitialPopulation = double(options.InitialPopulation);
    end
    % Also make sure that the default CrossoverFcn is set to
    % constraint-preserving crossover function
    if isempty(user_options.CrossoverFcn)
       options.CrossoverFcn = @crossoverscattered2;
       options.CrossoverFcnArgs = {};
    end
    % Warn if crossover function is known not to work with bitString and
    % what are supported ones for bitString.
    crossoverfcn = func2str(options.CrossoverFcn);
    if any(strcmpi(crossoverfcn, {'crossoverintermediate2','crossoverarithmetic2','crossoverheuristic2'}))
       warning(message('globaloptim:validate2:bitStringCrossoverFcn', crossoverfcn));
       options.CrossoverFcn = @crossoverscattered2;
       options.CrossoverFcnArgs = {};       
    end
end

% Remaining tests do not apply to custom population
if strcmpi(options.PopulationType,'custom')
    return;
end

if ~isnumeric(gaoptimget2(options,'InitialPopulation'))
    error(message('globaloptim:validate2:invalidInitialPopulation'));
end
if ~isnumeric(gaoptimget2(options,'InitialScores'))
    error(message('globaloptim:validate2:invalidInitialScores'));
end

% Make sure that initial population is consistent with GenomeLength
if ~isempty(options.InitialPopulation) && size(options.InitialPopulation,2) ~= gLength
    error(message('globaloptim:validate2:wrongSizeInitialPopulation'));
end


%-------------------------------------------------------------------------------

% Number of variables
function validNumberofVariables(GenomeLength)
valid =  isnumeric(GenomeLength) && isscalar(GenomeLength)&& (GenomeLength > 0) ...
         && (GenomeLength == floor(GenomeLength));
if ~valid
   error(message('globaloptim:validate2:validNumberofVariables:notValidNvars'));
end

%------------------------------------------------------------------------------
function options = rangeCorrection(nvars,options)
%rangeCorrection Check the size and consistency of PopInitRange

Range = options.PopInitRange;

lb = Range(1,:);
lb = lb(:);
lenlb = length(lb);
ub = Range(2,:);
ub = ub(:);
lenub = length(ub);

% Check maximum length
if lenlb > nvars
   lb = lb(1:nvars);   
elseif lenlb < nvars
   lb = [lb; lb(end)*ones(nvars-lenlb,1)];
end

if lenub > nvars
   ub = ub(1:nvars);
elseif lenub < nvars
   ub = [ub; ub(end)*ones(nvars-lenub,1)];
end
   
if any( lb > ub )
    count = nnz(lb > ub);
    error(message('globaloptim:validate2:infeasibleRange',count));
end

options.PopInitRange = [lb,ub]';
%------------------------------End of rangeCorrection --------------------------


