function [value,validvalue, errmsg, errid, validfield] = checkfield2(field, value, validStrings)
%

%checkfield2 Check validity of optimization options
%

%   [VALIDVALUE, ERRMSG, ERRID, VALIDFIELD] =
%   OPTIMOPTIONCHECKFIELD('field',V) checks the Contents2 of the specified
%   value V to be valid for the field 'field'.
%
%   [VALIDVALUE, ERRMSG, ERRID, VALIDFIELD] =
%   OPTIMOPTIONCHECKFIELD('field',V, VALIDSTRINGS) checks the Contents2 of
%   the specified value V to be valid for the field 'field'. In the case
%   where V can be a string, VALIDSTRINGS contains the possible strings
%   that V can take. VALIDSTRINGS can be a string or a cell array.

%   Copyright 2012-2015 The MathWorks, Inc.

% To check strings, checkfield2 specifies the known string in lower case.
% Convert validStrings to lower case if specified.
if nargin > 2
    validStrings = lower(validStrings);
end

% Some fields are checked in optimset/checkfield2: Display, MaxFunEvals, MaxIter,
% OutputFcn, TolFun, TolX. Some are checked in both (e.g., MaxFunEvals).
validfield = true;
switch field
    % Function or empty
    case {'OutputFcn','OutputFcns','PlotFcns','PlotFcn','CreationFcn', ...
          'FitnessScalingFcn','SelectionFcn','CrossoverFcn','MutationFcn', 'SearchFcn', ...
          'DistanceMeasureFcn', 'AcceptanceFcn', 'AnnealingFcn', 'TemperatureFcn', ...
          'SearchMethod'}
        if isempty(value)
            validvalue = true;
            errmsg = '';
            errid = '';
        else
            [validvalue, errmsg, errid] = functionOrCellArray(field,value);
        end
    case {'Display'}
        [validvalue, errmsg, errid] = stringsType(field,value,validStrings);
    case {'CheckGradients', 'HonorBounds', 'SpecifyConstraintGradient', ...
          'SpecifyObjectiveGradient', 'UseVectorized', 'UseCompletePoll', ...
          'UseCompleteSearch', 'AccelerateMesh', 'ScaleMesh'}
        [validvalue, errmsg, errid] = logicalType(field,value);
    case {'RelLineSrchBnd'}
        if isempty(value)
            validvalue = true;
            errmsg = '';
            errid = '';
        else
            [validvalue, errmsg, errid] = nonNegReal(field,value);
        end
    case {'TolFun','TolX','TolCon','TolPCG','ActiveConstrTol',...
          'DiffMaxChange','DiffMinChange','MaxTime','TimeLimit', ...
          'TolProjCGAbs', 'TolProjCG','TolGradCon','TolConSQP',...
          'TolGapAbs', 'InitDamping', 'ConstraintTolerance','EliteCount', ...
          'FunctionTolerance', 'OptimalityTolerance', 'StepTolerance', ...
          'AbsoluteGapTolerance', 'MeshTolerance','TolFunValue', 'CacheTol', ...
          'CacheSize', 'MaxMeshSize', 'TolBind', 'TolMesh'}
        % non-negative real scalar
        [validvalue, errmsg, errid] = nonNegReal(field,value);
    case {'TolFunLP', 'LPOptimalityTolerance'}
        % real scalar in the range [1e-10, 1e-1]
        [validvalue, errmsg, errid] = boundedReal(field,value,[1e-10, 1e-1]);
    case {'TolGapRel','RelObjThreshold','MinFractionNeighbors','MinNeighborsFraction', ...
          'ParetoFraction','CrossoverFraction','RelativeGapTolerance', ...
          'ObjectiveImprovementThreshold', 'MigrationFraction'}
        % real scalar in the range [0, 1]
        [validvalue, errmsg, errid] = boundedReal(field,value,[0,1]);
    case {'MeshContractionFactor', 'MeshContraction'}
        % real scalar > 1.0, < Inf
        [validvalue, errmsg, errid] = openRangeReal(field,value,[0.0, 1.0]);        
    case {'MeshExpansionFactor', 'MeshExpansion'}
        % real scalar > 1.0, < Inf
        [validvalue, errmsg, errid] = openRangeReal(field,value,[1.0, Inf]);
    case {'TolInteger','IntegerTolerance'}
        % real scalar in the range [1e-6, 1e-3]
        [validvalue, errmsg, errid] = boundedReal(field,value,[1e-6, 1e-3]);
    case {'ObjectiveLimit','FitnessLimit'}
        [validvalue, errmsg, errid] = realLessThanPlusInf(field,value);
    case {'SwarmSize'}
        [validvalue, errmsg, errid] = boundedInteger(field, value, [2,realmax]);
    case {'LargeScale','DerivativeCheck','Diagnostics','GradConstr','GradObj',...
          'Jacobian','Simplex','NoStopIfFlatInfeas','PhaseOneTotalScaling', ...
          'FunValCheck', 'Vectorized', 'MeshRotate', 'Cache', ...
          'CompletePoll', 'CompleteSearch', 'MeshAccelerator'}
        % off, on
        [validvalue, errmsg, errid] = stringsType(field,value,{'on';'off'});
    case {'PrecondBandWidth','MinAbsMax','GoalsExactAchieve', ...
          'RelLineSrchBndDuration', 'DisplayInterval', 'ReannealInterval', ...
          'RootLPMaxIter', 'RootLPMaxIterations', 'MaxFunEvals', 'MaxFunctionEvaluations', ...
          'MaxProjCGIter', 'MaxSQPIter', 'MaxPCGIter', 'MaxIter', ...
          'MaxIterations',  'EqualityGoalCount', 'Generations', 'MaxGenerations', ...
          'MaxStallGenerations', 'StallGenLimit', 'MigrationInterval', ...
          'PopulationSize', 'AbsoluteMaxObjectiveCount', 'HybridInterval'}
        % integer including inf
        [validvalue, errmsg, errid] = nonNegInteger(field,value);
    case {'StallIterLimit','MaxStallIterations'}
        % non-negative integer excluding inf
        [validvalue, errmsg, errid] = boundedInteger(field,value,[0,realmax]);
    case {'InitialSwarm','InitialSwarmMatrix','InitialPopulationMatrix', ...
          'InitialPopulation', 'InitialScoresMatrix', 'InitialScores'}
        % matrix
        [validvalue, errmsg, errid] = twoDimensionalMatrixType(field,value);
    case {'JacobPattern', 'HessPattern', 'InitialTemperature'}
        % matrix or default string
        [validvalue, errmsg, errid] = matrixType(field,value);
    case {'TypicalX'}
        % matrix or default string
        [validvalue, errmsg, errid] = matrixType(field,value);
        % If an array is given, check for zero values and warn
        if validvalue && isa(value,'double') && any(value(:) == 0)
            error('optimlib:options:checkfield2:zeroInTypicalX', ...
                getString(message('MATLAB:optimoptioncheckfield:zeroInTypicalX')));
        end
    case {'HessMult', 'HessFcn', 'JacobMult', 'HessianFcn', 'HessianMultiplyFcn', 'JacobianMultiplyFcn'}
        % function or empty
        if isempty(value)
            validvalue = true;
            errmsg = '';
            errid = '';
        else
            [validvalue, errmsg, errid] = functionType(field,value);
        end
    case {'HessUpdate'}
        % dfp, bfgs, steepdesc
        [validvalue, errmsg, errid] = stringsType(field,value,{'dfp' ; 'steepdesc';'bfgs'});
    case {'MeritFunction'}
        % singleobj, multiobj
        [validvalue, errmsg, errid] = stringsType(field,value,{'singleobj'; 'multiobj' });
    case {'InitialHessType'}
        % identity, scaled-identity, user-supplied
        [validvalue, errmsg, errid] = stringsType(field,value,{'identity' ; 'scaled-identity'; 'user-supplied'});
    case {'UseParallel'}
        % Logical scalar or specific strings
        [value,validvalue] = validateopts_UseParallel(value,false,true);
        if ~validvalue
          msgid = 'MATLAB:optimoptioncheckfield:NotLogicalScalar';
          errid = 'optimlib:options:checkfield2:NotLogicalScalar';
          errmsg = getString(message(msgid, field));
        else
          errid = '';
          errmsg = '';
        end

    case {'Algorithm'}
        % See options objects for the algorithms that are supported for
        % each solver.
        if iscell(value) && numel(value) == 2 && ...
                strcmpi(value{1}, 'levenberg-marquardt')

            % When setting options via optimset, users can specify the
            % Levenberg-Marquardt parameter, lambda, in the following way:
            % opts = optimset('Algorithm', {'levenberg-marquardt',lambda}). 
            %
            % For optimoptions2, we restrict the 'Algorithm' option to be a
            % string only. As such, we provide a helpful error if a user
            % tries to set the Levenberg-Marquardt parameter via a cell
            % array rather than using InitDamping in optimoptions2.
            validvalue = false;           
           msgObj = message('optimlib:options:checkfield2:levMarqAsCell', ...
                            num2str(value{2}), ...
                            addLink2('Setting the Levenberg-Marquardt parameter', ...
                            'optim', 'helptargets.map', 'lsq_set_initdamping', false));
            errid = 'optimlib:options:checkfield2:levMarqAsCell';
            errmsg = getString(msgObj);

        else
            [validvalue, errmsg, errid] = stringsType(field,value,validStrings);
        end
    case {'AlwaysHonorConstraints'}
        % none, bounds
        [validvalue, errmsg, errid] = ...
            stringsType(field,value,{'none' ; 'bounds'});
    case {'ScaleProblem'}
        % NOTE: ScaleProblem accepts logical values (documented) as well as
        % a set of strings (hidden). For Levenberg-Marquardt, it only
        % accepts the strings. Therefore, we check for logical values
        % separately, IF the solver is Levenberg-Marquardt, which we
        % determine from the set of valid strings passed in. Also note,
        % passing the valid set of strings IS a requirement.
        if ~any(strcmp(validStrings,'jacobian')) && islogical(value)
            validvalue = value;
            errmsg = '';
            errid = '';
            return
        end

        [validvalue, errmsg, errid] = stringsType(field,value,validStrings);
    case {'FinDiffType', 'FiniteDifferenceType'}
        % forward, central
        [validvalue, errmsg, errid] = stringsType(field,value,{'forward' ; 'central'});
    case {'FinDiffRelStep', 'FiniteDifferenceStepSize'}
        % Although this option is documented to be a strictly positive
        % vector, matrices are implicitly supported because linear indexing
        % is used. Therefore, posVectorType is called with a reshaped
        % value.
        value = value(:);
        [validvalue, errmsg, errid] = posVectorType(field, value);
    case {'Hessian'}
        if ~iscell(value)
            % If character string, has to be user-supplied, bfgs, lbfgs,
            % fin-diff-grads, on, off
            if isempty(validStrings)
                validStrings = {'user-supplied' ; 'bfgs'; 'lbfgs'; ...
                    'fin-diff-grads'; 'on' ; 'off'};
            end
            [validvalue, errmsg, errid] = ...
                stringsType(field,value,validStrings);
        else
            % If cell-array, has to be {'lbfgs',positive integer}
            [validvalue, errmsg, errid] = stringPosIntegerCellType(field,value,'lbfgs');
        end
    case {'HessianApproximation'}
        if ~iscell(value)
            % If character string, has to be bfgs, lbfgs, or finite-difference
            [validvalue, errmsg, errid] = ...
                stringsType(field,value,{'bfgs'; 'lbfgs'; 'finite-difference'});
        else
            % If cell-array, has to be {'lbfgs',positive integer}
            [validvalue, errmsg, errid] = stringPosIntegerCellType(field,value,'lbfgs');
        end        
    case {'SubproblemAlgorithm'}
        if ~iscell(value)
            % If character string, has to be 'ldl-factorization' or 'cg',
            [validvalue, errmsg, errid] = ...
                stringsType(field,value,{'ldl-factorization' ; 'cg'; 'factorization'});
        else
                % Either {'ldl-factorization',positive integer} or {'cg',positive integer}
                [validvalue, errmsg, errid] = stringPosRealCellType(field,value,{'ldl-factorization' ; 'cg'; 'direct'});
        end
    case {'InitialHessMatrix', 'InitialSwarmSpan'}
        % strictly positive matrix
        [validvalue, errmsg, errid] = posVectorType(field,value);
    case {'BranchingRule', 'BranchRule', 'Heuristics', 'NodeSelection', ...
          'CutGeneration', 'IntegerPreprocess','IPPreprocess', 'LPPreprocess', ...
          'Preprocess', 'RootLPAlgorithm'}
        [validvalue, errmsg, errid] = stringsType(field,value,validStrings);
    case {'InitBarrierParam', 'InitTrustRegionRadius', 'StallTimeLimit', ...
            'MaxStallTime', 'InitialMeshSize', 'InitialPenalty', 'PenaltyFactor'}
        % positive real
        [validvalue, errmsg, errid] = posReal(field,value);
    case {'PlotInterval', 'MaxNodes'}
        % positive integer
        [validvalue, errmsg, errid] = posInteger(field,value);        
    case {'SelfAdjustment','SelfAdjustmentWeight', ...
          'SocialAdjustment','SocialAdjustmentWeight'}
        % particleswarm
        [validvalue, errmsg, errid] = boundedReal(field,value,[-realmax,realmax]);
    case {'InitialPopulationRange', 'PopInitRange'}
        [validvalue, errmsg, errid] = rangeType(field,value);
    case {'InertiaRange'}
        % particleswarm
        [validvalue, errmsg, errid] = sameSignRange(field,value);
    case {'PresolveOps'}
        [validvalue, errmsg, errid] = nonNegIntegerVector(field,value);
    case {'CutGenMaxIter','CutMaxIterations'}
        % intlinprog
        [validvalue, errmsg, errid] = boundedInteger(field,value,[1, 50]);
    case {'MaxNumFeasPoints','MaxFeasiblePoints','LPMaxIter', ...
          'LPMaxIterations','HeuristicsMaxNodes'}
        % intlinprog
        [validvalue, errmsg, errid] = boundedInteger(field, value, [1, Inf]);
    case {'ObjectiveCutOff'}
        % intlinprog
        [validvalue, errmsg, errid] = realGreaterThanMinusInf(field,value);
    case {'PopulationType','DataType'}
        validValues = {'custom','double'};
        if strcmp(field,'PopulationType')
            validValues = [validValues {'doubleVector','bitString'}];
        end
        [validvalue, errmsg, errid] = stringsType(field,value,validValues);
    case {'NonlinearConstraintAlgorithm', 'NonlinConAlgorithm'}
        validValues = {'auglag','penalty'};
        [validvalue, errmsg, errid] = stringsType(field,value,validValues);
    case {'StallTest'}
        validValues = {'geometricWeighted','averageChange'};
        [validvalue, errmsg, errid] = stringsType(field,value,validValues);        
    case {'MigrationDirection'}
        validValues = {'both','forward'};
        [validvalue, errmsg, errid] = stringsType(field,value,validValues);                
    case {'PollOrderAlgorithm','PollingOrder'}
        validValues = {'random','success','consecutive'};
        [validvalue, errmsg, errid] = stringsType(field,value,validValues);
    case {'PollMethod'}        
        validValues = {'gpspositivebasisnp1', 'gpspositivebasis2n' ...
                'positivebasisnp1', 'positivebasis2n','madspositivebasisnp1', ...
                'madspositivebasis2n', 'gsspositivebasisnp1', 'gsspositivebasis2n'};
        [validvalue, errmsg, errid] = stringsType(field,value,validValues);
    % Function or empty
    case {'HybridFcn'}
        if isempty(value)
            validvalue = true;
            errmsg = '';
            errid = '';
        else
            [validvalue, errmsg, errid] = functionOrCellArray(field,value);
            
            if isempty(errid)
               % Extra checking for set of possible functions (validStrings)
               valueToTest = value;
               if iscell(valueToTest)
                   valueToTest = valueToTest{1};
               end
               if isa(valueToTest,'function_handle')
                   valueToTest = func2str(valueToTest);
               end               
               if ~any(strcmpi(valueToTest,validStrings))
                  % Format strings for error message
                  validStrings = formatCellArrayOfStrings(validStrings);
                  validvalue = false;
                  errid = 'optimlib:options:checkfield2:InvalidHybridFcn';
                  errmsg = getString(message('MATLAB:optimoptioncheckfield:notAStringsType','HybridFcn',validStrings));
               end
            end
        end        
    otherwise
        % External users should not get here. We throw an error to remind
        % internal callers that they need to add new options to this
        % function.
        validfield = false;
        validvalue = false;
        errid = 'optimlib:options:checkfield2:unknownField';
        errmsg = getString(message(errid, field));
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = nonNegReal(field,value,string)
% Any nonnegative real scalar or sometimes a special string
valid =  isreal(value) && isscalar(value) && (value >= 0) ;
if nargin > 2
    valid = valid || isequal(value,string);
end
if ~valid
    if ischar(value)
        msgid = 'MATLAB:optimoptioncheckfield:nonNegRealStringType';
        errid = 'optimlib:options:checkfield2:nonNegRealStringType';
    else
        msgid = 'MATLAB:optimoptioncheckfield:notAnonNegReal';
        errid = 'optimlib:options:checkfield2:notAnonNegReal';
    end
    errmsg = getString(message(msgid, field));
else
    errid = '';
    errmsg = '';
end
%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = nonNegInteger(field,value)
% Any nonnegative real integer scalar or sometimes a special string
valid =  isreal(value) && isscalar(value) && (value >= 0) && value == floor(value) ;
if ~valid
    msgid = 'MATLAB:optimoptioncheckfield:notANonNegInteger';
    errid = 'optimlib:options:checkfield2:notANonNegInteger';
    errmsg = getString(message(msgid, field));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = boundedInteger(field,value,bounds)
% Any positive real integer scalar or sometimes a special string
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    value == floor(value) && (value >= bounds(1)) && (value <= bounds(2));
if ~valid
    errid = 'optimlib:options:checkfield2:notABoundedInteger';
    errmsg = getString(message(errid, field, sprintf('[%6.3g, %6.3g]', bounds(1), bounds(2))));
else
    errid = '';
    errmsg = '';
end

%--------------------------------------------------------------------------------

function [valid, errmsg, errid] = sameSignRange(field,value)
% A two-element vector in ascending order; cannot mix positive and negative
% numbers.
valid = isnumeric(value) && isreal(value) && numel(value) == 2 && ...
    value(1) <= value(2) && (all(value>=0) || all(value<=0));
if ~valid
    errid = 'optimlib:options:checkfield2:notSameSignRange';
    errmsg = getString(message(errid, field));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = twoDimensionalMatrixType(field,value,strings)
% Any matrix
valid =  isa(value,'double') && ismatrix(value);
if nargin > 2
    valid = valid || any(strcmp(value,strings));
end
if ~valid
    if ischar(value)
        errid = 'optimlib:options:checkfield2:twoDimTypeStringType';
    else
        errid = 'optimlib:options:checkfield2:notATwoDimMatrix';
    end
    errmsg = getString(message(errid, field));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = matrixType(field,value)
% Any non-empty double (this "matrix" can have more 2 dimensions)
valid = ~isempty(value) && ismatrix(value) && isa(value,'double');
if ~valid
    msgid = 'MATLAB:optimoptioncheckfield:notAMatrix';
    errid = 'optimlib:options:checkfield2:notAMatrix';
    errmsg = getString(message(msgid, field));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = posVectorType(field,value)
% Any non-empty positive scalar or all positive vector
valid = ~isempty(value) && isa(value,'double') && isvector(value) && all(value > 0) ;
if ~valid
    msgid = 'MATLAB:optimoptioncheckfield:notAPosMatrix';
    errid = 'optimlib:options:checkfield2:notAPosMatrix';
    errmsg = getString(message(msgid, field));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = rangeType(field,value)
% A 2-row, double, all finite, non-empty array
valid = isa(value,'double') && isempty(value) || ...
    (size(value,1) == 2) && all(isfinite(value(:)));
if ~valid
    errid = 'optimlib:options:checkfield2:notARange';
    errmsg = getString(message(errid, field));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = openRangeReal(field,value,range)
% Any scalar
valid = isscalar(value) && isa(value,'double') && ~isempty(value) && ...
    (value > range(1)) && (value < range(2));
if ~valid
    errid = 'optimlib:options:checkfield2:notInAnOpenRangeReal';
    errmsg = getString(message(errid, field, sprintf('%6.3g',range(1)), sprintf('%6.3g',range(2))));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = nonNegIntegerVector(field,value)
% A vector of positive integers
valid = isnumeric(value) && isvector(value) && all(value >= 0) && ...
    all(round(value) - value == 0);
if ~valid
    msgid = 'optimlib:options:checkfield2:notANonNegIntVector';
    errid = 'optimlib:options:checkfield2:notANonNegIntVector';
    errmsg = getString(message(msgid, field));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = logicalType(field,value)
% Any function handle or string (we do not test if the string is a function name)
valid =  isscalar(value) && islogical(value);
if ~valid
    msgid = 'MATLAB:optimoptioncheckfield:NotLogicalScalar';
    errid = 'optimlib:options:checkfield2:NotLogicalScalar';
    errmsg = getString(message(msgid, field));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = functionType(field,value)
% Any function handle or string (we do not test if the string is a function name)
valid =  ischar(value) || isa(value, 'function_handle');
if ~valid
    msgid = 'MATLAB:optimoptioncheckfield:notAFunction';
    errid = 'optimlib:options:checkfield2:notAFunction';
    errmsg = getString(message(msgid, field));
else
    errid = '';
    errmsg = '';
end
%-----------------------------------------------------------------------------------------
function [valid, errmsg, errid] = stringsType(field,value,strings)
% One of the strings in cell array strings
valid =  ischar(value) && any(strcmpi(value,strings));

if ~valid
    % Format strings for error message
    allstrings = formatCellArrayOfStrings(strings);

    msgid = 'MATLAB:optimoptioncheckfield:notAStringsType';
    errid = 'optimlib:options:checkfield2:notAStringsType';
    errmsg = getString(message(msgid, field, allstrings));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------
function [valid, errmsg, errid] = boundedReal(field,value,bounds)
% Scalar in the bounds
valid =  isa(value,'double') && isscalar(value) && ...
    (value >= bounds(1)) && (value <= bounds(2));
if ~valid
    msgid = 'MATLAB:optimoptioncheckfield:notAboundedReal';
    errid = 'optimlib:options:checkfield2:notAboundedReal';
    errmsg = getString(message(msgid, field, sprintf('[%6.3g, %6.3g]', bounds(1), bounds(2))));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------
function [valid, errmsg, errid] = stringPosIntegerCellType(field,value,strings)
% A cell array that is either {strings,positive integer} or {strings}
valid = numel(value) == 1 && any(strcmp(value{1},strings)) || numel(value) == 2 && ...
    any(strcmp(value{1},strings)) && isreal(value{2}) && isscalar(value{2}) && value{2} > 0 && value{2} == floor(value{2});

if ~valid
    msgid = 'MATLAB:optimoptioncheckfield:notAStringPosIntegerCellType';
    errid = 'optimlib:options:checkfield2:notAStringPosIntegerCellType';
    errmsg = getString(message(msgid, field, strings));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------
function [valid, errmsg, errid] = stringPosRealCellType(field,value,strings)
% A cell array that is either {strings,positive real} or {strings}
valid = (numel(value) >= 1) && any(strcmpi(value{1},strings));
if (numel(value) == 2)
   valid = valid && isreal(value{2}) && (value{2} >= 0);
end

if ~valid
    % Format strings for error message
    allstrings = formatCellArrayOfStrings(strings);

    msgid = 'MATLAB:optimoptioncheckfield:notAStringPosRealCellType';
    errid = 'optimlib:options:checkfield2:notAStringPosRealCellType';
    errmsg = getString(message(msgid, field,allstrings));
else
    errid = '';
    errmsg = '';
end
%-----------------------------------------------------------------------------------------
function [valid, errmsg, errid] = posReal(field,value)
% Any positive real scalar or sometimes a special string
valid =  isnumeric(value) && isreal(value) && isscalar(value) && (value > 0) ;
if ~valid
    msgid = 'MATLAB:optimoptioncheckfield:nonPositiveNum';
    errid = 'optimlib:options:checkfield2:nonPositiveNum';
    errmsg = getString(message(msgid, field));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------
function [valid, errmsg, errid] = posInteger(field,value)
% Any positive real scalar or sometimes a special string
valid =  isnumeric(value) && isreal(value) && isscalar(value) && ...
    (value > 0) && value == floor(value);
if ~valid
    errid = 'optimlib:options:checkfield2:nonPositiveInteger';
    errmsg = getString(message(errid, field));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = realLessThanPlusInf(field,value,string)
% Any real scalar that is less than +Inf, or sometimes a special string
valid =  isnumeric(value) && isreal(value) && isscalar(value) && (value < +Inf);
if nargin > 2
    valid = valid || strcmpi(value,string);
end
if ~valid
    if ischar(value)
        msgid = 'MATLAB:optimoptioncheckfield:realLessThanPlusInfStringType';
        errid = 'optimlib:options:checkfield2:realLessThanPlusInfStringType';
    else
        msgid = 'MATLAB:optimoptioncheckfield:PlusInfReal';
        errid = 'optimlib:options:checkfield2:PlusInfReal';
    end
    errmsg = getString(message(msgid, field));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = realGreaterThanMinusInf(field,value)
% Any real scalar that is greater than -Inf
valid =  isnumeric(value) && isreal(value) && isscalar(value) && (value > -Inf);
if ~valid
    errid = 'optimlib:options:checkfield2:minusInfReal';
    errmsg = getString(message(errid, field));
else
    errid = '';
    errmsg = '';
end

%---------------------------------------------------------------------------------
function allstrings = formatCellArrayOfStrings(strings)
%formatCellArrayOfStrings converts cell array of strings "strings" into an
% array of strings "allstrings", with correct punctuation and "or"
% depending on how many strings there are, in order to create readable
% error message.

% To print out the error message beautifully, need to get the commas and
% "or"s in all the correct places while building up the string of possible
% string values.

% Add quotes around each string in the cell array
strings = cellfun(@(x) sprintf('''%s''', x), strings, 'UniformOutput', false);

% Create comma separated list from cell array. Note that strjoin requires
% the cell array to be a 1xN row vector.
allstrings = strjoin(strings(:)', ', ');

% Replace last comma with ', or ' or ' or ' depending on the length of the
% list. If there is only one string then there is no string match and 'or'
% is not inserted into the string.
numStrings = length(strings);
if numStrings > 2
    finalConjunction = ', or ';
elseif numStrings == 2
    finalConjunction = ' or ';
else
    % For one string, there is no comma. The following call to regexprep
    % does nothing in this case. As such, we can set finalConjunction
    % arbitrarily to an empty string.
    finalConjunction = '';
end
allstrings = regexprep(allstrings, ', ', finalConjunction, numStrings-1);

%--------------------------------------------------------------------------------

function [valid, errmsg, errid] = functionOrCellArray(field,value)
% Any function handle, string or cell array of functions
valid = ischar(value) || isa(value, 'function_handle') || iscell(value);
if ~valid
    msgid = 'MATLAB:optimset:notAFunctionOrCellArray';
    errid = 'optimlib:options:checkfield2:notAFunctionOrCellArray';
    errmsg = getString(message(msgid, field));
else
    errid = '';
    errmsg = '';
end
%---------------------------------------------------------------------------------

