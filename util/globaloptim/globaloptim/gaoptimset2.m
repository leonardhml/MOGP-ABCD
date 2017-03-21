function options = gaoptimset2(varargin)
%GAOPTIMSET2 Create/alter GA2 OPTIONS structure.
%   GAOPTIMSET2 returns a listing of the fields in the options structure as
%   well as valid parameters and the default parameter.
%   
%   OPTIONS = GAOPTIMSET2('PARAM',VALUE) creates a structure with the
%   default parameters used for all PARAM not specified, and will use the
%   passed argument VALUE for the specified PARAM.
%
%   OPTIONS = GAOPTIMSET2('PARAM1',VALUE1,'PARAM2',VALUE2,....) will create a
%   structure with the default parameters used for all fields not specified.
%   Those FIELDS specified will be assigned the corresponding VALUE passed,
%   PARAM and VALUE should be passed as pairs.
%
%   OPTIONS = GAOPTIMSET2(OLDOPTS,'PARAM',VALUE) will create a structure named 
%   OPTIONS.  OPTIONS is created by altering the PARAM specified of OLDOPTS to
%   become the VALUE passed.  
%
%   OPTIONS = GAOPTIMSET2(OLDOPTS,'PARAM1',VALUE1,'PARAM2',VALUE2,...) will
%   reassign those fields in OLDOPTS specified by PARAM1, PARAM2, ... to 
%   VALUE1, VALUE2, ...
%
%GAOPTIMSET2 PARAMETERS
%   
%   PopulationType      - The type of Population being entered
%                       [ 'bitstring' | 'custom' | {'doubleVector'} ]
%   PopInitRange        - Initial range of values a population may have
%                       [ Matrix  | {[-10;10]} ]
%   PopulationSize      - Positive scalar indicating the number of individuals
%                       [ positive scalar ]
%   EliteCount          - Number of best individuals that survive to next 
%                         generation without any change
%                       [ positive scalar | 0.05*PopulationSize ]
%   CrossoverFraction   - The fraction of genes swapped between individuals
%                       [ positive scalar | {0.8} ]
%   ParetoFraction      - The fraction of population on non-dominated front
%                       [ positive scalar | {0.35} ]
%   MigrationDirection  - Direction that fittest individuals from the various
%                         sub-populations may migrate2 to other sub-populations
%                       ['both' | {'forward'}]  
%   MigrationInterval   - The number of generations between the migration of
%                         the fittest individuals to other sub-populations
%                       [ positive scalar | {20} ]
%   MigrationFraction   - Fraction of those individuals scoring the best
%                         that will migrate2
%                       [ positive scalar | {0.2} ]
%   Generations         - Maximum number of generations allowed
%                       [ positive scalar ]
%   TimeLimit           - Maximum time (in seconds) allowed  
%                       [ positive scalar | {Inf} ]
%   FitnessLimit        - Minimum fitness function value desired 
%                       [ scalar | {-Inf} ]
%   StallGenLimit       - Number of generations over which cumulative
%                         change in fitness function value is less than TolFun
%                       [ positive scalar ]
%   StallTest           - Measure used to check for stalling
%                       [ 'geometricWeighted' | {'averageChange'} ]
%   StallTimeLimit      - Maximum time over which change in fitness function
%                         value is less than zero
%                       [ positive scalar | {Inf} ]
%   TolFun              - Termination tolerance on fitness function value
%                       [ positive scalar ]
%   TolCon              - Termination tolerance on constraints
%                       [ positive scalar | {1e-6} ]
%   InitialPopulation   - The initial population used in seeding the GA2
%                         algorithm; can be partial
%                       [ Matrix | {[]} ]
%   InitialScores       - The initial scores used to determine fitness; used
%                         in seeding the GA2 algorithm; can be partial
%                       [ column vector | {[]} ]
%   NonlinConAlgorithm  - The algorithm used to handle nonlinear constraints
%                         within the GA2 algorithm
%                          [ 'penalty' | {'auglag'} ]
%   InitialPenalty      - Initial value of penalty parameter
%                          [ positive scalar | {10} ]
%   PenaltyFactor       - Penalty update parameter
%                          [ positive scalar | {100} ]
%   CreationFcn         - Function used to generate initial population
%                       [ @gacreationlinearfeasible2 | @gacreationuniform2 ]
%   FitnessScalingFcn   - Function used to scale fitness scores.
%                       [ @fitscalingshiftlinear2 | @fitscalingprop2 | @fitscalingtop2 |
%                         {@fitscalingrank2} ]
%   SelectionFcn        - Function used in selecting parents for next generation
%                       [ @selectionremainder2 | @selectionuniform2 | 
%                         @selectionroulette2  |  @selectiontournament2 | 
%                         @selectionstochunif2 ]
%   CrossoverFcn        - Function used to do crossover
%                       [ @crossoverheuristic2 | @crossoverintermediate2 | 
%                         @crossoversinglepoint2 | @crossovertwopoint2 | 
%                         @crossoverarithmetic2 | @crossoverscattered2 ]
%   MutationFcn         - Function used in mutating genes
%                       [ @mutationuniform2 | @mutationadaptfeasible2 |
%                         @mutationgaussian2 ]
%   DistanceMeasureFcn  - Function used to measure average distance of
%                         individuals from their neighbors
%                       [ {@distancecrowding2} ]
%   HybridFcn           - Another optimization function to be used once GA2 
%                         has normally terminated (for whatever reason)
%                       [ @fminsearch | @patternsearch2 | @fminunc | @fmincon2 | {[]} ]
%   Display              - Level of display 
%                       [ 'off' | 'iter' | 'diagnose2' | {'final'} ]
%   OutputFcns          - Function(s) called in every generation. This is more   
%                         general than PlotFcns.
%
%   PlotFcns            - Function(s) used in plotting various quantities 
%                         during simulation
%                       [ @gaplotbestf2 | @gaplotbestindiv2 | @gaplotdistance2 | 
%                         @gaplotexpectation2 | @gaplotgenealogy2 | @gaplotselection2 |
%                         @gaplotrange2 | @gaplotscorediversity2  | @gaplotscores2 | 
%                         @gaplotstopping2 | @gaplotmaxconstr2 | @gaplotrankhist2 |
%                         @gaplotpareto2 | @gaplotspread2 | @gaplotparetodistance2 | 
%                         {[]} ]
%   PlotInterval        - The number of generations between plotting results
%                       [ positive scalar | {1} ]
%   Vectorized           - Objective function is vectorized and it can evaluate
%                         more than one point in one call 
%                       [ 'on' | {'off'} ]
%   UseParallel         - Use PARFOR to evaluate objective and nonlinear 
%                         constraint functions.
%                       [ logical scalar | true | {false} ]


%   Copyright 2003-2015 The MathWorks, Inc.

% Check to see if an optimoptions2 object has been passed in as the first
% input argument.
if nargin > 0 && isa(varargin{1}, 'optim.options.SolverOptions2')
    error('globaloptim:gaoptimset2:OptimOptionsFirstInput', ...
        getString(message('globaloptim:optimoptionsmsg:OptimOptionsFirstInput','GAOPTIMSET2')));
end

% Check to see if an optimoptions2 object has been passed in as the second
% input argument.
if nargin > 1 && isa(varargin{2}, 'optim.options.SolverOptions2')
    error('globaloptim:gaoptimset2:OptimOptionsSecondInput', ...
        getString(message('globaloptim:optimoptionsmsg:OptimOptionsSecondInput')));
end

if (nargin == 0) && (nargout == 0)
    fprintf('          PopulationType: [ ''bitstring''      | ''custom''    | {''doubleVector''} ]\n');
    fprintf('            PopInitRange: [ matrix           | {[-10;10]} ]\n');
    fprintf('          PopulationSize: [ positive scalar ]\n');
    fprintf('              EliteCount: [ positive scalar  | {0.05*PopulationSize} ]\n');
    fprintf('       CrossoverFraction: [ positive scalar  | {0.8} ]\n\n');
    fprintf('          ParetoFraction: [ positive scalar  | {0.35} ]\n\n');
    
    fprintf('      MigrationDirection: [ ''both''           | {''forward''} ]\n');
    fprintf('       MigrationInterval: [ positive scalar  | {20} ]\n');
    fprintf('       MigrationFraction: [ positive scalar  | {0.2} ]\n\n');
    
    fprintf('             Generations: [ positive scalar ]\n');
    fprintf('               TimeLimit: [ positive scalar  | {Inf} ]\n');
    fprintf('            FitnessLimit: [ scalar           | {-Inf} ]\n');
    fprintf('           StallGenLimit: [ positive scalar ]\n');
    fprintf('               StallTest: [ ''geometricWeighted'' | {''averageChange''} ]\n');
    fprintf('          StallTimeLimit: [ positive scalar  | {Inf} ]\n');
    fprintf('                  TolFun: [ positive scalar ]\n\n');
    fprintf('                  TolCon: [ positive scalar  | {1e-6} ]\n\n');
    
    fprintf('       InitialPopulation: [ matrix           | {[]} ]\n');
    fprintf('           InitialScores: [ column vector    | {[]} ]\n\n');
    
    fprintf('      NonlinConAlgorithm: [ ''penalty'' | {''auglag''} ]\n');
    fprintf('          InitialPenalty: [ positive scalar | {10} ]\n');
    fprintf('           PenaltyFactor: [ positive scalar | {100} ]\n\n');

    fprintf('             CreationFcn: [ function_handle  | @gacreationuniform2 | @gacreationlinearfeasible2 ]\n');
    fprintf('       FitnessScalingFcn: [ function_handle  | @fitscalingshiftlinear2  | @fitscalingprop2  | \n');
    fprintf('                            @fitscalingtop2   | {@fitscalingrank2} ]\n');
    fprintf('            SelectionFcn: [ function_handle  | @selectionremainder2    | @selectionuniform2 | \n');
    fprintf('                            @selectionroulette2 | @selectiontournament2   | @selectionstochunif2 ]\n');
    fprintf('            CrossoverFcn: [ function_handle  | @crossoverheuristic2  | @crossoverintermediate2 | \n'); 
    fprintf('                            @crossoversinglepoint2 | @crossovertwopoint2 | @crossoverarithmetic2 | \n');
    fprintf('                            @crossoverscattered2 ]\n');
    fprintf('             MutationFcn: [ function_handle  | @mutationuniform2 | @mutationadaptfeasible2 | \n');
    fprintf('                            @mutationgaussian2 ]\n');
    fprintf('      DistanceMeasureFcn: [ function_handle  | {@distancecrowding2} ]\n');
    fprintf('               HybridFcn: [ @fminsearch | @patternsearch2 | @fminunc | @fmincon2 | {[]} ]\n\n');
    
    fprintf('                 Display: [ ''off'' | ''iter'' | ''diagnose2'' | {''final''} ]\n');
    fprintf('              OutputFcns: [ function_handle  | {[]} ]\n');
    fprintf('                PlotFcns: [ function_handle  | @gaplotbestf2 | @gaplotbestindiv2 | @gaplotdistance2 | \n');
    fprintf('                            @gaplotexpectation2 | @gaplotgenealogy2 | @gaplotselection2 | @gaplotrange2 | \n');
    fprintf('                            @gaplotscorediversity2  | @gaplotscores2 | @gaplotstopping2  | \n'); 
    fprintf('                            @gaplotmaxconstr2 | @gaplotrankhist2 | @gaplotpareto2 | @gaplotspread2 | \n');
    fprintf('                            @gaplotparetodistance2 |{[]} ]\n');
    fprintf('            PlotInterval: [ positive scalar  | {1} ]\n\n');
        
    fprintf('              Vectorized: [ ''on''  | {''off''} ]\n\n');
    fprintf('             UseParallel: [ logical scalar | true | {false} ]\n');
    return; 
end     

numberargs = nargin; 

%Return options with default values and return it when called with one output argument
options=struct('PopulationType', [], ...
               'PopInitRange', [], ...
               'PopulationSize', [], ...
               'EliteCount', [], ...
               'CrossoverFraction', [], ...
               'ParetoFraction', [], ...               
               'MigrationDirection',[], ...
               'MigrationInterval',[], ...
               'MigrationFraction',[], ...
               'Generations', [], ...
               'TimeLimit', [], ...
               'FitnessLimit', [], ...
               'StallGenLimit', [], ...
               'StallTest',[], ...
               'StallTimeLimit', [], ...
               'TolFun', [], ...
               'TolCon', [], ...
               'InitialPopulation',[], ...
               'InitialScores', [], ...
               'NonlinConAlgorithm',[], ...
               'InitialPenalty', [], ...
               'PenaltyFactor', [], ...
               'PlotInterval',[], ...
               'CreationFcn',[], ...
               'FitnessScalingFcn', [], ...
               'SelectionFcn', [], ...
               'CrossoverFcn',[], ...
               'MutationFcn',[], ...
               'DistanceMeasureFcn',[], ...               
               'HybridFcn',[], ...
               'Display', [], ...
               'PlotFcns', [], ...
               'OutputFcns', [], ...
               'Vectorized',[], ...
               'UseParallel', []);   


% If we pass in a function name then return the defaults.
if (numberargs==1) && (ischar(varargin{1}) || isa(varargin{1},'function_handle') )
    if ischar(varargin{1})
        funcname = lower(varargin{1});
        if ~exist(funcname,'file')
            error(message('globaloptim:gaoptimset2:functionNotFound',funcname));
        end
    elseif isa(varargin{1},'function_handle')
        funcname = func2str(varargin{1});
    end
    try 
        optionsfcn = feval(varargin{1},'defaults');
    catch ME
        error(message('globaloptim:gaoptimset2:noDefaultOptions',funcname));
    end
    % To get output, run the rest of psoptimset2 as if called with gaoptimset2(options, optionsfcn)
    varargin{1} = options;
    varargin{2} = optionsfcn;
    numberargs = 2;
end

Names = fieldnames(options);
m = size(Names,1);
names = lower(Names);

i = 1;
while i <= numberargs
    arg = varargin{i};
    if ischar(arg)                         % arg is an option name
        break;
    end
    if ~isempty(arg)                      % [] is a valid options argument
        if ~isa(arg,'struct')
            error(message('globaloptim:gaoptimset2:invalidArgument', i));
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),Names{j,:}))
                val = arg.(Names{j,:});
            else
                val = [];
            end
            if ~isempty(val)
                if ischar(val)
                    val = deblank(val);
                end
                options.(Names{j,:}) = checkfield2(Names{j,:},val);
            end
        end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(numberargs-i+1,2) ~= 0
    error(message('globaloptim:gaoptimset2:invalidArgPair'));
end
expectval = 0;                          % start expecting a name, not a value
while i <= numberargs
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error(message('globaloptim:gaoptimset2:invalidArgFormat', i));
        end
        
        lowArg = lower(arg);
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error(message('globaloptim:gaoptimset2:invalidParamName', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                allNames = ['(' Names{j(1),:}];
                for k = j(2:length(j))'
                    allNames = [allNames ', ' Names{k,:}];
                end
                allNames = sprintf('%s).', allNames);
                error(message('globaloptim:gaoptimset2:ambiguousParamName',arg,allNames));
            end
        end
        expectval = 1;                      % we expect a value next
        
    else           
        if ischar(arg)
            arg = (deblank(arg));
        end
        options.(Names{j,:}) = checkfield2(Names{j,:},arg);
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error(message('globaloptim:gaoptimset2:invalidParamVal', arg));
end


%-------------------------------------------------
function value = checkfield2(field,value)
%checkfield2 Check validity of structure field contents2.
%   checkfield2('field',V) checks the contents2 of the specified
%   value V to be valid for the field 'field'. 
%

% empty matrix is always valid
if isempty(value)
    return
end

switch field
    case {'PopulationType','MigrationDirection'}
        if ~isa(value,'char') 
            error(message('globaloptim:gaoptimset2:checkfield2:NotAString','OPTIONS',field));
        end
        
    case {'FitnessScalingFcn','SelectionFcn','CrossoverFcn','MutationFcn',...
                'CreationFcn','HybridFcn','PlotFcns','OutputFcns','DistanceMeasureFcn'}
        if ~(iscell(value) ||  isa(value,'function_handle'))
            error(message('globaloptim:gaoptimset2:checkfield2:NotAFunctionOrCellArray','OPTIONS',field));
        end
        
    case {'ParetoFraction','CrossoverFraction','MigrationInterval','PlotInterval','TolCon', ...
                'TolFun','MigrationFraction','TimeLimit','StallTimeLimit','FitnessLimit','StallGenLimit'} 
        if ~(isa(value,'double'))
            if ischar(value)
                error(message('globaloptim:gaoptimset2:checkfield2:NotAPosRealNumButString','OPTIONS',field));
            else
                error(message('globaloptim:gaoptimset2:checkfield2:NotAPosRealNum','OPTIONS',field));
            end
        end
    case {'PopInitRange'}
        % PopInitRange must be an array of finite doubles with 2 rows
        if ~(isa(value,'double')) || (size(value,1) ~= 2) || any(~isfinite(value(:)))
            error(message('globaloptim:gaoptimset2:checkfield2:NotARangeType','OPTIONS',field));
        end
    case {'InitialPopulation','InitialScores','InitialPenalty','PenaltyFactor'}
       % The content is checked elsewhere.
    case {'Display'}
        validValues = {'off','none','iter','diagnose2','final'};
        checkValidStrings(field,value,validValues);
    case {'Vectorized'}
        if ~isa(value,'char') || ~any(strcmp(value,{'on','off'}))
            error(message('globaloptim:gaoptimset2:checkfield2:NotOnOrOff','OPTIONS',field,'off','on'));
        end
     case {'Generations'} % integer including inf or default string
        if ~(isscalar(value) && isa(value,'double') && value >= 0) && ...
                ~strcmpi(value, '200*numberOfVariables') && ...
                ~strcmpi(value, '100*numberOfVariables')
            if ischar(value)
                error(message('globaloptim:gaoptimset2:checkfield2:NotAPosNumericScalarButString','OPTIONS',field));                
            else
                error(message('globaloptim:gaoptimset2:checkfield2:NotAPosNumericScalar','OPTIONS',field));  
            end
        end
        
    case {'PopulationSize'} % integer including inf or default string
        if ~(isa(value,'double') && all(value(:) >= 0)) && ...
                ~strcmpi(value,'15*numberOfVariables') && ...
                ~strcmpi(value,'50 when numberOfVariables <= 5, else 200')
            if ischar(value)
                error(message('globaloptim:gaoptimset2:checkfield2:NotAPosNumericButString','OPTIONS',field));
            else
                error(message('globaloptim:gaoptimset2:checkfield2:NotAPosNumeric','OPTIONS',field));
            end
        end
    
    case 'EliteCount'
        if ~(isa(value,'double') && all(value(:) >= 0)) && ...
                ~strcmpi(value,'0.05*PopulationSize') 
            if ischar(value)
                error(message('globaloptim:gaoptimset2:checkfield2:NotAPosNumericButString','OPTIONS',field));
            else
                error(message('globaloptim:gaoptimset2:checkfield2:NotAPosNumeric','OPTIONS',field));
            end
        end
    case 'UseParallel'
        [value,valid] = validateopts_UseParallel(value,false,true);
        if  ~valid
            error(message('globaloptim:gaoptimset2:checkfield2:NotLogicalScalar','OPTIONS',field));
        end
    case 'StallTest'
        validValues = {'averageChange','geometricWeighted'};
        checkValidStrings(field,value,validValues);
    case 'NonlinConAlgorithm'
        validValues = {'auglag','penalty'};
        checkValidStrings(field,value,validValues);
    otherwise
        error(message('globaloptim:gaoptimset2:unknownOptionsField'))
end    

%--------------------------------------------------------------------------
function checkValidStrings(optionName,value,validValues)
    try
        stringSet2(optionName,value,validValues);
    catch ME
        error('globaloptim:gaoptimset2:checkfield2:NotAValidString', ...
            ME.message);
    end
        