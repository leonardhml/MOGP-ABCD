function options = saoptimset2(varargin)
%SAOPTIMSET2 Create/alter SIMULANNEALBND2 OPTIONS structure.  
%   OPTIONS = SAOPTIMSET2('PARAM1',VALUE1,'PARAM2',VALUE2, ...) creates an
%   optimization options structure OPTIONS in which the named parameters
%   have the specified values. Any unspecified parameters are set to
%   default value for that parameter. It is sufficient to type  only the
%   leading characters that uniquely identify the parameter. Case is
%   ignored for parameter names. NOTE: For values that are strings, correct
%   case and the complete string are required.
%   
%   OPTIONS = SAOPTIMSET2(OLDOPTS,'PARAM1',VALUE1, ...) creates a copy of
%   OLDOPTS with the named parameters altered with the specified values.
%   
%   OPTIONS = SAOPTIMSET2(OLDOPTS,NEWOPTS) combines an existing options
%   structure OLDOPTS with a new options structure NEWOPTS. Any parameters
%   in NEWOPTS with non-empty values overwrite the corresponding old
%   parameters in OLDOPTS. 
%   
%   SAOPTIMSET2 with no input arguments and no output arguments displays all
%   parameter names and their possible values, with defaults shown in {}.
%
%   OPTIONS = SAOPTIMSET2 (with no input arguments) creates an options
%   structure OPTIONS where all the fields are set to [].
%
%   OPTIONS = SAOPTIMSET2(@simulannealbnd2) creates an options structure with
%   all the parameter names and default values relevant to the function
%   'simulannealbnd2'.
%
%SAOPTIMSET2 PARAMETERS
%   TolFun               - Termination tolerance on function value
%                          [ positive scalar | {1e-6} ] 
%   MaxIter              - Maximum number of iterations allowed 
%                          [ positive scalar | {Inf} ]
%   MaxFunEvals          - Maximum number of function (objective)
% 	                        evaluations allowed 
%                          [ positive scalar | {3000*numberOfVariables} ]
%   TimeLimit            - Total time (in seconds) allowed for optimization
%                          [ positive scalar | {Inf} ]
%
%   ObjectiveLimit       - Minimum objective function value desired 
%                          [ scalar | {-Inf} ]
%   StallIterLimit       - Number of iterations over which average
%                          change in objective function value at current
%                          point is less than options.TolFun 
%                          [ positive scalar | {'500*numberOfVariables'} ]
%
%   DataType             - The type of decision variable 
%                          [ 'custom' | {'double'} ]
%
%   InitialTemperature   - Initial temperature at start
%                          [ positive scalar | {100} ]
%
%   ReannealInterval     - Interval for Reannealing 
%                          [ positive integer | {100} ]
%
%   AnnealingFcn         - Function used to generate new points 
%                          [ function_handle | @annealingboltz2 | 
%                            {@annealingfast2} ]
%
%   TemperatureFcn       - Function used to update temperature schedule 
%                          [ function_handle  | @temperatureboltz2 | 
%                            @temperaturefast2 | {@temperatureexp2} ]
%
%   AcceptanceFcn        - Function used to determine if a new point is
%                          accepted or not
%                          [ function_handle | {@acceptancesa2}
%           
%   HybridFcn            - Automatically run HybridFcn (another 
%                          optimization function) during or at the end of
%                          iterations of the solver
%                          [ @fminsearch | @patternsearch2 | @fminunc |
%                            @fmincon2 | {[]} ]
%
%   HybridInterval       - Interval (if not 'end' or 'never') at which
%                          HybridFcn is called 
%                          [ positive integer | 'never' | {'end'} ]
%
%   Display              - Controls the level of display 
%                          [ 'off' | 'iter' | 'diagnose2' | {'final'} ]
%
%   DisplayInterval      - Interval for iterative display 
%                          [ positive integer | {10} ]
%
%   OutputFcns           - Function(s) gets iterative data and can change
%                          options at run time
%                          [ function handle or cell array of function 
%                            handles | {[]} ]
%   PlotFcns             - Plot function(s) called during iterations 
%                          [ function handle or cell array of function 
%                            handles | @saplotbestf2 | @saplotbestx2 | 
%                            @saplotf2 | @saplotstopping2 | 
%                            @saplottemperature2 | {[]} ]
%   PlotInterval         - Interval at which PlotFcns are called
%                          [ positive integer {1} ]
%
% 	 See also SAOPTIMGET2.

%   Copyright 2006-2015 The MathWorks, Inc.

% Check to see if an optimoptions2 object has been passed in as the first
% input argument.
if nargin > 0 && isa(varargin{1}, 'optim.options.SolverOptions2')
    error('globaloptim:saoptimset2:OptimOptionsFirstInput', ...
        getString(message('globaloptim:optimoptionsmsg:OptimOptionsFirstInput','SAOPTIMSET2')));
end

% Check to see if an optimoptions2 object has been passed in as the second
% input argument.
if nargin > 1 && isa(varargin{2}, 'optim.options.SolverOptions2')
    error('globaloptim:saoptimset2:OptimOptionsSecondInput', ...
        getString(message('globaloptim:optimoptionsmsg:OptimOptionsSecondInput')));
end

if (nargin == 0) && (nargout == 0)
    fprintf('            AnnealingFcn: [ function_handle | @annealingboltz2 | {@annealingfast2} ]\n');
    fprintf('          TemperatureFcn: [ function_handle | @temperatureboltz2 | @temperaturefast2 | \n');
    fprintf('                           {@temperatureexp2} ]\n');
    fprintf('           AcceptanceFcn: [ function_handle | {@acceptancesa2} ]\n');
    fprintf('\n')
    fprintf('                  TolFun: [ non-negative scalar | {1e-6} ]\n');
    fprintf('          StallIterLimit: [ positive integer    | {''500*numberOfVariables''} ]\n');
    fprintf('\n')
    fprintf('             MaxFunEvals: [ positive integer    | {''3000*numberOfVariables''} ]\n');
    fprintf('               TimeLimit: [ positive scalar     | {Inf} ]\n');
    fprintf('                 MaxIter: [ positive integer    | {Inf} ]\n');
    fprintf('            ObjectiveLimit: [ scalar              | {-Inf} ]\n');
    fprintf('\n')
    fprintf('                 Display: [ ''off'' | ''iter'' | ''diagnose2'' | {''final''} ]\n');
    fprintf('         DisplayInterval: [ positive integer | {10} ]\n');
    fprintf('\n')
    fprintf('               HybridFcn: [ function_handle  | @fminsearch | @patternsearch2 | \n');
    fprintf('                            @fminunc | @fmincon2 | {[]} ]\n');
    fprintf('          HybridInterval: [ positive integer | ''never''     | {''end''} ]\n');
    fprintf('\n')
    fprintf('                PlotFcns: [ function_handle  | @saplotbestf2 | @saplotbestx2 | \n');
    fprintf('                            @saplotstopping2 | @saplottemperature2 | @saplotx2 | @saplotf2 | {[]} ]\n');
    fprintf('            PlotInterval: [ positive integer | {1} ]\n');
    fprintf('\n')
    fprintf('              OutputFcns: [ function_handle  | {[]} ]\n');
    fprintf('      InitialTemperature: [ positive scalar  | 100 ]\n');
    fprintf('        ReannealInterval: [ positive integer | {100} ]\n');
    fprintf('                DataType: [ ''custom'' | {''double''} ]\n');
    return; 
end

numberargs = nargin; 

% Return options with default values and return it when called with one
% output argument
options=struct('AnnealingFcn', [], ...
'TemperatureFcn', [], ...
'AcceptanceFcn', [], ...
'TolFun', [], ...
'StallIterLimit', [], ...
'MaxFunEvals', [], ...
'TimeLimit', [], ...
'MaxIter', [], ...
'ObjectiveLimit', [], ...
'Display', [], ...
'DisplayInterval', [], ...
'HybridFcn', [], ...
'HybridInterval', [], ...
'PlotFcns', [], ...
'PlotInterval', [], ...
'OutputFcns', [], ...
'InitialTemperature', [], ...
'ReannealInterval', [], ...
'DataType', []);

% If we pass in a function name then return the defaults.
if (numberargs==1) && (ischar(varargin{1}) || isa(varargin{1},'function_handle') )
    if ischar(varargin{1})
        funcname = lower(varargin{1});
        if ~exist(funcname)
            error(message('globaloptim:saoptimset2:functionNotFound',funcname));
        end
    elseif isa(varargin{1},'function_handle')
        funcname = func2str(varargin{1});
    end
    try 
        optionsfcn = feval(varargin{1},'defaults');
    catch
        error(message('globaloptim:saoptimset2:noDefaultOptions',funcname))
    end
    % To get output, run the rest of psoptimset2 as if called with psoptimset2(options, optionsfcn)
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
    if ischar(arg) % arg is an option name
        break;
    end
    if ~isempty(arg) % [] is a valid options argument
        if ~isa(arg,'struct')
            error(message('globaloptim:saoptimset2:invalidArgument', i));
        end
        argFieldnames = fieldnames(arg);
        
        for j = 1:m
            if any(strcmp(argFieldnames,Names{j,:}))
                val = arg.(Names{j,:});
            else
                val = [];
            end
            if ~isempty(val)
                if ischar(val)
                    val = lower(deblank(val));
                end
                checkfield2(Names{j,:},val);
                options.(Names{j,:}) = val;
            end
        end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(numberargs-i+1,2) ~= 0
    error(message('globaloptim:saoptimset2:invalidArgPair'));
end
expectval = 0;  % start expecting a name, not a value
while i <= numberargs
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error(message('globaloptim:saoptimset2:invalidArgFormat', i));
        end
        
        lowArg = lower(arg);
        
        j = strmatch(lowArg,names);
        if isempty(j) % if no matches
            error(message('globaloptim:saoptimset2:invalidParamName', arg));
        elseif length(j) > 1 % if more than one match
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
                error(message('globaloptim:saoptimset2:AmbiguousParamName',arg,allNames));
            end
        end
        expectval = 1;                      % we expect a value next
        
    else           
        if ischar(arg)
            arg = lower(deblank(arg));
        end
        checkfield2(Names{j,:},arg);
        options.(Names{j,:}) = arg;
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error(message('globaloptim:saoptimset2:invalidParamVal', arg));
end


%-------------------------------------------------
function checkfield2(field,value)
%checkfield2 Check validity of structure field contents2.
%   checkfield2('field',V) checks the contents2 of the specified
%   value V to be valid for the field 'field'. 
%

% empty matrix is always valid
if isempty(value)
    return
end

switch field
    %Functions                        
    case {'MaxIter', 'DisplayInterval', 'ReannealInterval', 'PlotInterval','TimeLimit','ObjectiveLimit','TolFun','InitialTemperature'}
        valid = isa(value,'double');
        if ~valid
            error(message('globaloptim:saoptimset2:checkfield2:NotAPosDouble',field));  
        end
        
    case {'HybridInterval'}
        valid = isa(value,'double') || (ischar(value) && any(strcmpi(value,{'end','never'})));
        if  ~valid
            error(message('globaloptim:saoptimset2:checkfield2:NotAPosIntEndOrNever',field,'end','never'));  
        end

    case {'MaxFunEvals'}
        valid = isa(value,'double') || (ischar(value) && strcmpi(value,'3000*numberofvariables'));
        if ~valid
            error(message('globaloptim:saoptimset2:checkfield2:NotAPosInteger',field));             
        end
        
    case {'StallIterLimit'}
        valid = (isnumeric(value) && (1 == length(value)) && (value > 0) && (value == floor(value))) || strcmpi(value,'500*numberofvariables');
        if ~valid
            error(message('globaloptim:saoptimset2:checkfield2:NotAPosInteger',field));
        end

    case {'Display'}
        valid = ischar(value) && any(strcmpi(value,{'off', 'none', 'final', 'iter', 'diagnose2'}));
        if ~valid
            error(message('globaloptim:saoptimset2:checkfield2:NotADisplayType','OPTIONS',field,'off','final','iter','diagnose2'));
        end

    case {'DataType'}
        valid = ischar(value) && any(strcmpi(value,{'custom', 'double'}));
        if ~valid
            error(message('globaloptim:saoptimset2:checkfield2:NotCustomOrDouble','OPTIONS',field,'custom','double'));
        end

    case {'TemperatureFcn', 'AnnealingFcn', 'AcceptanceFcn','PlotFcns', 'HybridFcn', 'OutputFcns'}
        valid = iscell(value) ||  isa(value,'function_handle');
        if ~valid
            error(message('globaloptim:saoptimset2:checkfield2:NotAFunctionOrCellArray','OPTIONS',field));
        end
end    

