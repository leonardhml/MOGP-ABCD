function [x,fval,exitFlag,output,population,scores] = ga2(fun,nvars,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,intcon,options)
%GA2    Constrained optimization using genetic algorithm.
%   GA2 attempts to solve problems of the following forms:
%       min F(X)  subject to:  A*X  <= B, Aeq*X  = Beq (linear constraints)
%        X                     C(X) <= 0, Ceq(X) = 0 (nonlinear constraints)
%                              LB <= X <= UB 
%                              X(i) integer, where i is in the index
%                              vector INTCON (integer constraints)
%
%   Note: If INTCON is not empty, then no equality constraints are allowed.
%   That is:-
%   * Aeq and Beq must be empty
%   * Ceq returned from NONLCON must be empty
%
%   X = GA2(FITNESSFCN,NVARS) finds a local unconstrained minimum X to the
%   FITNESSFCN using GA2. NVARS is the dimension (number of design
%   variables) of the FITNESSFCN. FITNESSFCN accepts a vector X of size
%   1-by-NVARS, and returns a scalar evaluated at X.
%
%   X = GA2(FITNESSFCN,NVARS,A,b) finds a local minimum X to the function
%   FITNESSFCN, subject to the linear inequalities A*X <= B. Linear
%   constraints are not satisfied when the PopulationType option is set to
%   'bitString' or 'custom'. See the documentation for details.
%
%   X = GA2(FITNESSFCN,NVARS,A,b,Aeq,beq) finds a local minimum X to the
%   function FITNESSFCN, subject to the linear equalities Aeq*X = beq as
%   well as A*X <= B. (Set A=[] and B=[] if no inequalities exist.) Linear
%   constraints are not satisfied when the PopulationType option is set to
%   'bitString' or 'custom'. See the documentation for details.
%
%   X = GA2(FITNESSFCN,NVARS,A,b,Aeq,beq,lb,ub) defines a set of lower and
%   upper bounds on the design variables, X, so that a solution is found in
%   the range lb <= X <= ub. Use empty matrices for lb and ub if no bounds
%   exist. Set lb(i) = -Inf if X(i) is unbounded below;  set ub(i) = Inf if
%   X(i) is unbounded above. Linear constraints are not satisfied when the
%   PopulationType option is set to 'bitString' or 'custom'. See the 
%   documentation for details.
%
%   X = GA2(FITNESSFCN,NVARS,A,b,Aeq,beq,lb,ub,NONLCON) subjects the
%   minimization to the constraints defined in NONLCON. The function
%   NONLCON accepts X and returns the vectors C and Ceq, representing the
%   nonlinear inequalities and equalities respectively. GA2 minimizes
%   FITNESSFCN such that C(X)<=0 and Ceq(X)=0. (Set lb=[] and/or ub=[] if
%   no bounds exist.) Nonlinear constraints are not satisfied when the
%   PopulationType option is set to 'bitString' or 'custom'. See the 
%   documentation for details.
%
%   X = GA2(FITNESSFCN,NVARS,A,b,Aeq,beq,lb,ub,NONLCON,options) minimizes
%   with the default optimization parameters replaced by values in OPTIONS.
%   OPTIONS can be created with the optimoptions2 function. See optimoptions2
%   for details. For a list of options accepted by GA2 refer to the
%   documentation.
%
%   X = GA2(FITNESSFCN,NVARS,A,b,[],[],lb,ub,NONLCON,INTCON) requires that
%   the variables listed in INTCON take integer values. Note that GA2 does
%   not solve problems with integer and equality constraints. Pass empty
%   matrices for the Aeq and beq inputs if INTCON is not empty.
%
%   X = GA2(FITNESSFCN,NVARS,A,b,[],[],lb,ub,NONLCON,INTCON,options)
%   minimizes with integer constraints and the default optimization
%   parameters replaced by values in OPTIONS. OPTIONS can be created with
%   the optimoptions2 function. See optimoptions2 for details.
%
%   X = GA2(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a structure
%   that has the following fields:
%       fitnessfcn: <Fitness function>
%            nvars: <Number of design variables>
%            Aineq: <A matrix for inequality constraints>
%            bineq: <b vector for inequality constraints>
%              Aeq: <Aeq matrix for equality constraints>
%              beq: <beq vector for equality constraints>
%               lb: <Lower bound on X>
%               ub: <Upper bound on X>
%          nonlcon: <Nonlinear constraint function>
%           intcon: <Index vector for integer variables>
%          options: <Options created with optimoptions2('ga2',...)>
%         rngstate: <State of the random number generator>
%
%   [X,FVAL] = GA2(FITNESSFCN, ...) returns FVAL, the value of the fitness
%   function FITNESSFCN at the solution X.
%
%   [X,FVAL,EXITFLAG] = GA2(FITNESSFCN, ...) returns EXITFLAG which
%   describes the exit condition of GA2. Possible values of EXITFLAG and the
%   corresponding exit conditions are
%
%     1 Average change in value of the fitness function over
%        options.MaxStallGenerations generations less than 
%        options.FunctionTolerance and constraint violation less than 
%        options.ConstraintTolerance.
%     3 The value of the fitness function did not change in
%        options.MaxStallGenerations generations and constraint violation 
%        less than options.ConstraintTolerance.
%     4 Magnitude of step smaller than machine precision and constraint
%        violation less than options.ConstraintTolerance. This exit 
%        condition applies only to nonlinear constraints.
%     5 Fitness limit reached and constraint violation less than
%        options.ConstraintTolerance. 
%     0 Maximum number of generations exceeded.
%    -1 Optimization terminated by the output or plot function.
%    -2 No feasible point found.
%    -4 Stall time limit exceeded.
%    -5 Time limit exceeded.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = GA2(FITNESSFCN, ...) returns a
%   structure OUTPUT with the following information:
%             rngstate: <State of the random number generator before GA2 started>
%          generations: <Total generations, excluding HybridFcn iterations>
%            funccount: <Total function evaluations>
%        maxconstraint: <Maximum constraint violation>, if any
%              message: <GA2 termination message>
%
%   [X,FVAL,EXITFLAG,OUTPUT,POPULATION] = GA2(FITNESSFCN, ...) returns the
%   final POPULATION at termination.
%
%   [X,FVAL,EXITFLAG,OUTPUT,POPULATION,SCORES] = GA2(FITNESSFCN, ...) returns
%   the SCORES of the final POPULATION.
%
%
%   Example:
%     Unconstrained minimization of 'rastriginsfcn2' fitness function of
%     numberOfVariables = 2
%      x = ga2(@rastriginsfcn2,2)
%
%     Display plotting functions while GA2 minimizes
%      options = optimoptions2('ga2','PlotFcn',...
%        {@gaplotbestf2,@gaplotbestindiv2,@gaplotexpectation2,@gaplotstopping2});
%      [x,fval,exitflag,output] = ga2(@rastriginsfcn2,2,[],[],[],[],[],[],[],options)
%
%   An example with inequality constraints and lower bounds
%    A = [1 1; -1 2; 2 1];  b = [2; 2; 3];  lb = zeros(2,1);
%    % Use mutation function which can handle constraints
%    options = optimoptions2('ga2','MutationFcn',@mutationadaptfeasible2);
%    [x,fval,exitflag] = ga2(@lincontest62,2,A,b,[],[],lb,[],[],options);
%
%     FITNESSFCN can also be an anonymous function:
%        x = ga2(@(x) 3*sin(x(1))+exp(x(2)),2)
%
%   If FITNESSFCN or NONLCON are parameterized, you can use anonymous
%   functions to capture the problem-dependent parameters. Suppose you want
%   to minimize the fitness given in the function myfit, subject to the
%   nonlinear constraint myconstr, where these two functions are
%   parameterized by their second argument a1 and a2, respectively. Here
%   myfit and myconstr are MATLAB file functions such as
%
%        function f = myfit(x,a1)
%        f = exp(x(1))*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 2*x(2) + a1);
%
%   and
%
%        function [c,ceq] = myconstr(x,a2)
%        c = [1.5 + x(1)*x(2) - x(1) - x(2);
%              -x(1)*x(2) - a2];
%        % No nonlinear equality constraints:
%         ceq = [];
%
%   To optimize for specific values of a1 and a2, first assign the values
%   to these two parameters. Then create two one-argument anonymous
%   functions that capture the values of a1 and a2, and call myfit and
%   myconstr with two arguments. Finally, pass these anonymous functions to
%   GA2:
%
%     a1 = 1; a2 = 10; % define parameters first
%     % Mutation function for constrained minimization
%     options = optimoptions2('ga2','MutationFcn',@mutationadaptfeasible2);
%     x = ga2(@(x)myfit(x,a1),2,[],[],[],[],[],[],@(x)myconstr(x,a2),options)
%
%   Example: Solving a mixed-integer optimization problem
%   An example of optimizing a function where a subset of the variables are
%   required to be integers:
%   
%   % Define the objective and call GA2. Here variables x(2) and x(3) will
%   % be integer. 
%   fun = @(x) (x(1) - 0.2)^2 + (x(2) - 1.7)^2 + (x(3) -5.1)^2; 
%   x = ga2(fun,3,[],[],[],[],[],[],[],[2 3])
%          
%   See also optimoptions2, FITNESSFUNCTION, GAOUTPUTFCNTEMPLATE2, PATTERNSEARCH2, @.

%   Copyright 2003-2015 The MathWorks, Inc.

% If the first arg is not a gaoptimset2, then it's a fitness function followed by a genome
% length. Here we make a gaoptimset2 from the args.
defaultopt = struct('PopulationType', 'doubleVector', ...
    'PopInitRange', [], ... 
    'PopulationSize', '50 when numberOfVariables <= 5, else 200', ... 
    'EliteCount', '0.05*PopulationSize', ...  
    'CrossoverFraction', 0.8, ...
    'MigrationDirection','forward', ...
    'MigrationInterval',20, ...
    'MigrationFraction',0.2, ...
    'Generations', '100*numberOfVariables', ...
    'TimeLimit', inf, ...
    'FitnessLimit', -inf, ...
    'StallTest', 'averageChange', ... 
    'StallGenLimit', 50, ...
    'StallTimeLimit', inf, ...
    'TolFun', 1e-6, ...
    'TolCon', 1e-3, ...
    'InitialPopulation',[], ...
    'InitialScores', [], ...
    'NonlinConAlgorithm', 'auglag', ...    
    'InitialPenalty', 10, ...
    'PenaltyFactor', 100, ...
    'PlotInterval',1, ...
    'CreationFcn',@gacreationuniform2, ...
    'FitnessScalingFcn', @fitscalingrank2, ...
    'SelectionFcn', @selectionstochunif2, ... 
    'CrossoverFcn',@crossoverscattered2, ...
    'MutationFcn',{{@mutationgaussian2 1  1}}, ...
    'HybridFcn',[], ...
    'Display', 'final', ...
    'PlotFcns', [], ...
    'OutputFcns', [], ...
    'Vectorized','off', ...
    'UseParallel', false);

% Check number of input arguments
try 
    narginchk(1,11);
catch ME
    error(message('globaloptim:ga2:numberOfInputs', ME.message));
end

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && isequal(fun,'defaults')
    x = defaultopt;
    return
end

if nargin < 11, options = [];
    if nargin < 10,  intcon = [];
        if nargin < 9,  nonlcon = [];
            if nargin < 8, ub = [];
                if nargin < 7, lb = [];
                    if nargin <6, beq = [];
                        if nargin <5, Aeq = [];
                            if nargin < 4, bineq = [];
                                if nargin < 3, Aineq = [];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% Is third argument a structure
if nargin == 3 && (isstruct(Aineq) || isa(Aineq, 'optim.options.SolverOptions2')) % Old syntax
    options = Aineq; Aineq = [];
end

% Is tenth argument a structure? If so, integer constraints have not been
% specified
if nargin == 10 && (isstruct(intcon) || isa(intcon, 'optim.options.SolverOptions2'))
    options = intcon;
    intcon = [];
end

% One input argument is for problem structure
if nargin == 1
    if isa(fun,'struct')
        [fun,nvars,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,intcon,rngstate,options] = separateOptimStruct(fun);
        % Reset the random number generators
        resetDfltRng2(rngstate);
    else % Single input and non-structure. 
        error(message('globaloptim:ga2:invalidStructInput'));
    end
end

% Prepare the options for the solver
options = prepareOptionsForSolver2(options, 'ga2');

% If fun is a cell array with additional arguments get the function handle
if iscell(fun)
    FitnessFcn = fun{1};
else
    FitnessFcn = fun;
end

% Only function handles or inlines are allowed for FitnessFcn
if isempty(FitnessFcn) ||  ~(isa(FitnessFcn,'inline') || isa(FitnessFcn,'function_handle'))
    error(message('globaloptim:ga2:needFunctionHandle'));
end

% We need to check the nvars here before we call any solver
valid =  isnumeric(nvars) && isscalar(nvars)&& (nvars > 0) ...
    && (nvars == floor(nvars));
if ~valid
    error(message('globaloptim:ga2:notValidNvars'));
end

% Set default PopInitRange for non-MINLP problems
defaultopt.PopInitRange = [-10;10];

% Specific checks and modification of options for mixed integer GA2
if ~isempty(intcon)   
    % Check whether the user has specified options that the mixed integer
    % solver will either ignore or error.
    gaminlpvalidateoptions2(options);    
    % Change the default options for PopulationSize and EliteCount here.
    defaultopt.PopulationSize = max(min(10*nvars, 100), 40);
    defaultopt.EliteCount = ceil(0.05*defaultopt.PopulationSize);
    % Adjust PopInitRange for MINLPs
    defaultopt.PopInitRange = [-1e4 + 1; 1e4 + 1];
end

user_options = options;
% Use default options if empty
if ~isempty(options) && ~isa(options,'struct')
        error(message('globaloptim:ga2:optionsNotAStruct'));
elseif isempty(options)
    options = defaultopt;
end
% Take defaults for parameters that are not in options structure
options = gaoptimset2(defaultopt,options);

% If a user doesn't specify PopInitRange, we want to set it to the
% bounds when we create the initial population. Need to store a flag
% that indicates whether the user has specified PopInitRange so we can
% do this in the creation function.
options.UserSpecPopInitRange = isa(user_options, 'struct') && ...
    isfield(user_options, 'PopInitRange') && ~isempty(user_options.PopInitRange);

% Check for non-double inputs
msg = isoptimargdbl2('GA2', {'NVARS','A',   'b',   'Aeq','beq','lb','ub'}, ...
                            nvars,  Aineq, bineq, Aeq,  beq,  lb,  ub);
if ~isempty(msg)
    error('globaloptim:ga2:dataType',msg);
end

% Introduce field to describe the objective type
options.MultiObjective = false;

[x,fval,exitFlag,output,population,scores,FitnessFcn,nvars,Aineq,bineq,Aeq,beq,lb,ub, ...
    NonconFcn,options,Iterate,type] = gacommon2(nvars,fun,Aineq,bineq,Aeq,beq,lb,ub, ...
                                               nonlcon,intcon,options,user_options);
                                           
if exitFlag < 0
    return;
end

% Turn constraints into right size if they are empty.
if isempty(Aineq)
    Aineq = zeros(0,nvars);
end
if isempty(bineq)
    bineq = zeros(0,1);
end
if isempty(Aeq)
    Aeq = zeros(0,nvars); 
end
if isempty(beq)
    beq = zeros(0,1);
end

if ~isempty(options.OutputFcns) || ~isempty(options.PlotFcns)
    % For calling an OutputFcn, make an options object (to be updated
    % later) that can be passed in
    options.OutputPlotFcnOptions  = optimoptions2(@ga2);
    options.OutputPlotFcnOptions  = copyForOutputAndPlotFcn(options.OutputPlotFcnOptions,options);
end

% Call appropriate single objective optimization solver
if ~isempty(intcon)   
    [x,fval,exitFlag,output,population,scores] = gaminlp2(FitnessFcn,nvars, ...
        Aineq,bineq,Aeq,beq,lb,ub,NonconFcn,intcon,options,output,Iterate);
else    
    switch (output.problemtype)
        case 'unconstrained'
            [x,fval,exitFlag,output,population,scores] = gaunc2(FitnessFcn,nvars, ...
                options,output,Iterate);
        case {'boundconstraints', 'linearconstraints'}
            [x,fval,exitFlag,output,population,scores] = galincon2(FitnessFcn,nvars, ...
                Aineq,bineq,Aeq,beq,lb,ub,options,output,Iterate);
        case 'nonlinearconstr'
            if strcmpi(options.NonlinConAlgorithm,'auglag')
                [x,fval,exitFlag,output,population,scores] = gacon2(FitnessFcn,nvars, ...
                    Aineq,bineq,Aeq,beq,lb,ub,NonconFcn,options,output,Iterate,type);
            else
                [x,fval,exitFlag,output,population,scores] = gapenalty2(...
                    FitnessFcn,nvars,Aineq,bineq,Aeq,beq,lb,ub, ...
                    NonconFcn,[],options,output,Iterate,type);
            end
    end
end
