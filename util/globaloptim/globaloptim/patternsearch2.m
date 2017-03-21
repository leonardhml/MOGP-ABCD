function [X,FVAL,EXITFLAG,OUTPUT] = patternsearch2(FUN,initialX,Aineq,Bineq,Aeq,Beq,LB,UB,nonlcon,options)
%PATTERNSEARCH2 Constrained optimization using pattern search2.
%   PATTERNSEARCH2 attempts to solve problems of the form:
%       min F(X)  subject to:  A*X  <= B, Aeq*X  = Beq (linear constraints)
%        X                     C(X) <= 0, Ceq(X) = 0 (nonlinear constraints)
%                              LB <= X <= UB
%
%   X = PATTERNSEARCH2(FUN,X0) starts at X0 and finds a local minimum X to
%   the function FUN. FUN accepts input X and returns a scalar function
%   value evaluated at X. X0 may be a scalar or vector.
%
%   X = PATTERNSEARCH2(FUN,X0,A,b) starts at X0 and finds a local minimum X
%   to the function FUN, subject to the linear inequalities A*X <= B.
%
%   X = PATTERNSEARCH2(FUN,X0,A,b,Aeq,beq) starts at X0 and finds a local
%   minimum X to the  function FUN, subject to the linear equalities
%   Aeq*X = Beq as well as A*X <= B. (Set A=[] and B=[] if no inequalities
%   exist.)
%
%   X = PATTERNSEARCH2(FUN,X0,A,b,Aeq,beq,LB,UB) defines a set of lower and
%   upper bounds on the design variables, X, so that a solution is found in
%   the range LB <= X <= UB. Use empty matrices for LB and UB if no bounds
%   exist. Set LB(i) = -Inf if X(i) is unbounded below;  set UB(i) = Inf if
%   X(i) is unbounded above.
%
%   X = PATTERNSEARCH2(FUN,X0,A,b,Aeq,beq,LB,UB,NONLCON) subjects the
%   minimization to the constraints defined in NONLCON. The function
%   NONLCON accepts X and returns the vectors C and Ceq, representing the
%   nonlinear inequalities and equalities respectively. PATTERNSEARCH2
%   minimizes FUN such that C(X)<=0 and Ceq(X)=0. (Set LB=[] and/or UB=[]
%   if no bounds exist.)
%
%   X = PATTERNSEARCH2(FUN,X0,A,b,Aeq,beq,LB,UB,NONLCON,options) minimizes
%   with the default optimization parameters replaced by values in OPTIONS. 
%   OPTIONS can be created with the optimoptions2 function. See optimoptions2
%   for details. For a list of options accepted by PATTERNSEARCH2 refer to
%   the documentation.
%
%   X = PATTERNSEARCH2(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure that has the following fields:
%      objective: <Objective function>
%             x0: <Starting point>
%          Aineq: <A matrix for inequality constraints>
%          bineq: <B vector for inequality constraints>
%            Aeq: <A matrix for equality constraints>
%            beq: <B vector for equality constraints>
%             lb: <Lower bound on X>
%             ub: <Upper bound on X>
%        nonlcon: <Nonlinear constraint function>
%        options: <Options created with optimoptions2('patternsearch2',...)>
%         solver: <solver name 'patternsearch2'>
%       rngstate: <State of the random number generator>
%   This syntax is specially useful if you export a problem from
%   PSEARCHTOOL and use it from the command line to call PATTERNSEARCH2.
%   NOTE: PROBLEM must have all the fields as specified above.
%
%   [X,FVAL] = PATTERNSEARCH2(FUN,X0,...) returns FVAL, the value of the
%   objective function FUN at the solution X.
%
%   [X,FVAL,EXITFLAG] = PATTERNSEARCH2(FUN,X0,...) returns EXITFLAG which
%   describes the exit condition of PATTERNSEARCH2. Possible values of
%   EXITFLAG and the corresponding exit conditions are
%
%     1  Magnitude of mesh size is less than specified tolerance and
%         constraint violation less than options.ConstraintTolerance.
%     2  Change in X less than the specified tolerance and
%         constraint violation less than options.ConstraintTolerance.
%     3  Change in FVAL less than the specified tolerance and
%         constraint violation less than options.ConstraintTolerance.
%     4  Magnitude of step smaller than machine precision and
%         constraint violation less than options.ConstraintTolerance.
%     0  Maximum number of function evaluations or iterations reached.
%    -1  Optimization terminated by the output or plot function.
%    -2  No feasible point found.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = PATTERNSEARCH2(FUN,X0,...) returns a
%   structure OUTPUT with the following information:
%          function: <Objective function>
%       problemtype: <Type of problem> (Unconstrained, Bound constrained or
%                     linear constrained)
%        pollmethod: <Polling technique>
%      searchmethod: <Search2 technique used>, if any
%        iterations: <Total iterations>
%         funccount: <Total function evaluations>
%          meshsize: <Mesh size at X>
%     maxconstraint: <Maximum constraint violation>, if any
%           message: <PATTERNSEARCH2 termination message>
%
%   Examples:
%    FUN can be a function handle (using @)
%      X = patternsearch2(@lincontest62, ...)
%    In this case, F = lincontest62(X) returns the scalar function
%    value F of the  function  evaluated at X.
%
%   An example with inequality constraints and lower bounds
%    A = [1 1; -1 2; 2 1];  b = [2; 2; 3];  lb = zeros(2,1);
%    [X,FVAL,EXITFLAG] = patternsearch2(@lincontest62,[0 0],A,b,[],[],lb);
%
%     FUN can also be an anonymous function:
%        X = patternsearch2(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[0 0])
%     returns X = [0;0].
%
%   If FUN or NONLCON are parameterized, you can use anonymous functions to
%   capture the problem-dependent parameters. Suppose you want to minimize
%   the objective given in the function myobj, subject to the nonlinear
%   constraint myconstr, where these two functions are parameterized by
%   their second argument a1 and a2, respectively. Here myfit and myconstr
%   are MATLAB file functions such as
%
%        function f = myobj(x,a1)
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
%   functions that capture the values of a1 and a2, and call myobj and
%   myconstr with two arguments. Finally, pass these anonymous functions to
%   PATTERNSEARCH2:
%
%     a1 = 1; a2 = 10; % define parameters first
%     options = optimoptions2('patternsearch2','Display','iter'); % Display iterative output
%     x = patternsearch2(@(x)myobj(x,a1),[1;2],[],[],[],[],[],[],@(x)myconstr(x,a2),options)
%
%   See also optimoptions2, GA2, PSOUTPUTFCNTEMPLATE2, SEARCHFCNTEMPLATE2, @.

%   Copyright 2003-2015 The MathWorks, Inc.

% Old syntax X = PATTERNSEARCH2(FUN,X0,A,b,Aeq,beq,LB,UB,OPTIONS) may work
% fine but users are encouraged to update to new syntax of PATTERNSEARCH2
% which takes ninth argument as NONLCON and OPTIONS is passed as tenth
% argument.
defaultopt = struct('TolMesh', 1e-6, ...
    'TolCon', 1e-6, ...
    'TolX', 1e-6 , ...
    'TolFun',1e-6 , ...
    'TolBind',1e-3, ...
    'MaxIter', '100*numberofvariables', ...
    'MaxFunEvals', '2000*numberofvariables', ...
    'TimeLimit', Inf, ...
    'MeshContraction', 0.5, ...
    'MeshExpansion', 2.0, ...
    'MeshAccelerator','off', ...
    'MeshRotate','on', ...
    'InitialMeshSize', 1.0, ...
    'ScaleMesh', 'on', ...
    'MaxMeshSize', inf, ...
    'InitialPenalty', 10, ...
    'PenaltyFactor', 100, ...
    'PollMethod', 'gpspositivebasis2n', ...
    'CompletePoll','off', ...
    'PollingOrder', 'consecutive', ...
    'SearchMethod', [], ...
    'CompleteSearch','off', ...
    'Display', 'final', ...
    'OutputFcns', [], ...
    'PlotFcns',[], ...
    'PlotInterval', 1, ...
    'Cache', 'off', ...
    'CacheSize',1e4, ...
    'CacheTol',eps, ...
    'Vectorized','off', ...
    'UseParallel',false ...
    );

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && isequal(FUN,'defaults')
    X = defaultopt;
    return
end

if nargin < 1
    errmsg = getString(message('MATLAB:narginchk:notEnoughInputs'));
elseif nargin > 10
    errmsg = getString(message('MATLAB:narginchk:tooManyInputs'));
else
    errmsg = '';
end
% At least 1 argument is needed.
if ~isempty(errmsg)
    error(message('globaloptim:patternsearch2:atLeastOneInputArg'));
end

if nargin < 10,  options = [];
    if nargin < 9,  nonlcon = [];
        if nargin < 8, UB = [];
            if nargin < 7, LB = [];
                if nargin <6, Beq = [];
                    if nargin <5, Aeq = [];
                        if nargin < 4, Bineq = [];
                            if nargin <3, Aineq= [];
                            end
                        end
                    end
                end
            end
        end
    end
end

% Ninth argument (nonlcon) could be options structure (old syntax from ver 1.0)
if (isstruct(nonlcon) || isa(nonlcon, 'optim.options.SolverOptions2')) && nargin < 10
    options = nonlcon;
    nonlcon = [];
end

% One input argument is for problem structure
if nargin == 1
    if isa(FUN,'struct')
        [FUN,initialX,Aineq,Bineq,Aeq,Beq,LB,UB,nonlcon,rngstate,options] = separateOptimStruct(FUN);
        % Reset the random number generators
        resetDfltRng2(rngstate);
    else % Single input and non-structure.
        error(message('globaloptim:patternsearch2:inputArg'));
    end
end

% Prepare the options for the solver
options = prepareOptionsForSolver2(options, 'patternsearch2');

% Check for non-double inputs
msg = isoptimargdbl2('PATTERNSEARCH2', {'X0','A',   'b',   'Aeq', 'beq','LB','UB'}, ...
                           initialX,  Aineq, Bineq, Aeq,   Beq,  LB,  UB);
if ~isempty(msg)
    error('globaloptim:patternsearch2:dataType',msg);
end

% If FUN is a cell array with additional arguments, handle them
if iscell(FUN)
    objFcnArg = FUN(2:end);
    FUN = FUN{1};
else
    objFcnArg = {};
end

% Only function_handle or inlines are allowed
if isempty(FUN) ||  ~(isa(FUN,'inline') || isa(FUN,'function_handle'))
    error(message('globaloptim:patternsearch2:needFunctionHandleObj'));
end

% If NONLCON is a cell array with additional arguments, handle them
if iscell(nonlcon)
    conFcnArg = nonlcon(2:end);
    nonlcon = nonlcon{1};
else
    conFcnArg = {};
end
% Constraint function must be a function_handle or inline
if ~isempty(nonlcon) && ~(isa(nonlcon,'inline') || isa(nonlcon,'function_handle'))
    error(message('globaloptim:patternsearch2:needFunctionHandleConstr'));
end

if(~isempty(initialX))
    X = initialX;
    Iterate.x = initialX(:);
    numberOfVariables = length(Iterate.x);
else
    error(message('globaloptim:patternsearch2:initialPoint'));
end

% Determine the 'type' of the problem
if ~isempty(nonlcon)
    type = 'nonlinearconstr';
    % Determine the sub-problem type for the constrained problem (used in ALPS)
    if (~isempty(Aeq)   && ~all(Aeq(:) == 0))   || ~isempty(Beq) || ... 
       (~isempty(Aineq) && ~all(Aineq(:) == 0)) || ~isempty(Bineq)
        subtype = 'linearconstraints';
    elseif (~isempty(LB) && any(isfinite(LB(:)))) || ...
           (~isempty(UB) && any(isfinite(UB(:))))
        subtype = 'boundconstraints';
    else
        subtype = 'unconstrained';
    end
    % If Aeq or Aineq is not empty, then problem has linear constraints.
elseif (~isempty(Aeq)   && ~all(Aeq(:) == 0))   || ~isempty(Beq) || ... 
       (~isempty(Aineq) && ~all(Aineq(:) == 0)) || ~isempty(Bineq)
    type = 'linearconstraints';
    % This condition satisfies bound constraints
elseif (~isempty(LB) && any(isfinite(LB(:)))) || ...
        (~isempty(UB) && any(isfinite(UB(:))))
    type = 'boundconstraints';
    % If all constraints fields are empty then it is unconstrained
else
    type = 'unconstrained';
end

% Initialize output structure
OUTPUT = struct('function',FUN,'problemtype',type,'pollmethod',[], 'maxconstraint', [], ...
    'searchmethod',[],'iterations',0,'funccount',0,'meshsize',[],'rngstate',[]);

% Store the random state
dflt = RandStream.getGlobalStream;
OUTPUT.rngstate = struct('state',{dflt.State}, 'type',{dflt.Type});

% If nonlinear constraints, then subtype is needed to process linear
% constraints (see function preProcessLinearConstr2)
if strcmp(type,'nonlinearconstr')
    type = subtype;
end


% Use default options if empty
if ~isempty(options) && ~isa(options,'struct')
    error(message('globaloptim:patternsearch2:optionsNotAStruct'));
elseif isempty(options)
    options = defaultopt;
end
% Save all user options
user_options = options;

% Validate2 all options
options = checkoptions2(options,defaultopt,numberOfVariables);

% For calling an OutputFcn or PlotFcn, make an options object (to be
% updated later) that can be passed in
if ~isempty(options.OutputFcns) || ~isempty(options.PlotFcns)
    options.OutputPlotFcnOptions = copyForOutputAndPlotFcn(...
        optim.options.PatternsearchOptions2, options);
end

% Bound correction
[LB,UB,msg,EXITFLAG] = checkbound2(LB,UB,numberOfVariables);
if EXITFLAG < 0
    FVAL     =    [];
    X(:)     =  Iterate.x;
    OUTPUT.message = msg;
    if options.Verbosity > 0
        fprintf('%s\n',msg)
    end
    return
end
% Linear constraints correction
[Iterate.x,Aineq,Bineq,Aeq,Beq,LB,UB,msg,EXITFLAG] = ...
    preProcessLinearConstr2(Iterate.x,Aineq,Bineq,Aeq,Beq,LB,UB,options.TolCon, ...
                           numberOfVariables,type,options.Verbosity);
if EXITFLAG < 0
    FVAL     =    [];
    X(:)     =  Iterate.x;
    OUTPUT.message = msg;
    if options.Verbosity > 0
        fprintf('%s\n',msg)
    end
    return
end

% Find the null space of equality constraints, required by polling method GSSDirections2
NullBasisAeq = eqconstrnullspace2(Aeq,numberOfVariables);

% Number of constraints
nineqcstr = size(Aineq,1);
neqcstr   = size(Aeq,1);
ncstr     = nineqcstr + neqcstr;

% Check the objective and constraint functions; Evaluate at the start point
[Iterate,OUTPUT.funccount] = poptimfcnchk2(FUN,nonlcon,initialX,Iterate, ...
    options.Vectorized,objFcnArg,conFcnArg);

% Print some diagnostic information if verbosity > 2
if options.Verbosity > 2
    psdiagnose2(FUN,nonlcon,initialX,nineqcstr,neqcstr,ncstr,user_options);
end
% Call appropriate private solver
switch(OUTPUT.problemtype)
    case 'unconstrained'
        [X,FVAL,EXITFLAG,OUTPUT] = pfminunc2(FUN,objFcnArg,initialX,numberOfVariables,Iterate, ...
            options,defaultopt,OUTPUT);
    case 'boundconstraints'
        [X,FVAL,EXITFLAG,OUTPUT] = pfminbnd2(FUN,objFcnArg,initialX,numberOfVariables,Iterate, ...
            LB,UB,options,defaultopt,OUTPUT);
    case 'linearconstraints'
        [X,FVAL,EXITFLAG,OUTPUT] = pfminlcon2(FUN,objFcnArg,initialX,numberOfVariables,Iterate, ...
            Aineq,Bineq,Aeq,Beq,NullBasisAeq,LB,UB,options,defaultopt,OUTPUT);
    case 'nonlinearconstr'
        [X,FVAL,EXITFLAG,OUTPUT] = pfmincon2(FUN,objFcnArg,initialX,numberOfVariables,Iterate, ...
            Aineq,Bineq,Aeq,Beq,NullBasisAeq,LB,UB,nonlcon,conFcnArg,options,defaultopt,OUTPUT,subtype);
end

