function [xcurr, fvec, LAMBDA, JACOB, EXITFLAG, OUTPUT, msgData] = ...
    snlsFixedVar2(funfcn, xstart, l, u, verb, options, defaultopt, ...
    fval, JACval, caller, Jstr, computeLambda, mtxmpy, ...
    detailedExitMsg, optionFeedback, finDiffFlags, varargin)
%

%snlsFixedVar2  Sparse nonlinear least squares solver with fixed variables
%
%   Locate a local solution to the box-constrained nonlinear least-squares
%   problem:
%
%              min { ||F(x)||^2 :  l <= x <= u }
%
%   where F:R^n -> R^m, m > n, and || || is the 2-norm. Also, at least one
%   element of x(i) is fixed, that is l(i) = x(i) = u(i).
%
%   This function first removes the fixed variables from the problem. The
%   reduced problem is passed to snls2 to be solved. Once the reduced
%   problem solution has been found, the fixed variables are added back to
%   the solution.
%
%   OUTPUTS:-
%      xcurr : x from reduced problem with fixed variables
%
%       fvec : Function value
%
%     LAMBDA : Free variables; LAMBDA from reduced problem
%              Fixed variables - Lower bounds; fvec'*JACOB
%              Fixed variables - Upper bounds; -fvec'*JACOB
%
%      JACOB : Free variables; JACOB from reduced problem
%              Fixed variables; JACOB for fixed variables
%
%   EXITFLAG : }
%     OUTPUT : } Output from snls2
%    msgData : }
%
%   NOTE: We will not add anything about equality constraints to the
%   diagnostic output

%   See also snls2

%   Copyright 2014-2015 The MathWorks, Inc.

% Save shape of xstart so we can reshape for user function and output
sizeX = size(xstart);

% Set up sizes
origNVar = numel(l);

% Index where bounds are equal, and those that are free
idxEqual = l == u;
idxFree = ~idxEqual;

% Immediately return if all bounds are equal
if all(idxEqual)
    [xcurr, fvec, LAMBDA, JACOB, EXITFLAG, OUTPUT, msgData] = ...
        i_cleanUpAllFixedVar(funfcn, xstart, l, u, origNVar, ...
        defaultopt, options, Jstr, finDiffFlags, fval, JACval, ...
        mtxmpy, detailedExitMsg, caller, optionFeedback, verb, ...
        computeLambda, idxEqual, idxFree, varargin{:});
    return
end

% Calculate initial Jacobian if required
gradflag =  strcmp(optimget(options,'Jacobian',defaultopt,'fast'),'on');
if ~gradflag 
    % Initial estimate of the full Jacobian wasn't passed to this function
    % and the Jacobian has not been specified (either directly or via a
    % multiply function). We need to calculate the Jacobian w.r.t the fixed
    % variables using finite differencing.
    [JACval, initJacobFuncCount] = i_sfdnls(xstart, Jstr, origNVar, ...
        fval, funfcn, l, u, options, finDiffFlags, caller, varargin{:});        
else
    % Not used finite differencing to calculate initial Jacobian
    initJacobFuncCount = 0;    
end

% Determine initial gradient of the nonlinear least squares function
initGrad = feval(mtxmpy, JACval, fval(:,1), -1, varargin{:});
initGradFixed = initGrad(idxEqual);

% Indicate whether the user has supplied a Jacobian multiply function
hasJacobMult = isfield(options, 'JacobMult') && ~isempty(options.JacobMult);

% Overwrite JacobMult if specified
origmtxmpy = mtxmpy;
if hasJacobMult
    options.JacobMult = @(J, Y, flag)i_evalJacobMultFcn(J, Y, flag, ...
        options.JacobMult,  origNVar, idxFree);
    mtxmpy = options.JacobMult;
end

% Restrict objective function to "free" variables.
funfcn = i_createObjectiveFcn(funfcn, l, idxFree, hasJacobMult);

% Remove fixed variables from Jacobian pattern
if ~isempty(Jstr)
    Jstr = Jstr(:, idxFree);
end

% Separate Jacobian into fixed and free variables
if ~hasJacobMult
    % Get Jacobian with respect to the fixed variables
    Jfixed = JACval(:, idxEqual);
    
    % Remove fixed variables from initial Jacobian
    JACval = JACval(:, idxFree);
else
    % If here, the user has decided to use a Jacobian multiply function and
    % JACval is Jinfo. Furthermore, the Jacobian is never fully formed so
    % we set Jfixed to be empty.
    Jfixed = [];
end

% Wrap output and plot functions
options = snlsFixedVarWrapOutputAndPlotFcns(options, idxEqual, l, initGradFixed);

% Restrict size dependent options to free variables
if isfield(options, 'TypicalX')
    options.TypicalX = options.TypicalX(idxFree);
end
if isfield(options, 'FinDiffRelStep')
    options.FinDiffRelStep = options.FinDiffRelStep(idxFree);
end

% Call snls2 with free variables only
[xcurr,fvec,LAMBDA,JACOB,EXITFLAG,OUTPUT,msgData] = snls2(funfcn,xstart(idxFree),l(idxFree),u(idxFree),verb,options, ...
    defaultopt,fval,JACval,caller,Jstr,computeLambda,mtxmpy,detailedExitMsg,optionFeedback,finDiffFlags,varargin{:});

% Add initial Jacobian function count to number of function evaluations
OUTPUT.funcCount = OUTPUT.funcCount + initJacobFuncCount;

% Insert fixed variable back into x
xOut = zeros(origNVar, 1);
xOut(idxEqual) = l(idxEqual);
xOut(idxFree) = xcurr;
xcurr = xOut;

% Reshape xcurr
xcurr = reshape(xcurr, sizeX);

% Construct Jacobian
JACOB = i_constructJacobian(JACOB, idxFree, idxEqual, Jfixed, ...
    origNVar, hasJacobMult);

% Compute Lambda if required
if computeLambda
    LAMBDA = i_computeLambda(LAMBDA, origNVar, origmtxmpy, JACOB, ...
        fvec, idxEqual, idxFree, varargin{:});
else
    LAMBDA = [];
end

function funfcn = i_createObjectiveFcn(funfcn, lb, idxFree, hasJacobMult)
%I_CREATEOBJECTIVEFCN Wrap user's objective function
%
%    FUNFCN = I_CREATEOBJECTIVEFCN(FUNFCN, LB, IDXFREE, HASJACOBMULT)
%    creates an anonymous function which wraps the user's objective
%    function. These functions allow the fixed variables to be added back
%    into x before the function is evaluated. One of three functions is
%    returned in funfcn{3}.
%
%    FUNFCN{1} | hasJacobMult | Function in funfcn{3} | Description
%    ----------------------------------------------------------------------
%    'fun'     |   n/a        | i_evalObjEqualBnds    | Function evaluation only
%    ----------------------------------------------------------------------
%    'fungrad' |   true       | i_evalObjAndJacInfo   | Function evaluation plus 
%              |              |                       | passes Jacobian
%              |              |                       | info out without
%              |              |                       | modification
%    ----------------------------------------------------------------------
%    'fungrad' |   false      | i_evalObjAndJac       | Function evaluation plus 
%              |              |                       | evaluation of
%              |              |                       | Jacobian over free
%              |              |                       | variables only.
%    ----------------------------------------------------------------------
%    'fun_then |   true       | i_evalObjEqualBnds    | Function evaluation
%    _grad'    |              |                       |
%              | ............ | ..................... | ...................
%              |              | i_evalJacInfo         | Passes Jacobian
%              |              | (funfcn{4})           | info out without
%              |              |                       | modification
%    ----------------------------------------------------------------------
%    'fun_then |   true       | i_evalObjEqualBnds    | Function evaluation
%    _grad'    |              |                       |
%              | ............ | ..................... | ...................
%              |   false      | i_evalJac (funfcn{4}) | Evaluation of
%              |              |                       | Jacobian over free
%              |              |                       | variables only.
%    ----------------------------------------------------------------------
%
%    NOTE: We will not consider bounds to be equal within some tolerance.
%    We'll just handle the case that causes the issue in snls2. That is,
%    bounds are equal if they satisfy lb==ub.

% TODO: Handle all values of funfcn{1}

% Index where bounds are equal
idxEqual = ~idxFree;

% Number of variables in the problem
Nvars = numel(lb);

% Use an internal function to calculate objective
switch funfcn{1}
    case 'fun'
        funfcn{3} = @(x, varargin)i_evalObjEqualBnds(x, funfcn{3}, idxEqual, Nvars, lb, varargin{:});
    case 'fungrad'
        if hasJacobMult
            funfcn{3} = @(x, varargin)i_evalObjAndJacInfo(x, funfcn{3}, idxEqual, Nvars, lb, varargin{:});
        else
            funfcn{3} = @(x, varargin)i_evalObjAndJac(x, funfcn{3}, idxEqual, Nvars, lb, varargin{:});
        end
    case 'fun_then_grad'
        funfcn{3} = @(x, varargin)i_evalObjEqualBnds(x, funfcn{3}, idxEqual, Nvars, lb, varargin{:});
        if hasJacobMult
            funfcn{4} = @(x, varargin)i_evalJacInfo(x, funfcn{4}, idxEqual, Nvars, lb, varargin{:});
        else
            funfcn{4} = @(x, varargin)i_evalJac(x, funfcn{4}, idxEqual, Nvars, lb, varargin{:});
        end
end

function y = i_evalObjEqualBnds(x, fcn, idx, Nvars, lb, varargin)
%I_EVALOBJEQUALBNDS Evaluate user's objective function
%
%    Y = I_EVALOBJEQUALBNDS(X, FCN, IDXFIXED, NVARS, LB) evaluates the
%    user's objective function, FCN, at the free variable point, X. Before
%    the function is called, the fixed variables are inserted back into X.

% Create x with fixed variables
xin = i_createFullX(x, idx, Nvars, lb);

% Evaluate the objective
y = fcn(xin, varargin{:});

function [y, J] = i_evalObjAndJac(x, fcn, idx, Nvars, lb, varargin)
%I_EVALOBJANDJAC Evaluate user's objective function and Jacobian
%
%    [Y, J] = I_EVALOBJANDJAC(X, FCN, IDXFIXED, NVARS, LB) evaluates the
%    user's objective function, FCN, at the free variable point, X. Before
%    the function is called, the fixed variables are inserted back into X.
%    Before the Jacobian, J, is returned, the fixed variable columns are
%    removed from J.

% Create full x
xFull = i_createFullX(x, idx, Nvars, lb);

% Evaluate objective function
[y, J] = fcn(xFull, varargin{:});

% Remove the fixed variables from the Jacobian
J = J(:, ~idx);

function [y, Jinfo] = i_evalObjAndJacInfo(x, fcn, idx, Nvars, lb, varargin)
%I_EVALOBJANDJACINFO Evaluate user's objective function and Jacobian info
%
%    [Y, J] = I_EVALOBJANDJACINFO(X, FCN, IDXFIXED, NVARS, LB) evaluates
%    the user's objective function, FCN, at the free variable point, X.
%    Before the function is called, the fixed variables are inserted back
%    into X. As we will pass Jinfo straight to the user's multiply
%    function, we do not edit it for the presence of fixed variables here.

% Create full x
xFull = i_createFullX(x, idx, Nvars, lb);

% Evaluate objective function and the Jacobian information for the Jacobian
% multiply function. As we will pass Jinfo straight to the user's multiply
% function, we do not edit it for the presence of fixed variables here.
[y, Jinfo] = fcn(xFull, varargin{:});

function J = i_evalJac(x, fcn, idx, Nvars, lb, varargin)
%I_EVALJAC Evaluate user's Jacobian function
%
%    [Y, J] = I_EVALJAC(X, FCN, IDXFIXED, NVARS, LB) evaluates the user's
%    Jacobian function, FCN, at the free variable point, X. Before the
%    function is called, the fixed variables are inserted back into X.
%    Before the Jacobian, J, is returned, the fixed variable columns are
%    removed from J.

% Create full x
xFull = i_createFullX(x, idx, Nvars, lb);

% Evaluate Jacobian function
J = fcn(xFull, varargin{:});

% Remove the fixed variables from the Jacobian
J = J(:, ~idx);

function Jinfo = i_evalJacInfo(x, fcn, idx, Nvars, lb, varargin)
%I_EVALJACINFO Evaluate user's Jacobian info function
%
%    JINFO = I_EVALJACINFO(X, FCN, IDXFIXED, NVARS, LB) evaluates the
%    user's Jacobian information function, FCN, at the free variable point,
%    X. Before the function is called, the fixed variables are inserted
%    back into X. As we will pass Jinfo straight to the user's multiply
%    function, we do not edit it for the presence of fixed variables here.

% Create full x
xFull = i_createFullX(x, idx, Nvars, lb);

% Evaluate objective function and the Jacobian information for the Jacobian
% multiply function. As we will pass Jinfo straight to the user's multiply
% function, we do not edit it for the presence of fixed variables here.
Jinfo = fcn(xFull, varargin{:});

function xin = i_createFullX(x, idx, Nvars, lb)
%I_CREATEFULLX Insert fixed variables back into X
%
%    XIN = I_CREATEFULLX(X, IDXFIXED, NVARS, LB) creates a 1-by-NVARS
%    vector for all the variables in the original problem. For the free
%    variables, XIN(~IDXFIXED) = X. For the fixed variables, XIN(IDXFIXED)
%    = LB.

xin = zeros(1, Nvars);
xin(idx) = lb(idx);
xin(~idx) = x;

function [xcurr, fvec, LAMBDA, JACOB, EXITFLAG, OUTPUT, msgData] = ...
    i_cleanUpAllFixedVar(funfcn, xcurr, lb, ub, n, defaultopt, options, ...
    Jstr, finDiffFlags, fvec, JACval, mtxmpy, detailedExitMsg, caller, ...
    optionFeedback, verb, computeLambda, idxEqual, idxFree, varargin)
%I_CLEANUPALLFIXEDVAR Create solution for all fixed variable case
%
%   [xcurr, ...] = I_CLEANUPALLFIXEDVAR(FUNFCN, XCURR, ...) creates a
%   solution when the problem contains purely fixed variables. In this case
%   there is no need to run snls2 as we can just set XCURR = LB. The other
%   outputs are calculated based on this solution.

% Set x to the fixed value
xcurr(:) = lb;

% Set exitflag
EXITFLAG = 1;

% Calculate Jacobian
gradflag =  strcmp(optimget(options,'Jacobian',defaultopt,'fast'),'on');
if ~gradflag % use sparse finite differencing
    [JACOB, numFinDiffFunEvals] = i_sfdnls(xcurr, Jstr, n, fvec, ...
        funfcn, lb, ub, options, finDiffFlags, caller, varargin{:});
else % user-supplied computation of J or dnewt
    JACOB = JACval;
    numFinDiffFunEvals = 0;
end

% Number of function evaluations taken is equal to that in finite
% differencing plus the function evaluation performed in lsqnonlin/curvefit
numFunEvals = numFinDiffFunEvals + 1;

% Form output structure
OUTPUT.iterations = 0;
OUTPUT.funcCount = numFunEvals; 
OUTPUT.algorithm = 'trust2-region-reflective';
OUTPUT.firstorderopt = 0;
OUTPUT.cgiterations = 0;
OUTPUT.stepsize = 0;

% Create exit message data cell array
tolFun = optimget(options,'TolFun',defaultopt,'fast');
msgData = {'snls2', 100, verb > 0, detailedExitMsg, caller, ...
    0, optionFeedback.TolFun, tolFun};

% Set Lagrange multipliers
if computeLambda
    LAMBDA = i_computeLambda([], n, mtxmpy, JACOB, fvec, idxEqual, ...
        idxFree, varargin{:});
else
    LAMBDA = [];
end

function w = i_evalJacobMultFcn(Jinfo, Y, flag, JacobMultFcn, Nvars, idxFree)
%I_EVALJACOBMULTFCN Evaluate user's Jacobian multiply function
%
%    W = I_EVALJACOBMULTFCN(JINFO, Y, FLAG, JACOBMULTFCN, NVARS, IDXFREE)
%    calls W = JACOBMULTFUN(JINFO, YM, FLAG) where YM may be modified for
%    the presence of fixed variables in the problem depending on FLAG.
%
%    If FLAG > 0 or FLAG < 0, entries for fixed variables are also be
%    removed from W. That is, 
%      
%    W = JACOBMULTFUN(JINFO, YM, FLAG);
%    W = W(IDXFREE, :) 

% Remove contribution from fixed variables from the multiply
% function
if flag > 0
    % w = J*Y
    
    % Create Y over all the variables. For fixed variables, we set
    % Y[idxFixed, :] = 0 as the variables should make no contribution
    % to the Jacobian vector product.
    nColY = size(Y, 2);
    Yfull = zeros(Nvars, nColY);
    Yfull(idxFree, :) = Y;
    
    % Call user's Jacobian multiply function
    w = JacobMultFcn(Jinfo, Yfull, flag);
    
elseif flag < 0
    % w = J'*Y;
    
    % Call user's Jacobian multiply function
    w = JacobMultFcn(Jinfo, Y, flag);
    
    % J' is nVar-by-nEq
    % Y is nEq-by-nColY
    % => w is nVar-by-nColY
    w = w(idxFree, :);
else
    % w = J'*J*Y;
    
    % Create Y over all the variables. For fixed variables, we set
    % Y[idxFixed] = 0 as the variables should make no contribution
    % to the Jacobian vector product.
    nColY = size(Y, 2);
    Yfull = zeros(Nvars, nColY);
    Yfull(idxFree, :) = Y;
    
    % Call user's Jacobian multiply function
    w = JacobMultFcn(Jinfo, Yfull, flag);
    
    % J' is nVar-by-nEq
    % J is nEq-by-nVar
    % Y is nVar-by-nColY
    % => w is (nVar-by-nEq)*(nEq-by-nVar)*(nVar-by-1) => (nVar-by-nColY)
    w = w(idxFree, :);
    
end

function JACOB = i_constructJacobian(JACOB, idxFree, idxEqual, Jfixed, ...
    origNvar, hasJacobMult)
%I_CONSTRUCTJACOBIAN Compute the Jacobian
%
%    JACOBALL = I_CONSTRUCTJACOBIAN(JACOB, ..., JFIXED, ...) constructs the
%    full Jacobian (i.e. with respect to free and fixed variables). First a
%    sparse NFUN-by-ORIGNVAR Jacobian, JACOBALL, is created. Then the
%    Jacobian with respect to the free variables, JACOB is inserted into
%    JACOBALL. If the Jacobian w.r.t the fixed variables, JFIXED, is
%    available, then that is immediately inserted into JACOBALL. Otherwise
%    the Jacobian w.r.t the fixed variables is estimated using sfdnls2 and
%    inserted into JACOBALL.

% Jacobian not defined for fixed variables and need to be added back in.
% However, if the user has passed a Jacobian multiply function, JACOB is
% returned back in a form specified by the user, so we just leave it.
if ~hasJacobMult
        
    % Number of functions
    nFun = size(JACOB, 1);
    
    % Jacobian over all the variables
    JACOBOut = sparse(nFun, origNvar);
    
    % Insert the Jacobian w.r.t the free variables
    JACOBOut(:, idxFree) = JACOB;
    
    % Insert the Jacobian w.r.t the fixed variables
    JACOBOut(:, idxEqual) = Jfixed;
    
    % Set Jacobian
    JACOB = JACOBOut;
    
end

function LAMBDA = i_computeLambda(LAMBDA, origNVar, mtxmpy, JACOB, ...
    fvec, idxEqual, idxFree, varargin)
%I_COMPUTELAMBDA Compute the Lagrange Multipliers (LMs)
%
%    LAMBDAALL = I_COMPUTELAMBDA(LAMBDAFREE, ...) constructs the full
%    Lagrange multipliers (i.e. with respect to free and fixed variables).
%    For each bound the following steps are performed:
%  
%    * Create a ORIGNVAR-by-1 zeros vector, L
%    * Set the Lagrange Multipliers for the free variables, L(IDXFREE) =
%    LAMBDA.<LOWER/UPPER>(idxFree). Note that if all the variables are
%    fixed, then we use the formula in next step to calculate the LMs for
%    the free variables.
%    * Set the Lagrange multipliers for the fixed variables. In this case, 
%
%    L(IDXFIXED & g >= 0) = (2*JACOB*fvec)[IDXFIXED & g >= 0] for lower bound LAMBDA
%    L(IDXFIXED & g >= 0) = 0 for upper bound LAMBDA
%    L(IDXFIXED & g < 0) = 0 for lower bound LAMBDA
%    L(IDXFIXED & g < 0) = -(2*JACOB*fvec)[IDXFIXED & g < 0] for upper bound LAMBDA
%
%    Note that this matches the definition for the Lagrange multipliers
%    returned in snls2. See the comments in snls2 for more background.

% Lagrange multipliers for fixed variables need to be added back in

% Determine gradient of the nonlinear least squares function with
% respect to the fixed variables. Note the factor of 2 as the gradient of
% the objective function is 2*J'*f.
gradFun = 2*feval(mtxmpy,JACOB,fvec(:,1),-1,varargin{:});

% Initialize Lagrange multipliers
LAMBDAOut.lower = zeros(origNVar,1);
LAMBDAOut.upper = zeros(origNVar,1);

% Lagrange multipliers for free variables
if ~isempty(LAMBDA)
    LAMBDAOut.lower(idxFree) = LAMBDA.lower;
    LAMBDAOut.upper(idxFree) = LAMBDA.upper;
end

% Lagrange multipliers for fixed variables
idxLowerActive = idxEqual & gradFun >= 0;
LAMBDAOut.lower(idxLowerActive) = gradFun(idxLowerActive);
idxUpperActive = idxEqual & gradFun < 0;
LAMBDAOut.upper(idxUpperActive) = -gradFun(idxUpperActive);

% Set lambda
LAMBDA = LAMBDAOut;

function [A, findiffevals] = i_sfdnls(xcurr, Jstr, origNvar, fvec, ...
    funfcn, l, u, options, finDiffFlags, caller, varargin)
%I_SFDNLS Estimate the Jacobian using sparse finite differencing

% Set up sizes structure, required for sfdnls2
sizes.xShape = size(xcurr);

% Determine coloring/grouping for sparse finite-differencing
p = colamd(Jstr)';
p = (origNvar+1)*ones(origNvar,1)-p;
group = color2(Jstr,p);

% Pass in user's objective function and not the wrapped function as
% we want the Jacobian with respect to all the variables
[A,findiffevals] = sfdnls2(xcurr,fvec,Jstr,group,[],funfcn,l,u,...
    options,sizes,finDiffFlags,varargin{:});

% Error if Jacobian is undefined
if any(~isfinite(nonzeros(A)))
    error('optimlib:snlsFixedVar2:DerivUndefAtX0', ...
        getString(message('optimlib:commonMsgs:FinDiffJacUndefAtX0', ...
        caller)));
end

% Return sparse Jacobian
A = sparse(A);
