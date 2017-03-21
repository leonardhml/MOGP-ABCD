function [xC,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = ...
    lsqncommon2(funfcn,xC,lb,ub,options,defaultopt,caller,initVals, ...
    sizes,flags,mtxmpy,optionFeedback,varargin)
%

%lsqncommon2 Solves non-linear least squares problems.
%   It contains all the setup code common to both LSQNONLIN and lsqcurvefit2 to 
%   call the different algorithms.

%   Copyright 1990-2015 The MathWorks, Inc.

lFinite = ~isinf(lb);
uFinite = ~isinf(ub);

[sizes.fRows,sizes.fCols] = size(initVals.F);
sizes.nFun = numel(initVals.F);
initVals.F = initVals.F(:);

diagnostics = strcmpi(optimget(options,'Diagnostics',defaultopt,'fast'),'on');

algorithm = optimget(options,'Algorithm',defaultopt,'fast');
initDamping = optimget(options,'InitDamping',defaultopt,'fast');
if ~iscell(algorithm)    
    initLMparam = initDamping; 
else
    initLMparam = algorithm{2}; % Initial Levenberg-Marquardt parameter
    algorithm = algorithm{1};   % Algorithm string
end

switch algorithm
    case 'trust2-region-reflective'
        trustRegion = true;
    case 'levenberg-marquardt'
        trustRegion = false;
    otherwise
        error(message('optimlib:lsqncommon2:InvalidAlgorithm', upper( caller )))
end

if flags.grad
    % check size of Jacobian
    [Jrows, Jcols] = size(initVals.J);
    if isempty(options.JacobMult) 
        % Not using 'JacobMult' so Jacobian must be correct size
        if Jrows ~= sizes.nFun || Jcols ~= sizes.nVar
            error(message('optimlib:lsqncommon2:InvalidJacSize', sizes.nFun, sizes.nVar))
        end
    end
else
    Jrows = sizes.nFun; 
    Jcols = sizes.nVar;   
end

% trust2-region-reflective and not enough equations -- switch to other methods
if trustRegion && sizes.nFun < sizes.nVar
    if (isempty(lb(lFinite)) && isempty(ub(uFinite)))
        warning(message('optimlib:lsqncommon2:SwitchToLineSearch'))
    end
    trustRegion = false; % Call levenbergMarquardt2
% Leveberg-Marquardt or line-search with bounds and enough equations,
% switch to trust2-region
elseif ~trustRegion && (~isempty(lb(lFinite)) || ~isempty(ub(uFinite))) && sizes.nFun >= sizes.nVar
    warning(message('optimlib:lsqncommon2:SwitchToLargeScale'))
    trustRegion = true;
end
% Can't handle this one:
if (~isempty(lb(lFinite)) || ~isempty(ub(uFinite))) && sizes.nFun < sizes.nVar
    error(message('optimlib:lsqncommon2:ProblemNotHandled'));
end

% Set confcn for diagnostics and derivative check
confcn = {''};

if diagnostics
    % Do diagnostics on information so far
    non_eq=0; non_ineq=0; lin_eq=0; lin_ineq=0;
    hessflag=false; constflag=false; gradconstflag=false;
    if trustRegion
        OUTPUT.algorithm = 'trust2-region-reflective';
    else
        OUTPUT.algorithm = 'levenberg-marquardt';
    end
    diagnose2(caller,OUTPUT,flags.grad,hessflag,constflag,gradconstflag,...
        xC(:),non_eq,non_ineq,lin_eq,lin_ineq,lb,ub,funfcn,confcn);
end

% Prepare options and flags for finitedifferences2
sizes.xShape = size(xC);
sizes.mNonlinEq = 0; sizes.mNonlinIneq = 0; % No nonlinear constraints

% Read options for finitedifferences2
options.FinDiffType = optimget(options,'FinDiffType',defaultopt,'fast');
options.GradObj = optimget(options,'Jacobian',defaultopt,'fast');
options.GradConstr = 'off';
DerivativeCheck = strcmp(optimget(options,'DerivativeCheck',defaultopt,'fast'),'on');
% Read in and error check option TypicalX
[typicalx,ME] = getNumericOrStringFieldValue2('TypicalX','ones(numberOfVariables,1)', ...
    ones(sizes.nVar,1),'a numeric value',options,defaultopt);
if ~isempty(ME)
    throw(ME)
end
checkoptionsize2('TypicalX', size(typicalx), sizes.nVar);
options.TypicalX = typicalx;
options.DiffMinChange = optimget(options,'DiffMinChange',defaultopt,'fast');
options.DiffMaxChange = optimget(options,'DiffMaxChange',defaultopt,'fast');
options = validateFinDiffRelStep2(sizes.nVar,options,defaultopt);
options.UseParallel     = optimget(options,'UseParallel',defaultopt,'fast');

% Create structure of flags for finitedifferences2
finDiffFlags.fwdFinDiff = strcmpi(options.FinDiffType,'forward'); % Check for forward fin-diff
finDiffFlags.scaleObjConstr = false; % No scaling
finDiffFlags.chkComplexObj = false;  % Don't check whether objective function values are complex
finDiffFlags.isGrad = false;         % Compute Jacobian, not gradient
finDiffFlags.hasLBs = false(sizes.nVar,1);
finDiffFlags.hasUBs = false(sizes.nVar,1);
% finDiffFlags.chkFunEval will be first used by validateFirstDerivatives2().
if trustRegion
    % Don't check function values during finite differences
    finDiffFlags.chkFunEval = false;
else
    % Check function values during finite differences
    finDiffFlags.chkFunEval = true;
end

if trustRegion  % Check for finite bounds
    if ~isempty(lb)
        finDiffFlags.hasLBs = isfinite(lb);
    end
    if ~isempty(ub)
        finDiffFlags.hasUBs = isfinite(ub);
    end
end

% Check derivatives
if DerivativeCheck && flags.grad          % user wants to check derivatives
    validateFirstDerivatives2(funfcn,confcn,xC, ...
        lb,ub,options,finDiffFlags,sizes,varargin{:});
end

% Execute algorithm
if trustRegion
    if ~flags.grad % provide sparsity of Jacobian if not provided.
        Jstr = optimget(options,'JacobPattern',defaultopt,'fast');
        if isempty(Jstr)  
            % Put this code separate as it might generate OUT OF MEMORY error
            Jstr = sparse(ones(Jrows,Jcols));
        elseif ischar(Jstr) 
            % options.JacobPattern is the default: 'sparse(ones(jrows,jcols))'
            Jstr = sparse(ones(Jrows,Jcols));
        else 
            % Pattern matrix  - other datatypes (cell-array, struct) are checked in optimset and its
            % helper functions
            checkoptionsize2('JacobPattern', size(Jstr), sizes.nVar, sizes.nFun);
        end
    else
        Jstr = [];
    end
    % Set MaxFunEvals appropriately for trust2-region-reflective
    defaultopt.MaxFunEvals = '100*numberOfVariables';
    
    if any(lb==ub)
        [xC,FVAL,LAMBDA,JACOB,EXITFLAG,OUTPUT,msgData]= ...
            snlsFixedVar2(funfcn,xC,lb,ub,flags.verbosity,options,defaultopt,initVals.F,initVals.J,caller, ...
            Jstr,flags.computeLambda,mtxmpy,flags.detailedExitMsg,optionFeedback,finDiffFlags,varargin{:});
    else
        [xC,FVAL,LAMBDA,JACOB,EXITFLAG,OUTPUT,msgData]= ...
            snls2(funfcn,xC,lb,ub,flags.verbosity,options,defaultopt,initVals.F,initVals.J,caller, ...
            Jstr,flags.computeLambda,mtxmpy,flags.detailedExitMsg,optionFeedback,finDiffFlags,varargin{:});
    end
else
    % Set MaxFunEvals appropriately for Levenberg-Marquardt
    defaultopt.MaxFunEvals = '200*numberOfVariables';
    
    [xC,FVAL,JACOB,EXITFLAG,OUTPUT,msgData] = ...
        levenbergMarquardt2(funfcn,xC,flags.verbosity,options,defaultopt,initVals.F,initVals.J, ...
        caller,initLMparam,flags.detailedExitMsg,optionFeedback,finDiffFlags,varargin{:});
    LAMBDA.upper = zeros(sizes.nVar,1); 
    LAMBDA.lower = zeros(sizes.nVar,1);
end

Resnorm = FVAL'*FVAL;
OUTPUT.message = createExitMsg2(msgData{:});
% Reset FVAL to original shapes
FVAL = reshape(FVAL,sizes.fRows,sizes.fCols);


