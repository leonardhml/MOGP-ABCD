 function [xCurrent,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = lsqcurvefit2(FUN,xCurrent,XDATA,YDATA,LB,UB,options,varargin)
%lsqcurvefit2 solves non-linear least squares problems.
%   lsqcurvefit2 attempts to solve problems of the form:
%   min  sum {(FUN(X,XDATA)-YDATA).^2}  where X, XDATA, YDATA and the
%    X                                  values returned by FUN can be 
%                                       vectors or matrices.
%
%   lsqcurvefit2 implements two different algorithms: trust2 region reflective and 
%   Levenberg-Marquardt. Choose one via the option Algorithm: for instance, to 
%   choose Levenberg-Marquardt, set 
%   OPTIONS = optimoptions2('lsqcurvefit2','Algorithm','levenberg-marquardt'), 
%   and then pass OPTIONS to lsqcurvefit2.  
%
%   X = lsqcurvefit2(FUN,X0,XDATA,YDATA) starts at X0 and finds coefficients
%   X to best fit the nonlinear functions in FUN to the data YDATA (in the 
%   least-squares sense). FUN accepts inputs X and XDATA and returns a
%   vector (or matrix) of function values F, where F is the same size as
%   YDATA, evaluated at X and XDATA. NOTE: FUN should return FUN(X,XDATA)
%   and not the sum-of-squares sum((FUN(X,XDATA)-YDATA).^2).
%   ((FUN(X,XDATA)-YDATA) is squared and summed implicitly in the
%   algorithm.) 
%
%   X = lsqcurvefit2(FUN,X0,XDATA,YDATA,LB,UB) defines a set of lower and
%   upper bounds on the design variables, X, so that the solution is in the
%   range LB <= X <= UB. Use empty matrices for LB and UB if no bounds
%   exist. Set LB(i) = -Inf if X(i) is unbounded below; set UB(i) = Inf if
%   X(i) is unbounded above.
%
%   X = lsqcurvefit2(FUN,X0,XDATA,YDATA,LB,UB,OPTIONS) minimizes with the
%   default parameters replaced by values in OPTIONS, an argument created
%   with the optimoptions2 function. See optimoptions2 for details. Use the
%   SpecifyObjectiveGradient option to specify that FUN also returns a
%   second output argument J that is the Jacobian matrix at the point X. If
%   FUN returns a vector F of m components when X has length n, then J is
%   an m-by-n matrix where J(i,j) is the partial derivative of F(i) with
%   respect to x(j). (Note that the Jacobian J is the transpose of the
%   gradient of F.)
%
%   X = lsqcurvefit2(PROBLEM) solves the non-linear least squares problem 
%   defined in PROBLEM. PROBLEM is a structure with the function FUN in 
%   PROBLEM.objective, the start point in PROBLEM.x0, the 'xdata' in 
%   PROBLEM.xdata, the 'ydata' in PROBLEM.ydata, the lower bounds in 
%   PROBLEM.lb, the upper bounds in PROBLEM.ub, the options structure in 
%   PROBLEM.options, and solver name 'lsqcurvefit2' in PROBLEM.solver. Use 
%   this syntax to solve at the command line a problem exported from 
%   OPTIMTOOL. 
%
%   [X,RESNORM] = lsqcurvefit2(FUN,X0,XDATA,YDATA,...) returns the value of 
%   the squared 2-norm of the residual at X: sum {(FUN(X,XDATA)-YDATA).^2}.
%
%   [X,RESNORM,RESIDUAL] = lsqcurvefit2(FUN,X0,...) returns the value of 
%   residual, FUN(X,XDATA)-YDATA, at the solution X. 
%
%   [X,RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit2(FUN,X0,XDATA,YDATA,...)
%   returns an EXITFLAG that describes the exit condition. Possible values
%   of EXITFLAG and the corresponding exit conditions are listed below. See
%   the documentation for a complete description.
%
%     1  lsqcurvefit2 converged to a solution.
%     2  Change in X too small.
%     3  Change in RESNORM too small.
%     4  Computed search direction too small.
%     0  Too many function evaluations or iterations.
%    -1  Stopped by output/plot function.
%    -2  Bounds are inconsistent.
%
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = lsqcurvefit2(FUN,X0,XDATA,
%   YDATA,...) returns a structure OUTPUT with the number of iterations 
%   taken in OUTPUT.iterations, the number of function evaluations in 
%   OUTPUT.funcCount, the algorithm used in OUTPUT.algorithm, the number of
%   CG iterations (if used) in OUTPUT.cgiterations, the first-order 
%   optimality (if used) in OUTPUT.firstorderopt, and the exit message in 
%   OUTPUT.message.
%
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = lsqcurvefit2(FUN,X0,XDATA,
%   YDATA,...) returns the set of Lagrangian multipliers, LAMBDA, at the 
%   solution: LAMBDA.lower for LB and LAMBDA.upper for UB.
%
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit2(FUN,
%   X0,XDATA,YDATA,...) returns the Jacobian of FUN at X.
%
%   Examples
%     FUN can be specified using @:
%        xdata = [5;4;6];          % example xdata
%        ydata = 3*sin([5;4;6])+6; % example ydata
%        x = lsqcurvefit2(@myfun, [2 7], xdata, ydata)
%
%   where myfun is a MATLAB function such as:
%
%       function F = myfun(x,xdata)
%       F = x(1)*sin(xdata)+x(2);
%
%   FUN can also be an anonymous function:
%       x = lsqcurvefit2(@(x,xdata) x(1)*sin(xdata)+x(2),[2 7],xdata,ydata)
%
%   If FUN is parameterized, you can use anonymous functions to capture the 
%   problem-dependent parameters. Suppose you want to solve the 
%   curve-fitting problem given in the function myfun, which is 
%   parameterized by its second argument c. Here myfun is a MATLAB file 
%   function such as
%
%       function F = myfun(x,xdata,c)
%       F = x(1)*exp(c*xdata)+x(2);
%
%   To solve the curve-fitting problem for a specific value of c, first 
%   assign the value to c. Then create a two-argument anonymous function 
%   that captures that value of c and calls myfun with three arguments. 
%   Finally, pass this anonymous function to lsqcurvefit2:
%
%       xdata = [3; 1; 4];           % example xdata
%       ydata = 6*exp(-1.5*xdata)+3; % example ydata     
%       c = -1.5;                    % define parameter 
%       x = lsqcurvefit2(@(x,xdata) myfun(x,xdata,c),[5;1],xdata,ydata)
%
%   See also optimoptions2, LSQNONLIN, FSOLVE, @, INLINE.

%   Copyright 1990-2015 The MathWorks, Inc.

% ------------Initialization----------------
defaultopt = struct(...
    'Algorithm','trust2-region-reflective',...
    'DerivativeCheck','off',...
    'Diagnostics','off',...
    'DiffMaxChange',Inf,...
    'DiffMinChange',0,...
    'Display','final',...
    'FinDiffRelStep', [], ...
    'FinDiffType','forward',...    
    'FunValCheck','off',...
    'InitDamping', 0.01, ...
    'Jacobian','off',...
    'JacobMult',[],... 
    'JacobPattern','sparse(ones(Jrows,Jcols))',...
    'MaxFunEvals',[],...
    'MaxIter',400,...
    'MaxPCGIter','max(1,floor(numberOfVariables/2))',...
    'OutputFcn',[],...
    'PlotFcns',[],...
    'PrecondBandWidth',Inf,...
    'ScaleProblem','none',...
    'TolFun', 1e-6,...
    'TolFunValue', 1e-6, ...    
    'TolPCG',0.1,...
    'TolX',1e-6,...
    'TypicalX','ones(numberOfVariables,1)',...
    'UseParallel',false );

% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(FUN,'defaults')
   xCurrent = defaultopt;
   return
end

if nargin < 7
    options = [];
    if nargin < 6
        UB = [];
        if nargin < 5
            LB = [];
        end
    end
end

problemInput = false;
if nargin == 1
    if isa(FUN,'struct')
        problemInput = true;
        [FUN,xCurrent,XDATA,YDATA,LB,UB,options] = separateOptimStruct(FUN);
    else % Single input and non-structure.
        error(message('optimlib:lsqcurvefit2:InputArg'));
    end
end

% Prepare the options for the solver
[options, optionFeedback] = prepareOptionsForSolver2(options, 'lsqcurvefit2');

if nargin < 4 && ~problemInput
    error(message('optimlib:lsqcurvefit2:NotEnoughInputs'));
end

% Check for non-double inputs
msg = isoptimargdbl2('lsqcurvefit2', {'X0','YDATA','LB','UB'}, ...
                                  xCurrent,YDATA,  LB,  UB);
if ~isempty(msg)
    error('optimlib:lsqcurvefit2:NonDoubleInput',msg);
end

caller = 'lsqcurvefit2';
[funfcn_x_xdata,mtxmpy_xdata,flags,sizes,funValCheck,xstart,lb,ub,EXITFLAG, ...
    Resnorm,FVAL,LAMBDA,JACOB,OUTPUT,earlyTermination] = ...
    lsqnsetup2(FUN,xCurrent,LB,UB,options,defaultopt,caller,nargout,length(varargin));
if earlyTermination
    return % premature return because of problem detected in lsqnsetup2()
end

xCurrent(:) = xstart; % reshape back to user shape before evaluation
funfcn = funfcn_x_xdata; % initialize user functions funfcn, which depend only on x
% Catch any error in user objective during initial evaluation only
switch funfcn_x_xdata{1}
    case 'fun'
        try
            initVals.F = feval(funfcn_x_xdata{3},xCurrent,XDATA,varargin{:});
        catch userFcn_ME
            optim_ME = MException('optimlib:lsqcurvefit2:InvalidFUN', ...
                getString(message('optimlib:lsqcurvefit2:InvalidFUN')));
            userFcn_ME = addCause(userFcn_ME,optim_ME);
            rethrow(userFcn_ME)
        end
        initVals.J = [];
        funfcn{3} = @objective; 
    case 'fungrad'
        try
            [initVals.F,initVals.J] = feval(funfcn_x_xdata{3},xCurrent,XDATA,varargin{:});
        catch userFcn_ME
            optim_ME = MException('optimlib:lsqcurvefit2:InvalidFUN', ...
                getString(message('optimlib:lsqcurvefit2:InvalidFUN')));
            userFcn_ME = addCause(userFcn_ME,optim_ME);
            rethrow(userFcn_ME)
        end
        funfcn{3} = @objectiveAndJacobian;
        funfcn{4} = @objectiveAndJacobian;
    case 'fun_then_grad'
        try
            initVals.F = feval(funfcn_x_xdata{3},xCurrent,XDATA,varargin{:});
        catch userFcn_ME
            optim_ME = MException('optimlib:lsqcurvefit2:InvalidFUN', ...
                getString(message('optimlib:lsqcurvefit2:InvalidFUN')));
            userFcn_ME = addCause(userFcn_ME,optim_ME);
            rethrow(userFcn_ME)
        end
        try    
            initVals.J = feval(funfcn_x_xdata{4},xCurrent,XDATA,varargin{:});
        catch userFcn_ME
            optim_ME = MException('optimlib:lsqcurvefit2:InvalidFUNJac', ...
                getString(message('optimlib:lsqcurvefit2:InvalidFUNJac')));
            userFcn_ME = addCause(userFcn_ME,optim_ME);
            rethrow(userFcn_ME)
        end
        funfcn{3} = @objective;
        funfcn{4} = @jacobian;
    otherwise
        error(message('optimlib:lsqcurvefit2:UndefCallType'))
end

% Check for non-double data typed values returned by user functions 
if ~isempty( isoptimargdbl2('lsqcurvefit2', {'F','J'}, initVals.F, initVals.J) )
    error('optimlib:lsqcurvefit2:NonDoubleFunVal',getString(message('optimlib:commonMsgs:NonDoubleFunVal','lsqcurvefit2')));
end

if ~isequal(size(initVals.F),size(YDATA))
    error(message('optimlib:lsqcurvefit2:YdataSizeMismatchFunVal'))
end
initVals.F = initVals.F - YDATA; % preserve initVals.F shape until after subtracting YDATA

% Pass functions that depend only on x: funfcn and jacobmult
[xCurrent,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = ...
    lsqncommon2(funfcn,xCurrent,lb,ub,options,defaultopt,caller,...
    initVals,sizes,flags,@jacobmult,optionFeedback,varargin{:}); 

%-------------------------------------------------------------------------------
% Nested functions that depend only on x and capture the constant values
% xdata and ydata, and also varargin
     function F = objective(x,varargin)
         F = feval(funfcn_x_xdata{3},x,XDATA,varargin{:});
         F = F - YDATA;
     end
     function [F,J] = objectiveAndJacobian(x,varargin)
         % Function value and Jacobian returned together
         [F,J] = feval(funfcn_x_xdata{3},x,XDATA,varargin{:});
         F = F - YDATA;
     end
     function J = jacobian(x,varargin)
         % Jacobian returned in separate function
         J = feval(funfcn_x_xdata{4},x,XDATA,varargin{:});
     end
     function W = jacobmult(Jinfo,Y,flag,varargin)
         W = feval(mtxmpy_xdata,Jinfo,Y,flag,XDATA,varargin{:});
     end
 end

