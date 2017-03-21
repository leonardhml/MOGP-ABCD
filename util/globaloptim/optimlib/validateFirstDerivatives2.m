function validateFirstDerivatives2(funfcn,confcn,x, ...
        lb,ub,options,finDiffFlags,sizes,varargin)
%

% validateFirstDerivatives2 Helper function that validates first derivatives of
% objective, nonlinear inequality, and nonlinear equality gradients against
% finite differences. The finite-difference calculation is done according to
% options.FinDiffType.
%
% This function assumes that objective and constraint functions, options,
% and flags for finitedifferences2 have been validated before calling.
% funfcn and confcn must be cell-array outputs of optimfcnchk2. lb and ub
% must be vectors of length number of variables. options is a structure
% that must contain non-empty fields: GradObj, GradConstr, TypicalX,
% FinDiffRelStep, DiffMinChange, DiffMaxChange, and ScaleProblem. 
% sizes is a structure that must contain the fields: 
% - nVar: the number of variables
% - nFun: the number of functions in a system of equations (or
% least-squares problem)
% - mNonlinEq: the number of nonlinear equality constraints
% - mNonlinIneq: the number of nonlinear inequality constraints
% - xShape: the original shape of the variable xCurrent as expected by
% the user-supplied functions funfcn and confcn.
% tol: Component-wise relative difference in gradients checked 
%            against this tolerance

%   Copyright 2007-2015 The MathWorks, Inc.

nVar = sizes.nVar;
if isempty(confcn)
    confcn = {''}; % Protection for the unconstrained solvers
end

% Component-wise relative difference in gradients checked against this tolerance
tol = 1e-6;

% --- Select point at which to evaluate the user-supplied derivatives ---

% We will first perturb by a random point in an effort to minimize the
% effects of a poorly chosen initial point (x0) and the likelihood of
% cancellations in the user-supplied derivatives.
randDir = 1 - 2*rand(nVar,1);  % A random direction (random numbers between +/- 1)
xDerivChk = x(:) + sqrt(tol)*max(abs(x(:)),1).*randDir;

% Second, we will attempt to center the perturbed x0 within the finite
% bounds. Variables without finite bounds will remain at the perturbed x0.
xDerivChk = shiftInitPtToInterior(nVar,xDerivChk,lb,ub,Inf);

% Now, evaluate the user objective and constraint functions at the chosen
% point.
[fval,grad,cIneq,cEq,JacCineqTrans,JacCeqTrans,sizes] = evalUserFuns(funfcn,confcn,...
                                                    xDerivChk,sizes,varargin{:});

mNonlinIneq = sizes.mNonlinIneq;  
mNonlinEq = sizes.mNonlinEq;

% Print marker to visually separate derivative check information from other
% optimization diagnostics.
separatorTxt = repmat('_',60,1);

% Create other common untranslated strings
derivCheck = 'CheckGradients';
indent = '  ';

% Pre-construct error message in case finite-difference estimate and
% user-supplied derivative do not match within the given tolerance

enableLinks = feature('hotlinks') && ~isdeployed;
if enableLinks
linkCmd = '<a href = "matlab: helpview([docroot ''/toolbox/optim/msg_csh/optim_msg_csh.map''],''%s'',''CSHelpWindow'');">';
cshTagStart = sprintf(linkCmd,'checkgradients_failed');
cshTagEnd = '</a>';
else
    cshTagStart = '';
    cshTagEnd = '';
end
 
derivDiffError = MException(message('optimlib:validateFirstDerivatives2:InvalidGrad', ...
    cshTagStart,derivCheck,cshTagEnd,options.FinDiffType,sprintf('%g',tol)));
                                         
fprintf('%s\n',separatorTxt);
fprintf(['   ',getString(message('optimlib:validateFirstDerivatives2:DerivCheckHeader',derivCheck))]);
        
if strcmpi(options.GradObj,'on') 
    % input to finitedifferences2()
    if finDiffFlags.isGrad
        grad = grad(:);
        grad_fd = zeros(nVar,1); 
    else
        grad_fd = zeros(sizes.nFun,nVar);
    end
    grad_fd = finitedifferences2(xDerivChk,funfcn{3},[],lb,ub,fval,[],[],1:nVar, ...
        options,sizes,grad_fd,[],[],finDiffFlags,[],varargin{:});

    % Vector of objective gradient relative error
    relGradError = full(abs(grad_fd - grad)./max(1.0,abs(grad)));    
    [maxDiff,i,j] = findRowColIndicesOfMaxElement(relGradError);
    fprintf(getString(message('optimlib:validateFirstDerivatives2:ObjFunHeader')));
    fprintf(getString(message('optimlib:validateFirstDerivatives2:RelMaxDiff', ...
        sprintf('%g',maxDiff) )));
    % Note, we can use any(any(...)) here as relGradError is guaranteed to
    % be a vector or a 2-d matrix
    if any(any(relGradError > tol)) 
        fprintf([indent,getString(message('optimlib:validateFirstDerivatives2:UsrObjGradVal', ...
            i,j,sprintf('%g',full(grad(i,j))) ))]);
        fprintf([indent,getString(message('optimlib:validateFirstDerivatives2:FinDiffObjGradVal', ...
            i,j,sprintf('%g',full(grad_fd(i,j))) ))]);
        % Display message that derivative check failed, with links to context-sensitive help.
        disp(getString(message('optimlib:validateFirstDerivatives2:DerivCheckFailed',cshTagStart,derivCheck,cshTagEnd)));
        fprintf('%s\n\n',separatorTxt);
        throw(derivDiffError);
    end
end

% If there are nonlinear constraints and their derivatives are provided,
% validate them
if strcmpi(options.GradConstr,'on') && sizes.mNonlinIneq + sizes.mNonlinEq > 0
    JacCineqTrans_fd = zeros(nVar,mNonlinIneq); % input to finitedifferences2()
    JacCeqTrans_fd = zeros(nVar,mNonlinEq); % input to finitedifferences2()
    
    [~,JacCineqTrans_fd,JacCeqTrans_fd] = finitedifferences2(xDerivChk,[],confcn{3}, ...
        lb,ub,[],cIneq(:),cEq(:),1:nVar,options,sizes,[],JacCineqTrans_fd, ...
        JacCeqTrans_fd,finDiffFlags,[],varargin{:});
    if sizes.mNonlinIneq > 0
        % Matrix of nonlinear inequality constraint gradient relative error
        % JacCineqTrans_fd is full so JaccIneqError will be full - store it as a full matrix        
        relJacCineqError = full(abs(JacCineqTrans - JacCineqTrans_fd))./max(1.0,abs(JacCineqTrans));
        [maxDiff,i,j] = findRowColIndicesOfMaxElement(relJacCineqError);
        fprintf(getString(message('optimlib:validateFirstDerivatives2:NonlinIneqHeader')));
        fprintf(getString(message('optimlib:validateFirstDerivatives2:RelMaxDiff', ...
            sprintf('%g',maxDiff) )));
        if any(any( relJacCineqError > tol ))
            fprintf([indent,getString(message('optimlib:validateFirstDerivatives2:UsrNonlinIneqGradVal', ...
                i,j,sprintf('%g',full(JacCineqTrans(i,j))) ))]);
            fprintf([indent,getString(message('optimlib:validateFirstDerivatives2:FinDiffNonlinIneqGradVal', ...
                i,j,sprintf('%g',full(JacCineqTrans_fd(i,j))) ))]);
            % Display message that derivative check failed, with links to context-sensitive help.
            disp(getString(message('optimlib:validateFirstDerivatives2:DerivCheckFailed',cshTagStart,derivCheck,cshTagEnd)));
            fprintf('%s\n\n',separatorTxt);
            throw(derivDiffError);
        end
    end
    
    if sizes.mNonlinEq > 0
        % Matrix of nonlinear equality constraint gradient relative error
        % JacCeqTrans_fd is full so JacCeqError will be full - store it as a full matrix
        relJacCeqError = full(abs(JacCeqTrans - JacCeqTrans_fd))./max(1.0,abs(JacCeqTrans)); 
        [maxDiff,i,j] = findRowColIndicesOfMaxElement(relJacCeqError);
        fprintf(getString(message('optimlib:validateFirstDerivatives2:NonlinEqHeader')))
        fprintf(getString(message('optimlib:validateFirstDerivatives2:RelMaxDiff', ...
            sprintf('%g',maxDiff) )));
        if any(any( relJacCeqError > tol ))
            fprintf([indent,getString(message('optimlib:validateFirstDerivatives2:UsrNonlinIneqGradVal', ...
                i,j,sprintf('%g',full(JacCeqTrans(i,j))) ))]);
            fprintf([indent,getString(message('optimlib:validateFirstDerivatives2:FinDiffNonlinIneqGradVal', ...
                i,j,sprintf('%g',full(JacCeqTrans_fd(i,j))) ))]);
            % Display message that derivative check failed, with links to context-sensitive help.
            disp(getString(message('optimlib:validateFirstDerivatives2:DerivCheckFailed',cshTagStart,derivCheck,cshTagEnd)));
            fprintf('%s\n\n',separatorTxt);
            throw(derivDiffError);
        end
    end
end

% Print marker to visually separate derivative check information from other
% optimization diagnostics.
fprintf(getString(message('optimlib:validateFirstDerivatives2:DerivCheckPassed',derivCheck)));
fprintf('%s\n\n',separatorTxt);

%-------------------------------------------------------------------------
function [maxVal,i,j] = findRowColIndicesOfMaxElement(A)
% Helper function that finds indices (i,j) of the maximum element of matrix A.
% It also returns the maximum element, maxVal

% Find max element by columns
[col_max_val,row_idx] = max(A,[],1);
% Find max element by rows
[maxVal,col_idx] = max(col_max_val);
i = row_idx(col_idx);
j = col_idx;

%-------------------------------------------------------------------------
function [fval,grad,cIneq,cEq,JacCineqTrans,JacCeqTrans,sizes] = evalUserFuns( ...
                                            funfcn,confcn,X,sizes,varargin)
%

% Reshape X into the size required by the user-functions
X = reshape(X,sizes.xShape);
                                    
% Create exceptions for possible failures in evaluating the objective or
% constraint functions
derivCheck = 'CheckGradients';
obj_ME = MException('optimlib:validateFirstDerivatives2:ObjectiveError', ...
    getString(message('optimlib:validateFirstDerivatives2:ObjectiveError',derivCheck)));
gradObj_ME = MException('optimlib:validateFirstDerivatives2:GradientError', ...
    getString(message('optimlib:validateFirstDerivatives2:GradientError',derivCheck)));
constr_ME = MException('optimlib:validateFirstDerivatives2:NonlconError', ...
    getString(message('optimlib:validateFirstDerivatives2:NonlconError',derivCheck)));
gradConstr_ME = MException('optimlib:validateFirstDerivatives2:NonlconFunOrGradError', ...
    getString(message('optimlib:validateFirstDerivatives2:NonlconFunOrGradError',derivCheck)));

% Default values in case there are no gradients
grad = zeros(0,1);
switch funfcn{1}
    case 'fun'
        try
            fval = feval(funfcn{3},X,varargin{:});
        catch userFcn_ME
            userFcn_ME = addCause(userFcn_ME,obj_ME);
            rethrow(userFcn_ME)
        end
    case {'fungrad','fungradhess'}
        try
            [fval,grad] = feval(funfcn{3},X,varargin{:});
        catch userFcn_ME
            userFcn_ME = addCause(userFcn_ME,obj_ME);
            rethrow(userFcn_ME)
        end
    case {'fun_then_grad','fun_then_grad_then_hess'}
        try
            fval = feval(funfcn{3},X,varargin{:});
        catch userFcn_ME
            userFcn_ME = addCause(userFcn_ME,obj_ME);
            rethrow(userFcn_ME)
        end
        try
            grad = feval(funfcn{4},X,varargin{:});
        catch userFcn_ME
            userFcn_ME = addCause(userFcn_ME,gradObj_ME);
            rethrow(userFcn_ME)
        end
end

JacCineqTrans = zeros(0,1);
JacCeqTrans = zeros(0,1);
% Evaluate constraints
switch confcn{1}
    case 'fun'
        try
            [cIneq,cEq] = feval(confcn{3},X,varargin{:});
        catch userFcn_ME
            userFcn_ME = addCause(userFcn_ME,constr_ME);
            rethrow(userFcn_ME)
        end
    case 'fungrad'
        try
            [cIneq,cEq,JacCineqTrans,JacCeqTrans] = feval(confcn{3},X,varargin{:});
        catch userFcn_ME
            userFcn_ME = addCause(userFcn_ME,constr_ME);
            rethrow(userFcn_ME)
        end
    case 'fun_then_grad'
        try
            [cIneq,cEq] = feval(confcn{3},X,varargin{:});
        catch userFcn_ME
            userFcn_ME = addCause(userFcn_ME,constr_ME);
            rethrow(userFcn_ME)
        end
        try
            [JacCineqTrans,JacCeqTrans] = feval(confcn{4},X,varargin{:});
        catch userFcn_ME
            userFcn_ME = addCause(userFcn_ME,gradConstr_ME);
            rethrow(userFcn_ME)
        end
    case ''
        % No nonlinear constraints. Reshaping of empty quantities is done later
        % in this file, where both cases, (i) no nonlinear constraints and (ii)
        % nonlinear constraints that have one type missing (equalities or
        % inequalities), are handled in one place
        cIneq = zeros(0,1);
        cEq = zeros(0,1);
end
cIneq = cIneq(:);
cEq = cEq(:);
sizes.mNonlinIneq = numel(cIneq);
sizes.mNonlinEq = numel(cEq);