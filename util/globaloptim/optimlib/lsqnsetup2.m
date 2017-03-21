function [funfcn,mtxmpy,flags,sizes,funValCheck,xstart,lb,ub,EXITFLAG,Resnorm,FVAL, ...
    LAMBDA,JACOB,OUTPUT,earlyTermination] = lsqnsetup2(FUN,xC,LB,UB,options, ...
    defaultopt,caller,callerNargout,lengthVarargin)
%

%lsqnsetup2 Helper function for lsqcurvefit2 and lsqnonlin.
% This function performs setup tasks common to both solvers.

%   Copyright 2008-2011 The MathWorks, Inc.

display = optimget(options,'Display',defaultopt,'fast');
flags.detailedExitMsg = ~isempty(strfind(display,'detailed'));
switch display
case {'off','none'}
    flags.verbosity = 0;
case {'iter','iter-detailed'}
    flags.verbosity = 2;
case {'final','final-detailed'}
    flags.verbosity = 1;
case 'testing'
    flags.verbosity = Inf;
otherwise
    flags.verbosity = 1;
end

funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');
flags.grad =  strcmp(optimget(options,'Jacobian',defaultopt,'fast'),'on');

if callerNargout > 5
   flags.computeLambda = true;
else 
   flags.computeLambda = false;
end

% Process objective function
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
    funfcn = lsqfcnchk2(FUN,caller,lengthVarargin,funValCheck,flags.grad);
else
    if strcmp(caller,'lsqnonlin')
        mexcptn = MException('optim:lsqnonlin:InvalidFunc', ...
            getString(message('optim:lsqnonlin:InvalidFunc')));
    else
        mexcptn = MException('optimlib:lsqcurvefit2:InvalidFunc', ...
            getString(message('optimlib:lsqcurvefit2:InvalidFunc')));
    end
    throwAsCaller(mexcptn);
end

% Process function JacobMult
mtxmpy = optimget(options,'JacobMult',defaultopt,'fast'); 
% Check if name clash
functionNameClashCheck2('JacobMult',mtxmpy,'atamult2','optimlib:lsqnsetup2:JacobMultNameClash');

% Use internal Jacobian-multiply function if user does not provide JacobMult function 
% or options.Jacobian is off
if isempty(mtxmpy) || (~strcmpi(funfcn{1},'fungrad') && ~strcmpi(funfcn{1},'fun_then_grad'))
    mtxmpy = @atamult2;
end

xstart = xC(:); % operate on column-vector xstart
sizes.nVar = numel(xstart);
[xstart,lb,ub,msg] = checkbounds2(xstart,LB,UB,sizes.nVar);
if ~isempty(msg)
    EXITFLAG = -2;
    [Resnorm,FVAL,LAMBDA,JACOB] = deal([]);
    OUTPUT.firstorderopt = [];
    OUTPUT.iterations = 0;
    OUTPUT.funcCount = 0;
    OUTPUT.cgiterations = [];
    OUTPUT.algorithm = ''; % Not known at this stage
    OUTPUT.message = msg;
    if flags.verbosity > 0
        disp(msg)
    end
    earlyTermination = true; % exit this helper function and return
else
    % All outputs need to be assigned to a value
    Resnorm = []; FVAL = []; LAMBDA = []; JACOB = []; OUTPUT = [];
    EXITFLAG = []; earlyTermination = false;
    
    % If components of initial x not within bounds, set those components
    % of initial point to a "box-centered" point
    xinitOutOfBounds_idx = xstart < lb | xstart > ub;
    if any(xinitOutOfBounds_idx)
        xstart = startx2(ub,lb,xstart,xinitOutOfBounds_idx);
    end
end


