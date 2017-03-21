function Population = gacreationnonlinearfeasible2(GenomeLength,~,options,ConstrFcn,varargin)
%GACREATIONNONLINEARFEASIBLE2 Creates the initial population for GA2.
%   POP = GACREATIONNONLINEARFEASIBLE2(NVARS,FITNESSFCN,OPTIONS,CONSTRFCN) 
%   creates an initial population that GA2 will then evolve into a solution.
%
%   Example:
%     options = optimoptions2('ga2','CreationFcn',@gacreationnonlinearfeasible2);
%
%    (Note: If calling gamultiobj2, replace 'ga2' with 'gamultiobj2')

%   Copyright 2015 The MathWorks, Inc.

totalPopulation = sum(options.PopulationSize);
userInitPop = options.InitialPopulation;
options.InitialPopulation = [];
initPopProvided = size(userInitPop,1);

% Parse optional inputs
thisFcn = 'GACREATIONNONLINEARFEASIBLE2';
parser = inputParser();
parser.FunctionName = thisFcn;
parser.addParameter('NumStartPts',0, ...
    @(x)validateattributes(x,{'numeric'},{'scalar','nonempty'},thisFcn,'NumStartPts') );
parser.addParameter('UseParallel',~options.SerialUserFcn, ...
    @(x)validateattributes(x,{'logical'},{'scalar','nonempty'},thisFcn,'UseParallel') );
parser.addParameter('SolverOpts',optimoptions2('fmincon2','Algorithm','sqp','TolCon',options.TolCon), ...
    @(x)validateattributes(x,{'optim.options.fmincon2'}, ...
    {'scalar','nonempty'},thisFcn,'SolverOpts') );
parser.parse(varargin{:});

individualsToCreate = totalPopulation - initPopProvided;
if parser.Results.NumStartPts ~= 0
    individualsToCreate = min(individualsToCreate,parser.Results.NumStartPts);
end

if individualsToCreate <= 0
    return;
end

A = []; b = []; Aeq = []; beq = []; lb = []; ub = []; 
if isfield(options,'LinearConstr')
    % Extract linear constraint data and append columns for the slack
    % variable of the minimum infeasibility formulation.
    b = options.LinearConstr.bineq;
    beq = options.LinearConstr.beq;
    A = options.LinearConstr.Aineq;
    Aeq = options.LinearConstr.Aeq;
    lb = options.LinearConstr.lb(:);
    ub = options.LinearConstr.ub(:);
end

% Create the rest of the initial points from gacreationuniform2
% NOTE: overwrite (or create) the "type" field as 'boundconstraints' to
% ensure points will be generated.
options.LinearConstr.type = 'boundconstraints';
options.PopulationSize = individualsToCreate;
initPts = gacreationuniform2(GenomeLength,[],options);

% Find out which points are already feasible

% Eval nonlinear constraints on initPts
if (isfield(options,'UserVectorized') && options.UserVectorized) || ...
   (~isfield(options,'UserVectorized') && strcmpi(options.Vectorized,'on'))
    [C,Ceq] = ConstrFcn(initPts);
else
    C = []; Ceq = [];
    for k = 1:individualsToCreate
        [thisC,thisCeq] = ConstrFcn(initPts(k,:));
        C   = [C;   thisC(:)']; 
        Ceq = [Ceq; thisCeq(:)'];
    end
end
isFeas = isTrialFeasible2(initPts,A,b,Aeq,beq,lb,ub,options.TolCon) & ...
         isNonlinearFeasible2(C,Ceq,options.TolCon);

% Add feasible points into Population
Population = Inf(totalPopulation,GenomeLength);
numFeas = sum(isFeas);
% Remove feasible points from initPts, and stuff the feasible initPts in
% the middle of Population so we can run through a contiguous array
% (helpful for the parfor).
Population(individualsToCreate-numFeas+1:individualsToCreate,:) = initPts(isFeas,:);
initPts(isFeas,:) = [];

% Only setup feasibility problem if needed
if numFeas < individualsToCreate   
    numPtsToRun = individualsToCreate-numFeas;
    
    % Reformulate the problem to find the minimizer of constraint violation
    maxInfeasObj = @(x) x(end);
    maxInfeasGrad = @(x) [zeros(GenomeLength,1); 1];
    A = [A  zeros(numel(b),1)];
    Aeq = [Aeq  zeros(numel(beq),1)];
    lb = [lb; 0];
    ub = [ub; Inf];
        
    % Turn display off, and GradObj on
    opts = parser.Results.SolverOpts;
    opts.Display = 'off';
    opts.GradObj = 'on';
    % Use OutputFcn to monitor time limit. Make sure we don't surpass TimeLimit.
    if isfinite(options.TimeLimit)
        opts.OutputFcn = @(a,b,c) toc(state.StartTime) >= options.TimeLimit;
    end
    
    warnstate = warning;
    warning off;
    % Solve a set of feasibility problems using fmincon2
    if ~parser.Results.UseParallel
        for k = 1:numPtsToRun
            Population(k,:) = callSolver({maxInfeasObj,maxInfeasGrad},initPts(k,:), ...
                A,b,Aeq,beq,lb,ub,ConstrFcn,opts);
        end
    else
        parfor k = 1:numPtsToRun
            Population(k,:) = callSolver({maxInfeasObj,maxInfeasGrad},initPts(k,:), ...
                A,b,Aeq,beq,lb,ub,ConstrFcn,opts);
        end
    end
    warning(warnstate);
end
% Fill in rest of population with user-provided Population, if any
if initPopProvided > 0
    Population(individualsToCreate+1:individualsToCreate+initPopProvided,:) = userInitPop;
end


%--------------------------------------------------------------------------
function [c,ceq] = minFeasConstr(x,origConFcn)
% This constraint wraps the original problem's constraints with an added
% slack variable; transforming the problem into a "minimum infeasibility"
% problem:
% min_{x,gamma}  gamma
% subject to:   A*x <= b
%               Aeq*x = beq
%               lb <= x <= ub
%               ceq(x) - gamma = 0
%               c(x) - gamma <= 0
%               gamma >= 0
%
gamma = x(end);
[c,ceq] = origConFcn(x(1:end-1));
c = c - gamma;
ceq = ceq - gamma;

    
%--------------------------------------------------------------------------
function x = callSolver(obj,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
% NOTE: It's expected that this function always returns row vectors.
try
    [x,~,flag] = fmincon2(obj,[x0 0],A,b,Aeq,beq,lb,ub,@(y)minFeasConstr(y,nonlcon),options);
    x = x(1:numel(x0)); % Discard slack
catch
    % fmincon2 errored, we can only return the initial point
    x = x0;
    return;
end
if flag < 0
    % fmincon2 didn't achieve feasibility. We should return the result if
    % possible. In all likelihood, fmincon2 should at least maintain the
    % same level of infeasibility as x0, if not improve.
    if isempty(x)
        x = x0;
    end
end

