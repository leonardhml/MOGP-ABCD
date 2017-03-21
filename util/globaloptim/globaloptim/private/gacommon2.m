function [x,fval,exitFlag,output,population,scores,FitnessFcn,nvars,Aineq,bineq,Aeq,beq,lb,ub, ...
    NonconFcn,options,Iterate,type] = gacommon2(nvars,fun,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,intcon,options,user_options,output)
%gacommon2 Common validation tasks for GA2 and GAMULTIOBJ2
%
%   This function is private to GA2 and GAMULTIOBJ2

%   Copyright 2007-2015 The MathWorks, Inc.

% Determine the population type
stringSet2('PopulationType',options.PopulationType,{'doubleVector','custom','bitString'});

% Determine the 'type' of the problem
if ~isempty(nonlcon)
    type = 'nonlinearconstr';
    % Determine the sub-problem type for the constrained problem (used in ALPS)
    if (~isempty(Aeq)   && ~all(Aeq(:) == 0))   || ~isempty(beq) || ... 
       (~isempty(Aineq) && ~all(Aineq(:) == 0)) || ~isempty(bineq)
        subtype = 'linearconstraints';
    elseif (~isempty(lb) && any(isfinite(lb))) || ...
           (~isempty(ub) && any(isfinite(ub)))
        subtype = 'boundconstraints';
    else
        subtype = 'unconstrained';
    end   
    % If Aeq or Aineq is not empty, then problem has linear constraints.
elseif (~isempty(Aeq)   && ~all(Aeq(:) == 0))   || ~isempty(beq) || ... 
       (~isempty(Aineq) && ~all(Aineq(:) == 0)) || ~isempty(bineq)
       type = 'linearconstraints';
    % This condition satisfies bound constraints
elseif (~isempty(lb) && any(isfinite(lb))) || ...
       (~isempty(ub) && any(isfinite(ub)))
    type = 'boundconstraints';
    % If all constraints fields are empty then it is unconstrained
else
    type = 'unconstrained';
end
% Make sure 'custom' or 'bitString' don't bring in unsupported constraints
% with them.
if strcmpi(options.PopulationType,'bitString') && ...
        ~strcmpi(type,'unconstrained')
    warning(message('globaloptim:gacommon2:bitStringConstraintsSupport'));
    type = 'unconstrained';
    Aineq = []; bineq = []; Aeq = []; beq = []; lb = []; ub = [];
elseif strcmpi(options.PopulationType,'custom') && ...
        ~any(strcmpi(type,{'unconstrained','boundconstraints'}))
    warning(message('globaloptim:gacommon2:customConstraintsSupport'));
    if (~isempty(lb) || ~isempty(ub))
        type = 'boundconstraints';
        Aineq = []; bineq = []; Aeq = []; beq = [];
    else
        type = 'unconstrained';
        Aineq = []; bineq = []; Aeq = []; beq = []; lb = []; ub = [];
    end
end
output.problemtype = type;
% If nonlinear constraints, then subtype is needed to process linear
% constraints (see function preProcessLinearConstr2)
if strcmp(type,'nonlinearconstr')
    type = subtype;
end

% Remember the random number state used
dflt = RandStream.getGlobalStream;
output.rngstate = struct('state',{dflt.State}, 'type',{dflt.Type});
% Initialize other fields of output
output.generations = 0;
output.funccount   = 0;
output.message   = '';
output.maxconstraint = [];

% Validate2 options and fitness function
[options,nvars,FitnessFcn,NonconFcn] = validate2(options,type,nvars,fun,nonlcon,user_options);

% Perform check on initial population, score, and range
options = initialPopulationCheck2(options);

if ~strcmp(output.problemtype,'unconstrained') 
    % Determine a start point
    if ~isempty(options.InitialPopulation)
        x = options.InitialPopulation(1,:);
    else
        x = randn(1,nvars);
    end
    Iterate.x = x(:);
else
    Iterate.x = [];
end
% Initialize output
fval = [];
x = [];
population = [];
scores = [];

% Bound correction
[lb,ub,msg,exitFlag] = checkbound2(lb,ub,nvars);
if exitFlag < 0
    output.message = msg;
    if options.Verbosity > 0
        fprintf('%s\n',msg)
    end
    return;
end
% Linear constraints correction
[Iterate.x,Aineq,bineq,Aeq,beq,lb,ub,msg,exitFlag] = preProcessLinearConstr2( ...
    Iterate.x,Aineq,bineq,Aeq,beq,lb,ub,options.TolCon,nvars,type,options.Verbosity);
if exitFlag < 0
    output.message = msg;
    if options.Verbosity > 0
        fprintf('%s\n',msg)
    end
    return;
end

% Verify that individuals in InitialPopulation are feasible
if ~isempty(options.InitialPopulation) && ~strcmp(type,'unconstrained')
    feasible = isTrialFeasible2(options.InitialPopulation,Aineq,bineq,Aeq,beq,lb,ub,options.TolCon);
    options.InitialPopulation(~feasible,:) = [];
    try % InitialScores may not be present
        options.InitialScores(~feasible,:) = [];
    catch
    end
end
% If initial population is empty at this point then we stuff the feasible point
% found before. Don't do this if integer constraints are present - we have
% found that this degrades the performance of some of the integer
% constrained problems.
if isempty(options.InitialPopulation) && ~isempty(Iterate.x) && isempty(intcon)
    options.InitialPopulation(1,:) = Iterate.x';
end

% Validate2 nonlinear constraints
[LinearConstr, Iterate,nineqcstr,neqcstr,ncstr] = constrValidate2(NonconFcn, ...
    Iterate,Aineq,bineq,Aeq,beq,lb,ub,intcon,type,options);
options.LinearConstr = LinearConstr;

% Add integer constraint field to the options structure.
options.IntegerVars = [];

% Make sure that bounds and PopInitRange are consistent 
options = checkPopulationInitRange2(lb,ub,options);

% Print some diagnostic information if asked for
if options.Verbosity > 2
    % Find numObj (for gamultiobj2) by making one call to fitness function
    % (if Initial Score is not present)
    % Determine who is the caller
    callStack = dbstack;
    [~,caller] = fileparts(callStack(2).file);
    if strcmp(caller,'gamultiobj2')
        if isempty(options.InitialScores)
            if isempty(options.InitialPopulation)
                % Do we have a feasible point
                options.InitialPopulation(1,:) =  randn(1,nvars);
            end
            % Evaluate InitialPopulation(1,:)
             score = FitnessFcn(options.InitialPopulation(1,:));
             options.InitialScores = score(:)';
        end
        numObj = size(options.InitialScores,2);
    else
        numObj = 1;
    end
    gadiagnose2(fun,nonlcon,nvars,nineqcstr,neqcstr,ncstr,numObj,user_options,caller);
end

