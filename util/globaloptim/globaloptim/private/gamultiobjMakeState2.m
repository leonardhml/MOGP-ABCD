function state = gamultiobjMakeState2(GenomeLength,FitnessFcn,ConstrFcn,~,options)
%gamultiobjMakeState2 Create an initial population and fitness scores

%   Copyright 2007-2015 The MathWorks, Inc.

% A variety of data used in various places
state.Generation = 0;		% current generation counter
state.StartTime = tic;	    % current clock time
state.StopFlag = []; 		% reason for termination
state.FunEval = 0;
state.Selection = [];       % selection indices
state.mIneq = 0;
state.mEq = 0;
state.mAll = 0;
popSize = sum(options.PopulationSize);

% If InitialPopulation is partly empty we will use the creation function to
% generate population (CreationFcn can utilize InitialPopulation)
if sum(options.PopulationSize) ~= size(options.InitialPopulation,1)
    state.Population = feval(options.CreationFcn,GenomeLength,FitnessFcn,options,options.CreationFcnArgs{:});
else % Initial population was passed in
    state.Population = options.InitialPopulation;
end

% Evaluate fitness function to get the number of objectives
try
    Score = FitnessFcn(state.Population(1,:));
catch userFcn_ME
    gads_ME = MException('globaloptim:gamultiobjMakeState2:fitnessCheck', ...
       getString(message('globaloptim:gamultiobjMakeState2:fitnessCheck')));
    userFcn_ME = addCause(userFcn_ME,gads_ME);
    rethrow(userFcn_ME)
end

if ~isempty(ConstrFcn)
% Evaluate nonlinear constraint function to get constraint counts and for
% error-checking purposes
    try
        [cineq,ceq] = ConstrFcn(state.Population(1,:));
    catch userFcn_ME
        gads_ME = MException('globaloptim:gamultiobjMakeState2:nonlconCheck', ...
           getString(message('globaloptim:gamultiobjMakeState2:nonlconCheck')));
        userFcn_ME = addCause(userFcn_ME,gads_ME);
        rethrow(userFcn_ME)
    end
    
    state.mIneq = numel(cineq);
    state.mEq = numel(ceq);
    state.mAll = state.mIneq + state.mEq;
end

state.FunEval = state.FunEval + 1;
Score = Score(:)';
numObj = numel(Score);
% Size of InitialScore and Score should match
if ~isempty(options.InitialScores) && numObj ~= size(options.InitialScores,2)
    error(message('globaloptim:gamultiobjMakeState2:initScoreSize','size(InitialScore,2)'));
end

% Calculate score for state.Population
totalPopulation = size(state.Population,1);
initScoreProvided = size(options.InitialScores,1);
state.C      = zeros(totalPopulation,state.mIneq);
state.Ceq    = zeros(totalPopulation,state.mEq);
state.isFeas = true(totalPopulation,1);
state.maxLinInfeas = zeros(totalPopulation,1);

if totalPopulation ~= initScoreProvided
    individualsToEvaluate = totalPopulation - initScoreProvided;
    state.Score  =  zeros(totalPopulation,numObj);

    if initScoreProvided > 0
        state.Score(1:initScoreProvided,:) = options.InitialScores(1:initScoreProvided,:);
    end
    
    if strcmpi(options.Vectorized, 'off')
        % Score remaining members of the population
        [Score,C,Ceq,isFeas] = objAndConVectorizer2(state.Population(initScoreProvided+1:end,:), ...
            FitnessFcn,ConstrFcn,numObj,state.mIneq,state.mEq,options.SerialUserFcn, ...
            true(individualsToEvaluate,1),options.TolCon);
    else       
        if state.mAll > 0
            [C,Ceq] = ConstrFcn(state.Population(initScoreProvided+1:end,:));
            % Make sure sizes of empties are correct
            if state.mEq == 0
                Ceq = reshape(Ceq,size(C,1),0);
            elseif state.mIneq == 0
                C = reshape(C,size(Ceq,1),0);
            end
            isFeas = isNonlinearFeasible2(C,Ceq,options.TolCon);
        else
            C = zeros(totalPopulation-initScoreProvided,0);
            Ceq = zeros(totalPopulation-initScoreProvided,0);
            isFeas = true(totalPopulation-initScoreProvided,1);
        end
        Score = FitnessFcn(state.Population(initScoreProvided+1:end,:));
        if size(Score,1) ~= individualsToEvaluate
           error(message('globaloptim:gamultiobjMakeState2:fitnessVectorizedCheck', ... 
                'Vectorized','on'));
        end         
    end
    
    state.Score(initScoreProvided+1:end,:)  = Score;
    state.C(initScoreProvided+1:end,:)      = C;
    state.Ceq(initScoreProvided+1:end,:)    = Ceq;
    state.isFeas(initScoreProvided+1:end,:) = isFeas;
    state.FunEval = state.FunEval+ individualsToEvaluate;          % number of function evaluations
else
    state.Score = options.InitialScores;    
end

% If the initial score (objective values) were given for any members of the
% population, then we must go back and evaluate the constraints on those
% points since it can't be given by the user.
if (initScoreProvided > 0) && (state.mAll > 0)
    fakeFitness = @(x) NaN(numObj,size(x,2));
    if strcmpi(options.Vectorized, 'off')
        [~,C,Ceq,isFeas] = objAndConVectorizer2(state.Population(1:initScoreProvided,:), ...
            fakeFitness,ConstrFcn,numObj,state.mIneq,state.mEq,options.SerialUserFcn, ...
            true(initScoreProvided,1),options.TolCon);     
    else
        [C,Ceq] = ConstrFcn(state.Population(1:initScoreProvided,:));
        % Make sure sizes of empties are correct
        if state.mEq == 0
            Ceq = reshape(Ceq,size(C,1),0);
        elseif state.mIneq == 0
            C = reshape(C,size(Ceq,1),0);
        end
        isFeas = isNonlinearFeasible2(C,Ceq,options.TolCon);
    end
    state.C(1:initScoreProvided,:) = C;
    state.Ceq(1:initScoreProvided,:) = Ceq;
    state.isFeas(1:initScoreProvided,:) = isFeas;
    state.FunEval = state.FunEval+ initScoreProvided;
end

% Partial population is allowed for 'doubleVector' and 'bitString'
% population type so make population of appropriate size
if ~strcmpi(options.PopulationType,'custom')
    lens = size(state.Population,1);
    npop = sum(options.PopulationSize);
    if npop > lens
        population = zeros(npop,GenomeLength);
        population(1:lens,:) = state.Population;
        population(lens+1:end,:) = repmat(state.Population(end,:),(npop-lens),1);
        scores = zeros(npop,numObj);
        C = zeros(npop,state.mIneq);
        Ceq = zeros(npop,state.mEq);
        scores(1:lens,:) = state.Score;
        C(1:lens,:) = state.C;
        Ceq(1:lens,:) = state.Ceq;
        scores(lens+1:end,:) = repmat(state.Score(end,:),(npop-lens),1);
        C(lens+1:end,:) = repmat(state.C(end,:),(npop-lens),1);
        Ceq(lens+1:end,:) = repmat(state.Ceq(end,:),(npop-lens),1);
        state.Population = population;
        state.Score = scores;
        state.C = C;
        state.Ceq = Ceq;
        state.isFeas = [state.isFeas; repmat(state.isFeas(end,:),(npop-lens),1)];
        state.maxLinInfeas = [state.maxLinInfeas; zeros(npop-lens,1)];
    else
        state.Population(npop+1:end,:) = [];
        state.Score(npop+1:end,:) = [];
        state.C(npop+1:end,:) = [];
        state.Ceq(npop+1:end,:) = [];
        state.isFeas(npop+1:end,:) = [];
        state.maxLinInfeas(npop+1:end,:) = [];
    end
end

% Update list of feasible population members accounting for linear constraints.
% NOTE: maxLinInfeas is initialized to zeros, so no need for "else" case.
if isfield(options,'LinearConstr') && options.LinearConstr.linconCheck
    [isLinFeas,state.maxLinInfeas] = isTrialFeasible2(state.Population,options.LinearConstr.Aineq, ...
        options.LinearConstr.bineq,options.LinearConstr.Aeq,options.LinearConstr.beq, ...
        options.LinearConstr.lb,options.LinearConstr.ub,options.TolCon);
    state.isFeas = state.isFeas & isLinFeas;
end

% Before computing rank, etc., we must set the score of infeasible
% population members to Inf, so as to prevent mutation/crossover on those
% members.
if state.mAll > 0 || (isfield(options,'LinearConstr') && options.LinearConstr.linconCheck)
    state.Score(~state.isFeas,:) = Inf(sum(~state.isFeas),size(state.Score,2));
end

% Get the rank and Distance measure of the population
[state.Population,state.Score,state.Rank,state.Distance,state.C,state.Ceq, ...
 state.isFeas,state.maxLinInfeas] = rankAndDistance2(state.Population,state.Score, ...
                state.C,state.Ceq,state.isFeas,state.maxLinInfeas,options,popSize);

% One step GA2 (selection, crossover and then mutation)
% How many crossover offspring will there be from each source?
nXoverKids = round(options.CrossoverFraction * popSize);
nMutateKids = popSize - nXoverKids;
% how many parents will we need to complete the population?
nParents = 2 * nXoverKids + nMutateKids;
% Selection.
parents = feval(options.SelectionFcn,[state.Rank,state.Distance],nParents,options,options.SelectionFcnArgs{:});
% If selection function does not return the correct parents, error here
if (length(parents) ~= nParents) || min(parents) <= 0 || max(parents) > popSize
   error(message('globaloptim:gamultiobjMakeState2:selectionFcnError', nParents, popSize));
end
% Shuffle to prevent locality effects. 
parents = parents(randperm(length(parents)));

xoverKids  = feval(options.CrossoverFcn,parents(1:(2 * nXoverKids)),options,GenomeLength, ...
     FitnessFcn,state.Score,state.Population,options.CrossoverFcnArgs{:});
mutateKids = feval(options.MutationFcn,parents((1 + 2 * nXoverKids):end),options,GenomeLength, ...
    FitnessFcn,state,state.Score,state.Population,options.MutationFcnArgs{:});

nextPopulation = [xoverKids ; mutateKids ];

if isfield(options,'LinearConstr') && options.LinearConstr.linconCheck
    [nextIsFeas,nextMaxLinInfeas] = isTrialFeasible2(nextPopulation,options.LinearConstr.Aineq, ...
        options.LinearConstr.bineq,options.LinearConstr.Aeq,options.LinearConstr.beq, ...
        options.LinearConstr.lb,options.LinearConstr.ub,options.TolCon);
else
    nextIsFeas = true(size(nextPopulation,1),1);
    nextMaxLinInfeas = zeros(size(nextPopulation,1),1);
end

if strcmpi(options.Vectorized, 'off')
    % Score remaining members of the population
    [nextScore,nextC,nextCeq,nextIsFeas] = objAndConVectorizer2(nextPopulation, ...
        FitnessFcn,ConstrFcn,numObj,state.mIneq,state.mEq,options.SerialUserFcn, ...
        nextIsFeas,options.TolCon);
else
    if state.mAll > 0
        [nextC,nextCeq] = ConstrFcn(nextPopulation);
        % Make sure sizes of empties are correct
        if state.mEq == 0
            nextCeq = reshape(nextCeq,size(nextC,1),0);
        elseif state.mIneq == 0
            nextC = reshape(nextC,size(nextCeq,1),0);
        end
        nextIsFeas = nextIsFeas & isNonlinearFeasible2(nextC,nextCeq,options.TolCon);
    else
        nextC = zeros(size(nextPopulation,1),0);
        nextCeq = zeros(size(nextPopulation,1),0);        
    end
    % Only evaluate feasible points
    nextScore = Inf(size(nextPopulation,1),numObj);
    nextScore(nextIsFeas,:) = FitnessFcn(nextPopulation(nextIsFeas,:));
end

state.FunEval =  state.FunEval + size(nextPopulation,1);
 
% Repeat for each sub-population (element of the PopulationSize vector)
offset = 0;
totalPop = options.PopulationSize;
for pop = 1:length(totalPop)
    populationSize =  totalPop(pop);
    subpopIndex = 1 + (offset:(offset + populationSize - 1));

    % Combine old and new population
    population = [state.Population(subpopIndex,:);nextPopulation(subpopIndex,:)];
    score_old = state.Score(subpopIndex,:);
    score = [state.Score(subpopIndex,:);nextScore(subpopIndex,:)];
    C = [state.C(subpopIndex,:); nextC(subpopIndex,:)];
    Ceq = [state.Ceq(subpopIndex,:); nextCeq(subpopIndex,:)];
    isFeas = [state.isFeas(subpopIndex,:); nextIsFeas(subpopIndex,:)];
    maxLinInfeas = [state.maxLinInfeas(subpopIndex,:); nextMaxLinInfeas(subpopIndex,:)];
    popSize = numel(subpopIndex);

    % Calculate the rank and distance measure of each individual 
    [state.Population(subpopIndex,:),state.Score(subpopIndex,:), ...
        state.Rank(subpopIndex,:),state.Distance(subpopIndex,:), ...
        state.C(subpopIndex,:),state.Ceq(subpopIndex,:), ...
        state.isFeas(subpopIndex,:),state.maxLinInfeas(subpopIndex,:)] = ...
        rankAndDistance2(population,score,C,Ceq,isFeas,maxLinInfeas,options,popSize);
    offset = offset + populationSize;

    % Calculate average distance and spread
    [state.AverageDistance(pop), state.Spread(state.Generation+1,pop)] = ...
        distanceAndSpread2(state.Distance(subpopIndex,:),state.Rank(subpopIndex,:), ...
                          state.Score(subpopIndex,:),score_old);
end

