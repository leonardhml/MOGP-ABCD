function Population = gacreationuniform2(GenomeLength,FitnessFcn,options)
%GACREATIONUNIFORM2 Creates the initial population for genetic algorithm.
%   POP = GACREATIONUNIFORM2(NVARS,FITNESSFCN,OPTIONS) Creates the
%   initial population that GA2 will then evolve into a solution.
%
%   Example:
%     options = optimoptions2('ga2','CreationFcn',@gacreationuniform2);
%
%    (Note: If calling gamultiobj2, replace 'ga2' with 'gamultiobj2')

%   Copyright 2003-2015 The MathWorks, Inc.
%

if strcmpi(options.PopulationType,'custom')
    error(message('globaloptim:gacreationuniform2:unknownPopulationType', options.PopulationType));
end

totalPopulation = sum(options.PopulationSize);
initPopProvided = size(options.InitialPopulation,1);
individualsToCreate = totalPopulation - initPopProvided;

if strcmpi(options.PopulationType,'doubleVector')
    % Initialize Population to be created
    Population = zeros(totalPopulation,GenomeLength);
    % Use initial population provided already
    if initPopProvided > 0
        Population(1:initPopProvided,:) = options.InitialPopulation;
    end
    % Create remaining population
    % problemtype is either 'unconstrained', 'boundconstraints', or
    % 'linearconstraints'. Nonlinear constrained algorithm 'ALGA' does not
    % create or use initial population of its own. It calls sub-problem
    % solvers (galincon2/gaunc2)
    if isfield(options,'LinearConstr')
        problemtype = options.LinearConstr.type;
    else
        problemtype = 'unconstrained';
    end
    % This function knows how to create initial population for
    % unconstrained and boundconstrained cases but does not know in
    % linearconstrained case.
    if ~strcmp(problemtype,'linearconstraints')
        range = options.PopInitRange;
        lowerBound = range(1,:);
        span = range(2,:) - lowerBound;
        Population(initPopProvided+1:end,:) = repmat(lowerBound,individualsToCreate,1) + ...
            repmat(span,individualsToCreate,1) .* rand(individualsToCreate,GenomeLength);
    else
        Population(initPopProvided+1:end,:) = []; 
    end
elseif strcmpi(options.PopulationType,'bitString')
    % Initialize Population to be created
    Population = ones(totalPopulation,GenomeLength);
    % Use initial population provided already
    if initPopProvided > 0
        Population(1:initPopProvided,:) = options.InitialPopulation;
    end
    % Create remaining population
    Population(initPopProvided+1:end,:) = double(0.5 > rand(individualsToCreate,GenomeLength));
end

if any(isnan(Population(:)))
    error(message('globaloptim:gacreationuniform2:populationIsNaN'));
elseif any(isinf(Population(:)))
    error(message('globaloptim:gacreationuniform2:populationIsInf'));
end

