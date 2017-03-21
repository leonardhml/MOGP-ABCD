function xoverKids  = crossoversinglepoint2(parents,options,GenomeLength,~,~,thisPopulation)
%CROSSOVERSINGLEPOINT2 Single point crossover.
%   XOVERKIDS = CROSSOVERSINGLEPPOINT(PARENTS,OPTIONS,GENOMELENGTH, ...
%   FITNESSFCN,SCORES,THISPOPULATION) creates the crossover children XOVERKIDS
%   of the given population THISPOPULATION using the available parents
%   PARENTS. A single crossover point for each child is chosen at random.
%   The child has the genes of the first parent up to this point and the genes
%   of the other parent after this point.
%
%   Example:
%    Create an options structure using CROSSOVERSINGLEPOINT2 as the crossover
%    function
%     options = optimoptions2('ga2', 'CrossoverFcn', @crossoversinglepoint2);
%    (Note: If calling gamultiobj2, replace 'ga2' with 'gamultiobj2') 

%   Copyright 2003-2015 The MathWorks, Inc.


% How many children to produce?
nKids = length(parents)/2;
% Extract information about linear constraints, if any
linCon = options.LinearConstr;
constr = ~isequal(linCon.type,'unconstrained');
% Allocate space for the kids
xoverKids = zeros(nKids,GenomeLength);

% To move through the parents twice as fast as thekids are
% being produced, a separate index for the parents is needed
index = 1;

for i=1:nKids
    % get parents
    parent1 = thisPopulation(parents(index),:);
    index = index + 1;
    parent2 = thisPopulation(parents(index),:);
    index = index + 1;

    % cut point is AFTER this index.
    xOverPoint = ceil(rand * (length(parent1) - 1));
    % make one child
    xoverKids(i,:) = [ parent1(1:xOverPoint),parent2(( xOverPoint + 1 ):  end )  ];
    % Make sure that offspring are feasible w.r.t. linear constraints
    if constr
        feasible = isTrialFeasible2(xoverKids(i,:)',linCon.Aineq,linCon.bineq,linCon.Aeq, ...
                            linCon.beq,linCon.lb,linCon.ub,options.TolCon);
        if ~feasible % Kid is not feasible
            % Children are arithmetic mean of two parents (feasible w.r.t
            % linear constraints)
            alpha = rand;
            xoverKids(i,:) = alpha*parent1 + (1-alpha)*parent2;
        end
    end
end