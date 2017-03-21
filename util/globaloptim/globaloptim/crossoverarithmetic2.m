function xoverKids  = crossoverarithmetic2(parents,options,GenomeLength,FitnessFcn,unused,thisPopulation)
%CROSSOVERARITHMETIC2 Crossover operator for constrained problems.

%   Copyright 2005-2007 The MathWorks, Inc.


% How many children to produce?
nKids = length(parents)/2;
% Allocate space for the kids
xoverKids = zeros(nKids,GenomeLength);
% To move through the parents twice as fast as the kids are
% being produced, a separate index for the parents is needed
index = 1;
% for each kid...
for i=1:nKids
    % get parents
    r1 = parents(index);
    index = index + 1;
    r2 = parents(index);
    index = index + 1;
    % Children are arithmetic mean of two parents
    alpha = rand;
    xoverKids(i,:) = alpha*thisPopulation(r1,:) + (1-alpha)*thisPopulation(r2,:);
end
