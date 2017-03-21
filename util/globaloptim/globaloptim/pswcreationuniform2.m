function swarm = pswcreationuniform2(problemStruct)
%PSWCREATIONUNIFORM2 Creates the initial positions for PARTICLESWARM2 algorithm.
%   SWARM = PSWCREATIONUNIFORM2(PROBLEMSTRUCT) 
%   Creates the initial positions that PARTICLESWARM2 will use.
%
%   Example:
%     options = optimoptions2('particleswarm2', 'CreationFcn',@pswcreationuniform2);

%   Copyright 2012-2015 The MathWorks, Inc.

nvars = problemStruct.nvars;
options = problemStruct.options;
% Determine finite bounds for the initial particles based on the problem's
% bounds and options.InitialSwarmSpan.
[lb,ub] = determinePositionInitBounds(problemStruct.lb, problemStruct.ub, ...
    options.InitialSwarmSpan);

numParticles = options.SwarmSize;
numInitPositions = size(options.InitialSwarm, 1);
numPositionsToCreate = numParticles - numInitPositions;

% Initialize particles to be created
swarm = zeros(numParticles,nvars);

% Use initial particles provided already
if numInitPositions > 0
    swarm(1:numInitPositions,:) = options.InitialSwarm;
end

% Create remaining particles, randomly sampling within lb and ub
span = ub - lb;
swarm(numInitPositions+1:end,:) = repmat(lb,numPositionsToCreate,1) + ...
    repmat(span,numPositionsToCreate,1) .* rand(numPositionsToCreate,nvars);

% Error if any values are not finite
if ~all(isfinite(swarm(:)))
    error(message('globaloptim:pswcreationuniform2:positionNotFinite'));
end
end

function [lb,ub] = determinePositionInitBounds(lb,ub,initialSwarmSpan)
% Update lb and ub using positionInitSpan, so that initial bounds are
% always finite
lbFinite = isfinite(lb);
ubFinite = isfinite(ub);
lbInf = ~lbFinite;
ubInf = ~ubFinite;

% If lb and ub are both finite, do not update the bounds.

% If lb & ub are both infinite, center the range around 0.
idx = lbInf & ubInf;
lb(idx) = -initialSwarmSpan(idx)/2;
ub(idx) = initialSwarmSpan(idx)/2;

% If only lb is finite, start the range at lb.
idx = lbFinite & ubInf;
ub(idx) = lb(idx) + initialSwarmSpan(idx);

% If only ub is finite, end the range at ub.
idx = lbInf & ubFinite;
lb(idx) = ub(idx) - initialSwarmSpan(idx);
end
