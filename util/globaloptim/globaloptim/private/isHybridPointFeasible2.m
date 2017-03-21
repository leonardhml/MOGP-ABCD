function [isFeas,maxconstraint] = isHybridPointFeasible2(hybridConInfo, hybridSolver, hybridOptions)
%isHybridPointFeasible2 Determine whether the point returned from the hybrid
%                      solver is feasible.
%
%   [isFeas,maxconstraint] = isHybridPointFeasible2(hybridConInfo,
%   hybridSolver) takes information regarding the constraints,
%   hybridConInfo, at the point returned from the hybrid solver,
%   hybridSolver. isHybridPointFeasible2 determines whether the hybrid
%   solver solution is feasible using hybridSolver's default constraint
%   tolerance. The largest constraint violation is returned in
%   maxconstraint.
%
%   The contents2 of hybridConInfo depend on hybridSolver:
%
%   hybridSolver            +  hybridConInfo
%   ------------------------|----------------------------------------
%   fmincon2, patternsearch2  |  Constraint violation returned from
%                           |  hybridSolver.
%   fgoalattain             |  {X, Aineq, bineq, Aeq, beq, lb, ub}
%                           |  where X is the solution returned by Hybrid
%                           |  solver and the remaining elements are the 
%                           |  linear constraints in the problem.
%
%   isFeas = isHybridPointFeasible2(hybridConViol, hybridSolver,
%   hybridOptions) uses the constraint tolerance specified in the
%   hybridOptions structure.

%   Copyright 2010 The MathWorks, Inc.

% If the hybrid options are not specified, then set them to be empty.
if nargin < 3
    hybridOptions = [];
end

% Get the constraint tolerance used by the hybrid function.
hybridTol = i_getHybridFcnTolCon(hybridSolver, hybridOptions);

% Test whether the point returned by the hybrid function is feasible w.r.t
% the hybrid function's constraint tolerance.
if strcmpi(hybridSolver, 'fgoalattain')
    % hybridConInfo contains the final point returned from fgoalattain and
    % the linear constraints in this case. We need to evaluate the linear
    % constraints rather than using the constraint violation as reported by
    % the solver. This is because a point returned from fgoalattain may be
    % feasible w.r.t. to the problem constraints, but infeasible w.r.t. the
    % constraint added by fgoalattain. Such a point may be an improvement
    % on that returned by gamultiobj2 and should be returned by the hybrid
    % function.
    %
    % Note that the only caller of fgoalattain as a hybrid solver is
    % gamultiobj2, which in turn only supports linear constraints.
    hybridConInfo{end+1} = hybridTol;
    [isFeas,maxconstraint] = isTrialFeasible2(hybridConInfo{:});
else
    % hybridConInfo contains either the maximum constraint violation
    % reported by the hybrid solver or it is empty in this case.
    isFeas = isempty(hybridConInfo) || hybridConInfo <= hybridTol;
    maxconstraint = hybridConInfo;
end

function TolCon = i_getHybridFcnTolCon(hybridFcn, hybridOpts)

switch hybridFcn
    case 'patternsearch2'
        
        % Ensure that the supplied options structure contains TolCon.
        hybridOpts = psoptimset2(hybridOpts);
        
        % Get the default options.
        defaultHybridOpts = patternsearch2('defaults');
        
        % Get TolCon
        TolCon = psoptimget2(hybridOpts, 'TolCon', defaultHybridOpts, 'fast');
        
    case {'fmincon2', 'fgoalattain'}
        
        % Ensure that the supplied options structure contains TolCon.
        hybridOpts = optimset(hybridOpts);
        
        % Get the default options.
        switch hybridFcn
            case 'fmincon2'
                defaultHybridOpts = fmincon2('defaults');
            case 'fgoalattain'
                defaultHybridOpts = fgoalattain('defaults');
        end
        
        % Get TolCon
        TolCon = optimget(hybridOpts, 'TolCon', defaultHybridOpts, 'fast');
                
end
