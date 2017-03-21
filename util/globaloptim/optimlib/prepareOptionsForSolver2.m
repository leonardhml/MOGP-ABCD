function [options, optionFeedback] = prepareOptionsForSolver2(options, solverName)
%

%prepareOptionsForSolver2 Prepare options for solver
%
%   [OPTIONS, OPTIONSFEEDBACK] = prepareOptionsForSolver2(OPTIONS,
%   SOLVERNAME) performs tasks to ensure the options are set up for use by
%   the solver. The following tasks are performed:-
%
%   * If a user has passed a SolverOptions2 object, first ensure that it is
%     a set of options for fseminf. 
%     
%   * If required by the solver, prepare strings to give feedback to users
%     on options they have or have not set. These are used in the exit
%     messages.
%
%   * If a user has passed a SolverOptions2 object, now extract the options
%     structure from the object.

%   Copyright 2012-2015 The MathWorks, Inc.

% If a user has passed a structure, we cannot tell whether a user wrote the
% code that passes the structure before or after 16a. In this case, we
% assume that the code was written before 16a and set TolFunValue to
% TolFun, if TolFunValue is not in options.
if isstruct(options) && isfield(options, 'TolFun') && ...
        ~isfield(options, 'TolFunValue')
    options.TolFunValue = options.TolFun;
end 

% If a user has passed a SolverOptions2 object, first ensure that it is a
% set of options for the solver 
if isa(options, 'optim.options.SolverOptions2')
    options = convertForSolver(options, solverName);
end

% If required, prepare strings to give feedback to users on options they
% have or have not set. These are used in the exit messages.
if nargout > 1
    optionFeedback = createOptionFeedback2(options);
end

% If a user has passed a SolverOptions2 object, now extract the options
% structure from the object
if isa(options, 'optim.options.SolverOptions2')
    options = extractOptionsStructure(options);
end
