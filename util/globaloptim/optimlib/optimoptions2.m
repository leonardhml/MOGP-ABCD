function options = optimoptions2(solver, varargin)
%optimoptions2  Create/modify optimization options
%
%   OPTIONS = optimoptions2(SOLVER) creates optimization options, OPTIONS,
%   with the option parameters set to the default values relevant to the
%   optimization solver named in SOLVER, for example 'fmincon2'. For a list
%   of valid optimization solvers, enter
%     doc optimoptions2
%   to open the documentation page for optimoptions2.
%
%   OPTIONS = optimoptions2(SOLVER,'PARAM1',VALUE1,...) creates default
%   optimization options for SOLVER with the named parameters altered with
%   the specified values. OPTIONS can also be created using optimoptions2
%   and setting the parameters via dot notation. For example:
%
%   % Create and set parameters using optimoptions2
%   options = optimoptions2('fminunc','StepTolerance',0.01,'Display','iter')
%
%   % Create same options using optimoptions2 and dot notation
%   options = optimoptions2('fminunc'); 
%   options.StepTolerance = 0.01;
%   options.Display = 'iter'
%
%   OPTIONS = optimoptions2(OLDOPTS,'PARAM1',VALUE1,...) creates a copy of
%   OLDOPTS with the named parameters altered with the specified values.
%
%   OPTIONS = optimoptions2(SOLVER, OTHEROPTS) first creates default
%   optimization options for SOLVER. The values of any parameters which are
%   common to OTHEROPTS and OPTIONS are copied to OPTIONS.
%
%   Note to see the parameters for a specific function, check the
%   documentation page for that function. For instance, enter
%     doc fmincon2
%   to open the reference page for fmincon2.
%
%      Example - setting options for a optimization solver.
%       % Create optimization options with the default parameters for
%       % fmincon2
%       options = optimoptions2('fmincon2');
%
%       % Set OptimalityTolerance to 1e-3
%       options = optimoptions2(options, 'OptimalityTolerance', 1e-3); 
%
%       % Set the Display option to 'iter' and StepTolerance to 1e-4
%       options.Display = 'iter';
%       options.StepTolerance = 1e-4
%
%       % Run the fmincon2 solver with these options:
%       [x, fval, exitflag] = fmincon2(@(x)peaks(x(1), x(2)), ...
%           [1 2], [1 1], 1, [], [], [-5 -5], [5 5], [], options)
%
%     Example - creating from another set of options
%
%      % Create some options for lsqnonlin
%      lsqOpts = optimoptions2('lsqnonlin', 'MaxFunctionEvaluations', 1000)
%
%      % Use the above options for a call to lsqnonlin
%      x = lsqnonlin(@(x) sin(3*x),[1 4], [], [], lsqOpts)
%
%      % Now we want to add a nonlinear constraint to the above
%      % least squares problem. We will need to call fmincon2, but we
%      % would like to use the same optimization options.
%
%      % First, create options for fmincon2 from the lsqnonlin options.
%      fmcOpts = optimoptions2('fmincon2', lsqOpts)
%     
%      % Next, set the algorithm to 'interior-point'.
%      fmcOpts.Algorithm = 'interior-point';
%
%      % Call fmincon2
%      [x, fval] = fmincon2(@(x) (sin(3*x(1)))^2 + (sin(3*x(2)))^2, [1 4], ...
%             [], [], [], [], [], [], @(x)deal((x(1) - 1)^2 - x(2)), fmcOpts)

%   Copyright 2012-2015 The MathWorks, Inc.

% Error if no input is given
if nargin == 0
    error(message('optimlib:optimoptions2:NotEnoughInputs'));
end

% Determine which solver we're creating options for. Note that invalid
% values of solver are caught in createSolverOptions2.
if isa(solver, 'optim.options.SolverOptions2')
    solverName = solver.SolverName;
    % If the user specified an old solver options object as the first
    % input, we need to pass this to the solver options object constructor.
    varargin = [{solver}, varargin]; 
elseif isa(solver, 'function_handle')
    solverName = func2str(solver);
else
    solverName = solver;
end

% Create the solver options object
try
    options = optim.options.createSolverOptions2(solverName, varargin{:});
catch ME
    switch ME.identifier
        case 'MATLAB:InputParser:ParamMustBeChar'
            [linkToSolverCshStart, linkToSolverCshEnd] = ...
                i_createLinkToSolverCsh(solverName);
            error(message('optimlib:optimoptions2:ParameterMustBeChar', ...
                upper(solverName), linkToSolverCshStart, linkToSolverCshEnd));
        case 'MATLAB:InputParser:UnmatchedParameter'
            [linkToSolverCshStart, linkToSolverCshEnd] = ...
                i_createLinkToSolverCsh(solverName);
            % Get the parameter that hasn't been matched
            idxQuote = regexp(ME.message, '''');
            error(message('optimlib:optimoptions2:UnmatchedParameter', ...
                ME.message(idxQuote(1):idxQuote(2)), ...
                upper(solverName), upper(solverName), linkToSolverCshStart, linkToSolverCshEnd));
        case 'MATLAB:InputParser:ParamMissingValue'
            % Get the parameter that hasn't been matched
            idxQuote = regexp(ME.message, '''');
            error(message('optimlib:optimoptions2:ParameterMissingValue', ...
                ME.message(idxQuote(1):idxQuote(2))));            
        otherwise
            throw(ME);
    end
end

function [docFileTagStart, docFileTagEnd] = i_createLinkToSolverCsh(solverName)
%i_createLinkToSolverCsh Create hyperlink tags to a solver doc page
%
%   [DOCFILETAGSTART, DOCFILETAGEND] = I_CREATELINKTOSOLVERCSH(SOLVERNAME)
%   creates hyperlink tags to the options section of the function reference
%   documentation page for the named solver. For example, for fmincon2
%
%   docFileTagStart = <a href = "matlab: helpview('<DOCROOT>/toolbox/optim/helptargets.map','fmincon_opts')">
%   docFileTagEnd = </a>
%
%   If links are not enabled in the MATLAB session, then the start and end
%   tags are empty strings.

% Determine toobox
if isempty(strfind(getfield(functions(str2func(solverName)),'file'),'globaloptim'))
    toolboxName = 'optim';
else
    toolboxName = 'gads';
end

[~,docFileTagStart, docFileTagEnd] = addLink2('',toolboxName,'helptargets.map', ...
                                             sprintf('%s_opts',solverName),false);


