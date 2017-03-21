function problemStructure = createOptimProblem2(solvername,varargin)
%createOptimProblem2 Create an optimization problem.
%   createOptimProblem2 creates a structure that contains the solver
%   name, solver options and problem description.
%
%   PROBLEM = createOptimProblem2(SOLVERNAME, 'PROP1', VAL1, ... 
%   'PROP2', VAL2, ...) specifies a set of property-value pairs relevant
%   for the Optimization Toolbox solver, SOLVERNAME, to construct the
%   optimization problem structure, PROBLEM. SOLVERNAME is a required
%   argument and can be either one of 'fmincon2', 'fminunc', 'lsqnonlin' or
%   'lsqcurvefit2'. The property names, can be one or more of the problem
%   structure field names: 'objective', 'x0', 'xdata', 'ydata', 'Aineq',
%   'bineq', 'Aeq', 'beq', 'lb', 'ub', 'nonlcon', 'options'. The field
%   names specified  
%        (1) need to be relevant to the chosen solver, SOLVERNAME,
%        (2) can be entered in any order,
%        (3) can consist of only non-empty input arguments for the chosen
%        solver.
%
%   The full input argument lists for the supported solvers and
%   createOptimProblem2 calls to create the corresponding problem structures
%   are:
%
%       fmincon2(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
%
%       createOptimProblem2('fmincon2','objective', FUN, 'x0', X0, ...
%       'Aineq', A, 'bineq', b, 'Aeq', Aeq, 'beq', beq, 'lb', LB, ...
%       'ub', UB, 'nonlcon', NONLCON, 'options', OPTIONS)
%
%       FMINUNC(FUN,X0,OPTIONS)
%
%       createOptimProblem2('fminunc','objective', FUN, 'x0', X0, ...
%       'options', OPTIONS)
%
%       LSQNONLIN(FUN,X0,LB,UB,OPTIONS)
%
%       createOptimProblem2('lsqnonlin','objective', FUN, 'x0', X0, ...
%       'lb', LB, 'ub', UB, 'options', OPTIONS)
%
%       lsqcurvefit2(FUN,X0,XDATA,YDATA,LB,UB,OPTIONS)
%
%       createOptimProblem2('lsqcurvefit2','objective', FUN, 'x0', X0, ...
%       'xdata', XDATA, 'ydata', YDATA, 'lb', LB, 'ub', UB, ...
%       'options', OPTIONS)
%   
%   Examples:
%
%   Create an optimization problem structure for solver fmincon2:
%       problem = createOptimProblem2('fmincon2', 'objective', @myfun, ...
%        'x0', [1;1], 'lb', [-3;-3], 'ub', [3;3]);
%
%       Here myfun is a function such as
%  
%           function f = myfun(x)      
%           f = sin(x(1))*x(2)^3;  
%
%   Create an optimization problem structure for solver lsqnonlin:
%       problem = createOptimProblem2('lsqnonlin', ...
%        'objective', @myfun, 'x0', [2 3 4]);
%
%       Here myfun is a function such as
%  
%           function f = myfun(x)      
%           f = sin(x);  
%
%   See also MULTISTART2, GLOBALSEARCH2, fmincon2, FMINUNC, LSQNONLIN, lsqcurvefit2

%   Copyright 2009-2014 The MathWorks, Inc.

useValues.options = [];
problemProp = varargin;

if nargin < 1
    error(message('globaloptim:createOptimProblem2:AtLeastOneInput'));
else
    if isa(solvername,'char')
    % First argument should be the name of these valid solvers
        validSolvers = {'fmincon2','fminunc','lsqnonlin','lsqcurvefit2'};
        if ~any(strcmp(solvername,validSolvers))
            error(message('globaloptim:createOptimProblem2:InvalidSolver'));
        end
        useValues.solver = solvername;
    else
            error(message('globaloptim:createOptimProblem2:SolverNameNotString'));
    end
end
if mod(length(problemProp),2) ~= 0
    error(message('globaloptim:createOptimProblem2:NotAPair'));
end
while length(problemProp) >=2
    prop = problemProp{1};
    if ~ischar(prop)
        error(message('globaloptim:createOptimProblem2:PropertyNameNotString'));
    end    
    val = problemProp{2};
    problemProp = problemProp(3:end);
    propName = i_isExactFieldName(prop);
    if ~isempty(propName)
        useValues.(propName) = val;
    else
        % No match at all
        error(message('globaloptim:createOptimProblem2:InvalidPropertyName', prop));
    end
end

% It is easy for users to incorrectly pass a string for options instead
% of a structure. In this case we will just use the default values but
% issue a warning telling the user that they have passed a string for
% options. Note that we are issuing a warning rather than an error
% since that's what we do when users pass strings for options via the
% solver's options input.
if ~isempty(useValues.options) && ischar(useValues.options)
    warning(message('globaloptim:createOptimProblem2:OptionsPassedAsString'));
end

problemStructure = createProblemStruct(useValues.solver,[],useValues);
% Get the solver options and update the changed ones thru
% useValues.options
createdOptions = false;
if ~isempty(useValues.options)
    if isa(useValues.options, 'optim.options.SolverOptions2')
        if strcmp(useValues.options.SolverName, useValues.solver)
            problemStructure.options = useValues.options;
        else
            warning(message('globaloptim:createOptimProblem2:WrongSolverOptions', ...
                upper(useValues.options.SolverName), ...
                upper(useValues.solver), ...
                upper(useValues.options.SolverName)));
            problemStructure.options = optimoptions2(useValues.solver, ...
                useValues.options);
        end
        createdOptions = true;
    elseif isstruct(useValues.options)
        % If options are specified as a structure, we assume that an option
        % set to empty means use the default value. This is the optimset
        % standard. 
        %
        % Also, almost all options do not allow empty as a valid value. For
        % those options that do allow empty as a valid value, empty is the
        % default value. As such, the option will not be set in the try
        % statement below or will be set to the default value.
        % 
        % We will remove options with empty values from the list of options
        % to be set. This means that any options that are on the
        % deprecation path will not be set to their default value, which
        % throws a warning.        
        allOptions = fieldnames(useValues.options);
        allValues = struct2cell(useValues.options);
        idxNotEmpty = ~cellfun('isempty', allValues);
        allOptions = allOptions(idxNotEmpty);
        problemStructure.options = optimoptions2(useValues.solver);
        for i = 1:length(allOptions)
            try %#ok
                problemStructure.options.(allOptions{i}) = ...
                    useValues.options.(allOptions{i});
            end
        end
        createdOptions = true;
    else
        createdOptions = false;
    end
end

if ~createdOptions
    problemStructure.options = optimoptions2(useValues.solver);
end

function correctPName = i_isExactFieldName(prop)
% Check whether an identifier matches exactly to a field name.

correctIdents = {'objective','x0','Aineq','bineq','Aeq','beq','lb','ub', ...
    'nonlcon','options','xdata','ydata'};
correctPName = '';
if any(strcmp(prop,correctIdents))
    correctPName = prop;
end
