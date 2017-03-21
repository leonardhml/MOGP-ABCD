classdef (Sealed) MultiStart2 < AbstractGlobalSolver2
%MultiStart2 A multi-start global optimization solver. 
%   A MultiStart2 object is a solver that attempts to locate multiple local
%   solutions (possibly the global one among them) to a given problem by
%   starting from different start points. Create this object by typing its
%   constructor, MultiStart2 and supplying parameter values (if different
%   from default). 
%
%   MS = MultiStart2 constructs a new multi-start optimization solver
%   with its properties set to defaults.
%
%   MS = MultiStart2(PROP, VAL, ...) specifies a set of property-value
%   pairs that are applied to the multi-start optimization solver before
%   creating it.
%
%   MS = MultiStart2(OLDMS, PROP, VAL, ...) creates a copy of the
%   MultiStart2 solver OLDMS with the named properties altered with the
%   specified values.
%
%   MS = MultiStart2(GS) constructs a new MultiStart2 solver and copies
%   the common parameter values in the GlobalSearch2 solver GS into the new
%   solver MS.
%
%   MultiStart2 properties:
%       Display             - detail of display
%       FunctionTolerance   - minimum distance between two separate 
%                             objective function values
%       XTolerance          - minimum distance between two separate points 
%       MaxTime             - time allowed to run the solver
%       UseParallel         - perform calculation in parallel mode
%       StartPointsToRun    - points that the solver will start from
%       OutputFcn           - user-defined output functions
%       PlotFcn             - functions to plot solver progress
%
%   MultiStart2 method:
%       run                 - run a local solver from multiple start points 
%                             for a given optimization problem
%
%
%   Typical workflow to run the MultiStart2 solver:
%   ==============================================
%   1. Set up the PROBLEM structure
%       PROBLEM = createOptimProblem2('fmincon2','objective',...)
%   2. Construct the MultiStart2 solver
%       MS = MultiStart2
%   3. Run the solver from NUMPOINTS random start points
%       run(MS,PROBLEM,NUMPOINTS)
%
%   Example:
%      Run fmincon2 from 20 random start points for the optimization problem 
%         minimize x(1)^2 + 4*sin(5*x(2)); subject to 
%                   (x(1)-1)^2 + (x(2)-1)^2 <= 25, 
%                   -5 <= x(1) <= 5 and
%                   -5 <= x(2) <= 5.
%
%      Specify the first constraint in a MATLAB file function such as
%         function [c,ceq] = mycon(x)
%         c = (x(1)-1)^2 + (x(2)-1)^2 - 25;
%         ceq = [];
%
%      Implement the typical workflow
%         opts = optimoptions2('fmincon2','Algorithm','sqp')
%         problem = createOptimProblem2('fmincon2','objective', ...
%         @(x) x(1)^2 + 4*sin(5*x(2)),'x0',[3 3],'lb',[-5 -5], ...
%         'ub',[5 5],'nonlcon',@mycon,'options',opts)
%         ms = MultiStart2
%         [x,f] = run(ms,problem,20)
%
%   See also GLOBALSEARCH2

%   Copyright 2009-2015 The MathWorks, Inc.

properties
%USEPARALLEL Perform calculation in parallel mode
%   The UseParallel property determines whether MultiStart2 will distribute
%   the calls to the Optimization Toolbox solver among available
%   processors. When UseParallel is set to true (logical) and MATLAB sessions for
%   parallel computation are available (e.g. via the PARPOOL command),
%   the calls to the Optimization toolbox solver in MultiStart2 will run in
%   parallel. Any requests for parallellization inside the Optimization
%   Toolbox solver will be ignored. 
%
%   See also MULTISTART2, RUN  
    UseParallel = false;

end
   methods
       function ms = MultiStart2(varargin)
%MultiStart2 Construct a multi-start based global optimization solver
%
%   MS = MultiStart2 constructs a new multi-start optimization solver
%   with its properties set to defaults.
%
%   MS = MultiStart2(PROP, VAL, ...) specifies a set of property-value
%   pairs that are applied to the multi-start optimization solver before
%   creating it.
%
%   MS = MultiStart2(OLDMS, PROP, VAL, ...) creates a copy of the
%   MultiStart2 solver OLDMS with the named properties altered with the
%   specified values.
%
%   MS = MultiStart2(GS) constructs a new MultiStart2 solver and copies
%   the common parameter values in the Global Search2 solver GS into the new
%   solver MS.
%
%   The MultiStart2 solver provides the means to call a specified 
%   Optimization Toolbox solver with multiple start points.
%  
%   See also MULTISTART2, RUN, GLOBALSEARCH2

           ms = ms@AbstractGlobalSolver2(varargin{:});
           % Record the version of the class
           ms.Version = 2;
       end
       function obj = set.UseParallel(obj,value)
           % error check for UseParallel property; throw error; convert
           % values to logical, if needed.
           obj.UseParallel = validateopts_UseParallel(value,true,true);
       end
       function [x,fval,exitflag,output,solutions] = ...
               run(obj,problem,startPointSets)
%RUN Run a local solver from multiple start points for a given optimization
%problem
%   RUN accepts Optimization Toolbox functions fmincon2, fminunc, lsqnonlin
%   and lsqcurvefit2 as local solvers.
%
%   X = RUN(MS,PROBLEM,NUMPOINTS) instructs the multi-start solver, MS, to
%   call the local solver specified in PROBLEM.solver from PROBLEM.x0 and
%   (NUMPOINTS - 1) random start points. PROBLEM is an optimization problem
%   structure that can be created by the function CREATEOPTIMPROBLEM2.
%   PROBLEM.solver can be one of the following: 'fmincon2', 'fminunc',
%   'lsqnonlin' or 'lsqcurvefit2'. RUN returns the point, X, that achieved
%   the lowest objective function value (if available). NUMPOINTS is a
%   positive integer.
%
%   X = RUN(MS,PROBLEM,STARTPOINTSETS) calls the specified Optimization
%   Toolbox solver from start points provided in STARTPOINTSETS.
%   STARTPOINTSETS is a start point set object (e.g. RANDOMSTARTPOINTSET2)
%   or a cell array of these objects. See the help for RANDOMSTARTPOINTSET2
%   for more information on how to create these and other start point sets.
%
%   [X,FVAL] = RUN(MS,PROBLEM,...) returns the value of the objective
%   function for PROBLEM at the solution X.
%
%   [X,FVAL,EXITFLAG] = RUN(MS,PROBLEM,...) returns an EXITFLAG that
%   describes the exit condition of multi-start algorithm. Possible values
%   of EXITFLAG and the corresponding exit conditions are listed below.
%   
%     2  At least one local minimum located. Some runs of local solver
%     converged.
%     1  At least one local minimum located. All runs of local solver
%     converged.
%     0  No local minimum located. Local solver called at least once and at
%     least one local minima call ran out of iterations.
%    -1  Stopped by the output or plot function.
%    -2  No feasible local minimum found.
%    -5  MaxTime limit exceeded.
%    -8  No solution found.
%   -10  Failures in the user provided functions encountered.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = RUN(MS,PROBLEM,...) also returns a
%   structure OUTPUT with the number of function evaluations in
%   OUTPUT.funcCount, the total number of local solver calls in
%   OUTPUT.localSolverTotal, the number of local solver runs with positive
%   exitflag in OUTPUT.localSolverSuccess, with zero exitflag 
%   OUTPUT.localSolverIncomplete, with negative exitflag
%   OUTPUT.localSolverNoSolution and the exit message in OUTPUT.message.
%
%   [X,FVAL,EXITFLAG,OUTPUT,SOLUTIONS] = RUN(MS,PROBLEM,...) also
%   returns a vector of solutions, SOLUTIONS. This contains output from the
%   Optimization Toolbox solver for each distinct local minimum found by
%   RUN. The vector SOLUTIONS is sorted by the objective function values.
%
%   See also MULTISTART2, CREATEOPTIMPROBLEM2, RANDOMSTARTPOINTSET2, CUSTOMSTARTPOINTSET2

       % MultiStart2.run requires exactly three input arguments    
       if nargin ~= 3
           error(message('globaloptim:MultiStart2:run:InvalidNumInputArgs'))
       end
       
       % First argument must be a MultiStart2 object
       if ~isa(obj, 'MultiStart2')
           error(message('globaloptim:MultiStart2:run:FirstArgNotObject'));
       end
       
       % MultiStart2 object must be scalar. This check also stops empty
       % objects being passed to the run method.
       if ~isscalar(obj)
           error(message('globaloptim:MultiStart2:run:ObjectNotScalar'));
       end
        
       % Check whether problem is a valid optim structure, has the x0 field
       % and employs the valid solvers
       requiredFields = {'solver','options','x0'};
       validValues = {{}, {}, {}};
       validSolvers = {'fmincon2','fminunc','lsqnonlin','lsqcurvefit2'};
       obj.checkProblem(problem,requiredFields,validValues,validSolvers);       
       dimX0 = numel(problem.x0);
       
       if isnumeric(startPointSets) && isreal(startPointSets) && ...
               isscalar(startPointSets) && (startPointSets >= 1) && ...
               startPointSets == floor(startPointSets)
           % If a positive integer is provided use one less number of
           % points in RandomStartPointSet2, since x0 is considered as a
           % start point.           
           numPoints = startPointSets - 1;
           if numPoints == 0
               % User wants to run a single point: that will be x0.
               startPointSets = {CustomStartPointSet2(problem.x0(:)')};
           else
               startPointSets = {CustomStartPointSet2(problem.x0(:)'), ...
                   RandomStartPointSet2('NumStartPoints',numPoints)};
           end
       else  
           % Check if the startPointSets has a single StartPointSet or a cell
           % array of sets
           if ~iscell(startPointSets) && ~isa(startPointSets,'AbstractStartPointSet2')
               error(message('globaloptim:MultiStart2:run:InvalidStartPointSetsInput'));
           end
           if isa(startPointSets,'AbstractStartPointSet2')
               startPointSets = {startPointSets};
           end        
           % Check whether the dimensions agree among themselves and with
           % problem
           for i = 1:length(startPointSets)
               if isa(startPointSets{i},'CustomStartPointSet2')
                   if startPointSets{i}.DimStartPoints ~= dimX0
                       error(message('globaloptim:MultiStart2:run:InvalidDimStartPointSet'));
                   end
               elseif ~isa(startPointSets{i},'AbstractStartPointSet2')
                        error(message('globaloptim:MultiStart2:run:NotAStartPointSet'));                 
               end
           end
       end       
            
       % Override local solver display if not set by user.
       problem = obj.overrideLocalSolverDisplay(problem);
       
       % For local solver fminunc: if options are passed in as an optimoptions2 object, 
       % check if there are warnings to display related to fminunc's algorithm default 
       % transition. (We create the warnings here before the options are converted into a
       % structure.)
       if isa(problem.options, 'optim.options.SolverOptions2') && strcmp(problem.solver,'fminunc') 
           [transitionMsgID,transitionMsgTxt] = algorithmDefaultTransitionMessages(problem.options);
           if ~isempty(transitionMsgID)
               warning(transitionMsgID,transitionMsgTxt)
           end
       end
       
       % Extract the local options structure from problem, if local
       % options is specified as a SolverOptions2 object
       if isa(problem.options, 'optim.options.SolverOptions2')
           problem.options = extractOptionsStructure(problem.options);
       end
       
       % map obj properties to msoptions
       msoptions.Display = obj.Display;
       msoptions.StartPointsToRun = obj.StartPointsToRun;
       msoptions.TolFun = obj.FunctionTolerance;
       msoptions.TolX = obj.XTolerance;
       msoptions.MaxTime = obj.MaxTime;
       msoptions.UseParallel = obj.UseParallel;
       msoptions.OutputFcns = obj.OutputFcn;
       msoptions.PlotFcns = obj.PlotFcn;
       % Call the driver
       if nargout == 5
           [x,fval,exitflag,output,solutions] = ...
               fmultistart2(problem,startPointSets,msoptions);
       else
           [x,fval,exitflag,output] = ...
               fmultistart2(problem,startPointSets,msoptions);
       end
       end
   end
   methods (Static)
       function obj = loadobj(obj)
           % Version 2 has two new properties OutputFcn and PlotFcn that
           % are set to empty by default.
           obj.Version = 2;
       end
   end
end % MultiStart2 class definition


