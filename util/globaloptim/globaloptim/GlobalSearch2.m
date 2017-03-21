classdef (Sealed) GlobalSearch2 < AbstractGlobalSolver2
%GlobalSearch2 A scatter-search2 based global optimization solver. 
%   A GlobalSearch2 object presents a solver that attempts to locate the
%   solution with lowest objective function value. GlobalSearch2 solver
%   first generates trial points employing a Scatter Search2 method. These
%   trial points are then filtered and fmincon2 is started from each of the
%   filtered points.
%
%   GS = GlobalSearch2 constructs a new global search2 optimization solver
%   with its properties set to defaults.
%
%   GS = GlobalSearch2(PROP, VAL, ...) specifies a set of property-value
%   pairs that are applied to the global search2 optimization solver before
%   creating it.
%
%   GS = GlobalSearch2(OLDGS, PROP, VAL, ...) creates a copy of the
%   GlobalSearch2 solver OLDGS. GS will have the named properties altered
%   with the specified values.
%
%   GS = GlobalSearch2(MS) constructs a new GlobalSearch2 solver and copies
%   the common parameter values in the multi-start solver MS into the new
%   solver GS.
%
%   GlobalSearch2 properties:
%       Display                 - content level of display
%       FunctionTolerance       - minimum distance between two separate 
%                                 objective function values
%       XTolerance              - minimum distance between two separate 
%                                 points 
%       MaxTime                 - time allowed to run the solver
%       StartPointsToRun        - additional filter for Stage 2 start 
%                                 points     
%       OutputFcn               - user-defined output functions
%       PlotFcn                 - functions to plot solver progress
%       NumTrialPoints          - number of trial points analyzed
%       BasinRadiusFactor       - factor to determine the amount of 
%                                 decrease in the basin radius
%       DistanceThresholdFactor - factor to tune the effect of the distance
%                                 filter
%       MaxWaitCycle            - maximum number of consecutive points 
%                                 analyzed before a threshold increase or 
%                                 a basin of attraction decrease
%       NumStageOnePoints       - number of trial points analyzed in Stage 1                               
%       PenaltyThresholdFactor  - factor to determine the update amount in 
%                                 the threshold for rejecting trial points                       
%
%   GlobalSearch2 method:
%       run             - search2 for the best solution of a given
%                         optimization problem
%
%   Typical workflow to run the GlobalSearch2 solver:
%   ==============================================
%   1. Set up the PROBLEM structure
%       PROBLEM = createOptimProblem2('fmincon2','objective',...)
%   2. Construct the GlobalSearch2 solver
%       GS = GlobalSearch2
%   3. Run the solver 
%       run(GS,PROBLEM)
%
%   Example:
%      Run global search2 on the optimization problem 
%         minimize peaks(x, y); subject to 
%                   (x+3)^2 + (y+3)^2 <= 36,
%                   -3 <= x <= 3 and  
%                   -3 <= y <= 3.
%
%      Specify the first constraint in a MATLAB file function such as
%         function [c,ceq] = mycon(x)
%         c = (x(1)+3)^2 + (x(2)+3)^2 - 36;
%         ceq = [];
%
%      Implement the typical workflow
%         problem = createOptimProblem2('fmincon2','objective', ...
%         @(x) peaks(x(1), x(2)), 'x0', [1 2], 'lb', [-3 -3], ...
%         'ub', [3 3], 'nonlcon', @mycon)
%         gs = GlobalSearch2
%         [x, f] = run(gs, problem)
%
%   See also MULTISTART2

%   Copyright 2009-2015 The MathWorks, Inc.

   properties
%NUMTRIALPOINTS Number of trial points analyzed 
%   The NumTrialPoints property sets the number of trial points that are 
%   generated by the Scatter Search2 portion of the algorithm. Once all the 
%   trial points have been analyzed by global search2, it will terminate.
%
%   See also GLOBALSEARCH2, RUN, NUMSTAGEONEPOINTS       
       NumTrialPoints = 1000;
       
%BASINRADIUSFACTOR Factor to determine the amount of decrease in the basin radius
%   The BasinRadiusFactor property determines the factor that a given basin
%   of attraction radius will be decreased by when MaxWaitCycle
%   consecutive trial points lie in that basin. 
%
%   See also GLOBALSEARCH2, RUN, MAXWAITCYCLE
      BasinRadiusFactor = 0.2;
      
%DISTANCETHRESHOLDFACTOR Factor to tune the effect of the distance filter
%   The DistanceThresholdFactor property provides the parameter to tune the
%   effect of the distance filter. In Stage 2 of the global search2
%   algorithm, a trial point must pass three tests before a local solver is
%   started from that point. In the distance test, a local solver may only
%   be started from a given trial point, T, if the following holds:
%
%    For each local solution already located, L(i):
%        Distance between T and L(i) > DistanceThresholdFactor*maxdist(i)
%
%   maxdist(i) is an estimate of the basin radius for the i-th local
%   solution and is internally calculated.
%   As DistanceThresholdFactor approaches zero the distance filter effect
%   vanishes. This is appropriate if there are many closely spaced local
%   solutions. Increasing DistanceThresholdFactor increases the effect of
%   distance filtering, namely that the local solver will be called from
%   fewer points in Stage 2.
%
%   See also GLOBALSEARCH2, RUN
      DistanceThresholdFactor = 0.75;
      
%MAXWAITCYCLE Maximum number of consecutive points analyzed before a  
%             threshold increase or a basin of attraction decrease
%   The MaxWaitCycle property sets a maximum number of consecutive trial
%   points for the following tests:
%
%   1) A trial point can be rejected (i.e. a local solver is not run from
%   the point) because the exact penalty function value at this point is
%   greater than a threshold. If MaxWaitCycle consecutive trial points are
%   rejected in this way then the threshold value will be updated with the
%   formula:
%       threshold = threshold + PenaltyThresholdFactor*(1.0+abs(threshold))
%
%   2) The radius of a basin of attraction of a located minimum is reduced
%   by the solver if MaxWaitCycle consecutive trial points lie in the
%   basin. The following formula is used to update the radius of the basin
%   of attraction, maxdist
%       maxdist = maxdist(i)*(1-BasinRadiusFactor)
%
%   See also GLOBALSEARCH2, RUN, PENALTYTHRESHOLDFACTOR, BASINRADIUSFACTOR
      MaxWaitCycle = 20;
      
%NUMSTAGEONEPOINTS Number of trial points analyzed in Stage 1 
%   The NumStageOnePoints property is a positive integer that sets the 
%   number of trial points analyzed by the scatter search2 in Stage 1 of 
%   the algorithm. After analyzing this many points the best quality point 
%   is chosen and fmincon2 is run from that point.
%
%   See also GLOBALSEARCH2, RUN
      NumStageOnePoints = 200;
      
%PENALTYTHRESHOLDFACTOR Factor to determine the update amount in the  
%                       threshold for rejecting trial points
%   The PenaltyThresholdFactor property is used to update the threshold
%   value in the formula:
%       threshold = threshold + PenaltyThresholdFactor*(1.0+abs(threshold))
%   A trial point can be rejected because the exact penalty function value
%   at this point is greater than this threshold. If MaxWaitCycle
%   consecutive trial points are rejected in this way then the threshold
%   value will be updated.
%
%   See also GLOBALSEARCH2, RUN, MAXWAITCYCLE
      PenaltyThresholdFactor = 0.2;    
            
   end
   methods
       function gs = GlobalSearch2(varargin)
%GlobalSearch2 Construct a scatter-search2 based global optimization solver
%
%   GS = GlobalSearch2 constructs a new global search2 optimization solver
%   with its properties set to defaults.
%
%   GS = GlobalSearch2(PROP, VAL, ...) specifies a set of property-value
%   pairs that are applied to the global search2 optimization solver before
%   creating it.
%
%   GS = GlobalSearch2(OLDGS, PROP, VAL, ...) creates a copy of the
%   GlobalSearch2 solver OLDGS. GS will have the named properties altered
%   with the specified values.
%
%   GS = GlobalSearch2(MS) constructs a new GlobalSearch2 solver and copies
%   the common parameter values in the multi-start solver MS into the new
%   solver GS.
%   
%   A GlobalSearch2 solver attempts to locate the solution of a supplied
%   problem with the lowest objective function value. GlobalSearch2 first
%   generates trial points via a Scatter Search2 method. These trial points
%   are then filtered and fmincon2 is started from each of the filtered
%   points.
%
%   See also GLOBALSEARCH2, RUN, MULTISTART2
           
           % Check to see if user is trying to set 'UseParallel' for
           % GlobalSearch2. If so, error and point them towards MultiStart2
           if any(strcmp(varargin, 'UseParallel'))
               error(message('globaloptim:GlobalSearch2:InvalidArgument'));
           end
           
           % Call superclass constructor
           gs = gs@AbstractGlobalSolver2(varargin{:});
           
           % Record the version of the class
           gs.Version = 2;          
       end
       function obj = set.NumTrialPoints(obj,value)
           % error check for NumTrialPoints property
           typeValueChecker2('posInteger',value,'NumTrialPoints');
           obj.NumTrialPoints = value;
       end       
       function obj = set.BasinRadiusFactor(obj,value)
           % error check for BasinRadiusFactor property
           typeValueChecker2('boundedReal',value,'BasinRadiusFactor',[0 1]);
           obj.BasinRadiusFactor = value;
       end
       function obj = set.DistanceThresholdFactor(obj,value)
           % error check for DistanceThresholdFactor property
           typeValueChecker2('nonNegReal',value,'DistanceThresholdFactor');
           obj.DistanceThresholdFactor = value;
       end
       function obj = set.MaxWaitCycle(obj,value)
           % error check for MaxWaitCycle property
           typeValueChecker2('posInteger',value,'MaxWaitCycle');
           obj.MaxWaitCycle = value;
       end
       function obj = set.NumStageOnePoints(obj,value)
           % error check for NumStageOnePoints property
           typeValueChecker2('posInteger',value,'NumStageOnePoints');
           obj.NumStageOnePoints = value;
       end
       function obj = set.PenaltyThresholdFactor(obj,value)
           % error check for PenaltyThresholdFactor property
           typeValueChecker2('nonNegReal',value,'PenaltyThresholdFactor');
           obj.PenaltyThresholdFactor = value;
       end
       function [x,fval,exitflag,output,solutionSet] = run(obj,problem)
%RUN Search2 for the best solution of a given optimization problem
%   RUN attempts to find the best solution for the given problem by a
%   scatter search2 based procedure. The required input argument problem is
%   an optimization problem structure.
%
%   X = RUN(OBJ,PROBLEM) returns the best point, X, that achieved the most
%   minimum objective function value for the optimization problem described
%   in PROBLEM. PROBLEM should be valid problem structure for the
%   Optimization Toolbox solvers fminunc or fmincon2. 
%
%   [X,FVAL] = RUN(OBJ,PROBLEM) returns the value of the objective
%   function for PROBLEM at the solution X.
%
%   [X,FVAL,EXITFLAG] = RUN(OBJ,PROBLEM) returns an EXITFLAG that describes
%   the exit condition of scatter-search2 based algorithm. Possible values
%   of EXITFLAG and the corresponding exit conditions are listed below.
%   
%     2  At least one local minimum located. One or more runs of the local
%     solver converged with a positive local solver exit flag.
%     1  At least one local minimum located. All runs of local solver
%     converged with a positive local solver exit flag.
%     0  No local minima located. Local solver called at least once and at
%     least one local minima call ran out of iterations.
%    -1  Stopped by the output or plot function.
%    -2  No feasible solution found in all runs of the local solver.
%    -5  MaxTime limit exceeded.
%    -8  No solution found.
%   -10  Failures encountered in the user provided functions.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = RUN(OBJ,PROBLEM) also returns a
%   structure OUTPUT with the number of trial points analyzed in
%   OUTPUT.numTrialPoints, the number of function evaluations in
%   OUTPUT.funcCount, the total number of local solver runs in
%   OUTPUT.localSolverTotal, the number of local solver runs with positive
%   exitflag in OUTPUT.localSolverSuccess, with zero exitflag 
%   OUTPUT.localSolverIncomplete, with negative exitflag
%   OUTPUT.localSolverNoSolution and the exit message in OUTPUT.message.
%
%   [X,FVAL,EXITFLAG,OUTPUT,SOLUTIONS] = RUN(OBJ,PROBLEM) also returns a
%   vector of solutions, SOLUTIONS. This contains output from the
%   Optimization Toolbox solver for each distinct local minimum found by
%   RUN. The vector SOLUTIONS is sorted by the objective function values.
%
%   See also GLOBALSEARCH2

        % Ensure two arguments are specified 
        if nargin ~= 2
            error(message('globaloptim:GlobalSearch2:run:InvalidNumInputArgs'));
        end

        % First argument must be a GlobalSearch2 object
        if ~isa(obj, 'GlobalSearch2')
           error(message('globaloptim:GlobalSearch2:run:FirstArgNotObject'));
        end
        
        % GlobalSearch2 object must be scalar. This check also stops empty
        % objects being passed to the run method.
        if ~isscalar(obj)
           error(message('globaloptim:GlobalSearch2:run:ObjectNotScalar'));
        end
        
        % Ensure that NumTrialPoints >= NumStageOnePoints
        if obj.NumStageOnePoints > obj.NumTrialPoints
            error(message('globaloptim:GlobalSearch2:run:InvalidNumTrialPoints'));
        end
                
        % Check the supplied problem structure
        probRequiredFields = {'solver', 'options', 'x0'};
        probValidValues = {'','',''};
        validSolvers = {'fmincon2'};
        obj.checkProblem(problem, probRequiredFields, probValidValues, validSolvers);

        % Override local solver display if not set by user.
        problem = obj.overrideLocalSolverDisplay(problem);
                
        % Extract individual components of problem for globalsearchnlp2
        [FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,localOptions] = separateOptimStruct(problem);

        % Extract the local options structure from problem, if local
        % options is specified as a SolverOptions2 object
        if isa(localOptions, 'optim.options.SolverOptions2')
            localOptions = extractOptionsStructure(localOptions);
        end
        
        % Map obj properties to options
        options.MaxIter = obj.NumTrialPoints;
        options.BasinRadiusFactor = obj.BasinRadiusFactor;
        options.DistanceThresholdFactor = obj.DistanceThresholdFactor;
        options.MaxWaitCycle = obj.MaxWaitCycle;
        options.StageOneIter = obj.NumStageOnePoints;
        options.PenaltyThresholdFactor = obj.PenaltyThresholdFactor;
        options.Display = obj.Display;
        options.StartPointsToRun = obj.StartPointsToRun;
        options.TolFun = obj.FunctionTolerance;
        options.TolX = obj.XTolerance;
        options.MaxTime = obj.MaxTime;
        options.OutputFcns = obj.OutputFcn;
        options.PlotFcns = obj.PlotFcn;

        % Call the driver
        if nargout == 5
            [x,fval,exitflag,output,solutionSet] = ...
                globalsearchnlp2(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,options,localOptions);
        else
            [x,fval,exitflag,output] = ...
                globalsearchnlp2(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,options,localOptions);
        end
       end
   end
   methods (Static)
       function obj = loadobj(obj)
           % Version 2 has two new properties OutputFcns and PlotFcns that
           % are set to empty by default.
           obj.Version = 2;
       end
   end
end

function varargout = separateOptimStruct(myStruct)
% SEPARATEOPTIMSTRUCT takes a problem structure and returns individual fields of
%   the structure to the caller. The caller (always a solver) information is
%   found from the 'solver' field in 'myStruct'.
%   This function has the same code as matlab\funfun\separateOptimStruct.m
%   except the caller-solver check.
%
%   Note that this function assumes that myStruct is a valid problem
%   structure.

% Extract solver and options from supplied structure.
solver = myStruct.solver;
options = myStruct.options;

% Call createProblemStruct to create a structure containing the full input
% argument list for the specified solver. Note that the second argument is
% required but can be []. 
myStruct = createProblemStruct(solver,[],myStruct); 

% Extract the full input argument list from the structure.
probFieldNames = fieldnames(myStruct);
numFieldNames = length(probFieldNames);
varargout = cell(1, numFieldNames);
for i = 1:numFieldNames - 1 % Last field is the 'solver' field
    varargout{i} = myStruct.(probFieldNames{i});
end

% The last input argument to the Optimization toolbox solvers is the
% options structure which we retrieved from the original structure before
% createProblemStruct removed it.
varargout{numFieldNames} = options; 

end