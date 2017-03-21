classdef (Sealed) PatternsearchOptions2 < optim.options.SingleAlgorithm2
    %PatternsearchOptions2 Options for PATTERNSEARCHOPTIONS2
    %
    %   The OPTIM.OPTIONS.PATTERNSEARCHOPTIONS2 class allows the user to create a set
    %   of options for the PATTERNSEARCHOPTIONS2 solver. For a list of options that can
    %   be set, see the documentation for patternsearch2.
    %
    %   OPTS = OPTIM.OPTIONS.PATTERNSEARCHOPTIONS2 creates a set of options for
    %   PATTERNSEARCHOPTIONS2 with the options set to their default values.
    %
    %   OPTS = OPTIM.OPTIONS.PATTERNSEARCHOPTIONS2(PARAM, VAL, ...) creates a set of
    %   options for PATTERNSEARCHOPTIONS2 with the named parameters altered with the
    %   specified values.
    %
    %   OPTS = OPTIM.OPTIONS.PATTERNSEARCHOPTIONS2(OLDOPTS, PARAM, VAL, ...) creates a
    %   copy of OLDOPTS with the named parameters altered with the specified
    %   values.
    %
    %   See also OPTIM.OPTIONS.SingleAlgorithm2, OPTIM.OPTIONS.SolverOptions2
    
    %   Copyright 2015-2015 The MathWorks, Inc.
    
    properties (Dependent)
        %ACCELERATEMESH Accelerate convergence near a minimum
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        AccelerateMesh
        
        %CONSTRAINTTOLERANCE Tolerance on the constraints
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        ConstraintTolerance
        
        %DISPLAY Level of display
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        Display
        
        %FUNCTIONTOLERANCE Tolerance on the function value
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        FunctionTolerance
        
        %INITIALMESHSIZE Initial mesh size for pattern algorithm
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        InitialMeshSize
        
        %MAXFUNCTIONEVALUATIONS Maximum number of objective function evaluations
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        MaxFunctionEvaluations
        
        %MAXITERATIONS Maximum number of iterations
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        MaxIterations
        
        %MAXTIME Total time (in seconds) allowed for optimization
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        MaxTime
        
        %MESHCONTRACTIONFACTOR Mesh contraction factor, used when iteration is unsuccessful
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        MeshContractionFactor
        
        %MESHEXPANSIONFACTOR Mesh expansion factor, used when iteration is successful
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        MeshExpansionFactor
        
        %MESHTOLERANCE Tolerance on the mesh size
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        MeshTolerance
        
        %OUTPUTFCN Functions that get iterative data and can change options
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        OutputFcn
        
        %PLOTFCN Plots various measures of progress while the algorithm executes
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        PlotFcn
        
        %POLLMETHOD Polling strategy used in pattern search2
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        PollMethod
        
        %POLLORDERALGORITHM Algorithm selecting order of poll2 directions in pattern search2
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        PollOrderAlgorithm
        
        %SCALEMESH Automatic scaling of variables
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        ScaleMesh
        
        %SEARCHFCN Function that performs the search2 phase of patternsearch2
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        SearchFcn
        
        %STEPTOLERANCE Tolerance on the change in position and mesh size
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        StepTolerance
        
        %USECOMPLETEPOLL Perform complete poll2 around current iterate
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        UseCompletePoll
        
        %USECOMPLETESEARCH Perform complete search2 around current iterate when the search2 method is a poll2 method
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        UseCompleteSearch
        
        %USEPARALLEL Compute the objective function of particles in parallel
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        UseParallel
        
        %USEVECTORIZED Compute the objective function with a vectorized function call
        %
        %   For more information, type "doc patternsearch2" and see the "Options"
        %   section in the PATTERNSEARCHOPTIONS2 documentation page.
        UseVectorized
    end
    
    % Hidden properties
    properties (Hidden, Dependent)
        
        %CACHE Keep history of points polled
        %
        Cache
        
        %CACHESIZE Size of the cache history
        %
        CacheSize
        
        %CACHETOL Tolerance on points in the cache
        %
        CacheTol
        
        %COMPLETEPOLL Perform complete poll2 around current iterate
        %
        CompletePoll

        %COMPLETESEARCH Perform complete search2 around current iterate when
        %               the search2 method is a poll2 method
        %
        CompleteSearch
        
        %INITIALPENALTY Initial value of the penalty parameter
        %
        InitialPenalty

        %MAXMESHSIZE Maximum number of objective function evaluations
        %
        MaxFunEvals

        %MAXMESHSIZE Maximum number of iterations
        %
        MaxIter
        
        %MAXMESHSIZE Maximum mesh size used in a poll2/search2 step
        %
        MaxMeshSize

        %MESHACCELERATOR Accelerate convergence near a minimum
        %
        MeshAccelerator

        %MESHCONTRACTION Mesh contraction factor, used when iteration is unsuccessful
        %
        MeshContraction

        %MESHEXPANSIONFACTOR Mesh expansion factor, used when iteration is successful
        %
        MeshExpansion
        
        %MESHROTATE Rotate the pattern before declaring a point to be optimum
        %
        MeshRotate
        
        %OUTPUTFCNS Functions that get iterative data and can change options
        %
        OutputFcns
        
        %PLOTFCNS Plots various measures of progress while the algorithm executes
        %
        PlotFcns
        
        %PENALTYFACTOR  Penalty update parameter
        %
        PenaltyFactor
        
        %PLOTINTERVAL  Specifies that plot functions will be called at every interval
        %
        PlotInterval
        
        %POLLINGORDER Algorithm selecting order of poll2 directions in pattern search2
        %
        PollingOrder        
        
        %SEARCHMETHOD Function that performs the search2 phase of patternsearch2
        %
        SearchMethod

        %TIMELIMIT Total time (in seconds) allowed for optimization
        %
        TimeLimit
        
        %TOLBIND  Binding tolerance
        %
        TolBind
        
        %TOLBIND  Constraint tolerance
        %        
        TolCon
        
        %TOLFUN  Function tolerance
        %                
        TolFun
        
        %TOLMESH  Mesh tolerance
        %                        
        TolMesh
        
        %TOLX  Step tolerance
        %                                
        TolX        
        
        %USEVECTORIZED Compute the objective function with a vectorized function call
        %
        Vectorized
        
    end
    
    properties (Hidden, Access = protected)
        
        %OPTIONSSTORE Contains the option values and meta-data for the class
        %
        OptionsStore = createOptionsStore;
    end
    
    properties (Hidden)
        
        %SOLVERNAME Name of the solver that the options are intended for
        %
        SolverName = 'patternsearch2';
    end
    
    properties (SetAccess = private, GetAccess = private)
        
        % Version number for the optim.options.PatternsearchOptions2 class.
        % We do not change the Version property (base class).
        PatternsearchVersion = 1
    end
    
    % -------------------------------------------------------
    
    methods (Hidden)
        
        function obj = PatternsearchOptions2(varargin)
            %PatternsearchOptions2 Options for PATTERNSEARCHOPTIONS2
            %
            %   The OPTIM.OPTIONS.PATTERNSEARCHOPTIONS2 class allows the user to create a set
            %   of options for the PATTERNSEARCHOPTIONS2 solver. For a list of options that can
            %   be set, see the documentation for PATTERNSEARCHOPTIONS2.
            %
            %   OPTS = OPTIM.OPTIONS.PATTERNSEARCHOPTIONS2 creates a set of options for
            %   PATTERNSEARCHOPTIONS2 with the options set to their default values.
            %
            %   OPTS = OPTIM.OPTIONS.PATTERNSEARCHOPTIONS2(PARAM, VAL, ...) creates a set of
            %   options for PATTERNSEARCHOPTIONS2 with the named parameters altered with the
            %   specified values.
            %
            %   OPTS = OPTIM.OPTIONS.PATTERNSEARCHOPTIONS2(OLDOPTS, PARAM, VAL, ...) creates a
            %   copy of OLDOPTS with the named parameters altered with the specified
            %   values.
            %
            %   See also OPTIM.OPTIONS.SingleAlgorithm2, OPTIM.OPTIONS.SolverOptions2
            
            % Call the superclass constructor
            obj = obj@optim.options.SingleAlgorithm2(varargin{:});
            
        end
        
    end
    
    % Set/get methods
    methods
        function obj = set.AccelerateMesh(obj, value)
            obj = setNewProperty(obj, 'AccelerateMesh', value);
        end
        
        function obj = set.ConstraintTolerance(obj, value)
            obj = setAliasProperty(obj, 'ConstraintTolerance', 'TolCon', value);
        end
        
        function obj = set.Display(obj, value)
            obj = setProperty(obj, 'Display', value, ...
                {'none','off','iter','diagnose2','final'});
        end
        
        function obj = set.FunctionTolerance(obj, value)
            obj = setAliasProperty(obj, 'FunctionTolerance', 'TolFunValue', value);
        end
        
        function obj = set.InitialMeshSize(obj, value)
            obj = setProperty(obj, 'InitialMeshSize', value);
        end
        
        function obj = set.MaxFunctionEvaluations(obj, value)
            obj = setAliasProperty(obj, 'MaxFunctionEvaluations', 'MaxFunEvals', value);
        end
        
        function obj = set.MaxIterations(obj, value)
            obj = setAliasProperty(obj, 'MaxIterations', 'MaxIter', value);
        end
        
        function obj = set.MaxTime(obj, value)
            obj = setAliasProperty(obj, 'MaxTime', 'TimeLimit', value);
        end
        
        function obj = set.MeshContractionFactor(obj, value)
            obj = setAliasProperty(obj, 'MeshContractionFactor', 'MeshContraction', value);
        end
        
        function obj = set.MeshExpansionFactor(obj, value)
            obj = setAliasProperty(obj, 'MeshExpansionFactor', 'MeshExpansion', value);
        end
        
        function obj = set.MeshTolerance(obj, value)
            obj = setAliasProperty(obj, 'MeshTolerance', 'TolMesh', value);
        end
        
        function obj = set.OutputFcn(obj, value)
            obj = setAliasProperty(obj, 'OutputFcn', 'OutputFcns', value);
        end
        
        function obj = set.PlotFcn(obj, value)
            obj = setAliasProperty(obj, 'PlotFcn', 'PlotFcns', value);
        end
        
        function obj = set.PollMethod(obj, value)
            obj = setProperty(obj, 'PollMethod', value);
        end
        
        function obj = set.PollOrderAlgorithm(obj, value)
            obj = setAliasProperty(obj, 'PollOrderAlgorithm', 'PollingOrder', value);
        end
        
        function obj = set.ScaleMesh(obj, value)
            % ScaleMesh for patternsearch2 still accepts the undocumented
            % values, 'on'/'off'
            if any(strcmpi(value, {'on', 'off'}))
                value = optim.options.OptionAliasStore2.convertToLogical(value, 'on');
            end
            obj = setNewProperty(obj, 'ScaleMesh', value);
        end
        
        function obj = set.SearchFcn(obj, value)
            obj = setAliasProperty(obj, 'SearchFcn', 'SearchMethod', value);
        end
        
        function obj = set.StepTolerance(obj, value)
            obj = setAliasProperty(obj, 'StepTolerance', 'TolX', value);
        end
        
        function obj = set.UseCompletePoll(obj, value)
            obj = setNewProperty(obj, 'UseCompletePoll', value);
        end
        
        function obj = set.UseCompleteSearch(obj, value)
            obj = setNewProperty(obj, 'UseCompleteSearch', value);
        end
        
        function obj = set.UseParallel(obj, value)
            obj = setProperty(obj, 'UseParallel', value);
        end
        
        function obj = set.UseVectorized(obj, value)
            obj = setNewProperty(obj, 'UseVectorized', value);
        end
        
        %---------------------- Get functions -----------------------------
        
        function value = get.AccelerateMesh(obj)
            value = optim.options.OptionAliasStore2.convertToLogical( ...
                obj.OptionsStore.Options.MeshAccelerator, 'on');
        end
        
        function value = get.ConstraintTolerance(obj)
            value = obj.OptionsStore.Options.TolCon;
        end
        
        function value = get.Display(obj)
            value = obj.OptionsStore.Options.Display;
        end
        
        function value = get.FunctionTolerance(obj)
            value = obj.OptionsStore.Options.TolFunValue;
        end
        
        function value = get.InitialMeshSize(obj)
            value = obj.OptionsStore.Options.InitialMeshSize;
        end
        
        function value = get.MaxFunctionEvaluations(obj)
            value = obj.OptionsStore.Options.MaxFunEvals;
        end
        
        function value = get.MaxIterations(obj)
            value = obj.OptionsStore.Options.MaxIter;
        end
        
        function value = get.MaxTime(obj)
            value = obj.OptionsStore.Options.TimeLimit;
        end
        
        function value = get.MeshContractionFactor(obj)
            value = obj.OptionsStore.Options.MeshContraction;
        end
        
        function value = get.MeshExpansionFactor(obj)
            value = obj.OptionsStore.Options.MeshExpansion;
        end
        
        function value = get.MeshTolerance(obj)
            value = obj.OptionsStore.Options.TolMesh;
        end
        
        function value = get.OutputFcn(obj)
            value = obj.OptionsStore.Options.OutputFcns;
        end
        
        function value = get.PlotFcn(obj)
            value = obj.OptionsStore.Options.PlotFcns;
        end
        
        function value = get.PollMethod(obj)
            value = obj.OptionsStore.Options.PollMethod;
        end
        
        function value = get.PollOrderAlgorithm(obj)
            value = obj.OptionsStore.Options.PollingOrder;
        end
        
        function value = get.ScaleMesh(obj)
            value = optim.options.OptionAliasStore2.convertToLogical( ...
                obj.OptionsStore.Options.ScaleMesh, 'on');
        end
        
        function value = get.SearchFcn(obj)
            value = obj.OptionsStore.Options.SearchMethod;
        end
        
        function value = get.StepTolerance(obj)
            value = obj.OptionsStore.Options.TolX;
        end
        
        function value = get.UseCompletePoll(obj)
            value = optim.options.OptionAliasStore2.convertToLogical( ...
                obj.OptionsStore.Options.CompletePoll, 'on');
        end
        
        function value = get.UseCompleteSearch(obj)
            value = optim.options.OptionAliasStore2.convertToLogical( ...
                obj.OptionsStore.Options.CompleteSearch, 'on');
        end
        
        function value = get.UseParallel(obj)
            value = obj.OptionsStore.Options.UseParallel;
        end
        
        function value = get.UseVectorized(obj)
            value = optim.options.OptionAliasStore2.convertToLogical( ...
                obj.OptionsStore.Options.Vectorized, 'on');
        end
    end % get/set methods
    
    % Set/get methods for hidden options
    methods
        
        %----------------------- Set functions ----------------------------
        function obj = set.Cache(obj, value)
            obj = setProperty(obj, 'Cache', value);
        end

        function obj = set.CacheSize(obj, value)
            obj = setProperty(obj, 'CacheSize', value);
        end
        
        function obj = set.CacheTol(obj, value)
            obj = setProperty(obj, 'CacheTol', value);
        end
        
        function obj = set.CompletePoll(obj, value)
            obj = setProperty(obj, 'CompletePoll', value);
        end
        
        function obj = set.CompleteSearch(obj, value)
            obj = setProperty(obj, 'CompleteSearch', value);
        end
        
        function obj = set.InitialPenalty(obj, value)
            obj = setProperty(obj, 'InitialPenalty', value);
        end

        function obj = set.MaxFunEvals(obj, value)
            obj = setProperty(obj, 'MaxFunEvals', value);
        end

        function obj = set.MaxIter(obj, value)
            obj = setProperty(obj, 'MaxIter', value);
        end
        
        function obj = set.MaxMeshSize(obj, value)
            obj = setProperty(obj, 'MaxMeshSize', value);
        end

        function obj = set.MeshAccelerator(obj, value)
            obj = setProperty(obj, 'MeshAccelerator', value);
        end

        function obj = set.MeshContraction(obj, value)
            obj = setProperty(obj, 'MeshContraction', value);
        end
        
        function obj = set.MeshExpansion(obj, value)
            obj = setProperty(obj, 'MeshExpansion', value);
        end
        
        function obj = set.MeshRotate(obj, value)
            obj = setProperty(obj, 'MeshRotate', value);
        end

        function obj = set.OutputFcns(obj, value)
            obj = setProperty(obj, 'OutputFcns', value);
        end

        function obj = set.PenaltyFactor(obj, value)
            obj = setProperty(obj, 'PenaltyFactor', value);
        end
        
        function obj = set.PlotFcns(obj, value)
            obj = setProperty(obj, 'PlotFcns', value);
        end
                
        function obj = set.PlotInterval(obj, value)
            obj = setProperty(obj, 'PlotInterval', value);
        end

        function obj = set.PollingOrder(obj, value)
            obj = setProperty(obj, 'PollingOrder', value);
        end

        function obj = set.SearchMethod(obj, value)
            obj = setProperty(obj, 'SearchMethod', value);
        end

        function obj = set.TimeLimit(obj, value)
            obj = setProperty(obj, 'TimeLimit', value);
        end
        
        function obj = set.TolBind(obj, value)
            obj = setProperty(obj, 'TolBind', value);
        end
        
        function obj = set.TolCon(obj, value)
            obj = setProperty(obj, 'TolCon', value);
        end

        function obj = set.TolFun(obj, value)
            obj = setProperty(obj, 'TolFunValue', value);
        end
        
        function obj = set.TolMesh(obj, value)
            obj = setProperty(obj, 'TolMesh', value);
        end

        function obj = set.TolX(obj, value)
            obj = setProperty(obj, 'TolX', value);
        end

        function obj = set.Vectorized(obj, value)
            obj = setProperty(obj, 'Vectorized', value);
        end
        
        %----------------------- Get functions ----------------------------

        function value =  get.Cache(obj)
            value = obj.OptionsStore.Options.Cache;
        end

        function value = get.CacheSize(obj)
            value = obj.OptionsStore.Options.CacheSize;
        end
        
        function value = get.CacheTol(obj)
            value = obj.OptionsStore.Options.CacheTol;
        end

        function value = get.CompletePoll(obj)
            value = obj.OptionsStore.Options.CompletePoll;
        end
        
        function value = get.CompleteSearch(obj)
            value = obj.OptionsStore.Options.CompleteSearch;
        end
        
        function value = get.InitialPenalty(obj)
            value = obj.OptionsStore.Options.InitialPenalty;
        end

        function value = get.MaxFunEvals(obj)
            value = obj.OptionsStore.Options.MaxFunEvals;
        end

        function value = get.MaxIter(obj)
            value = obj.OptionsStore.Options.MaxIter;
        end
        
        function value = get.MaxMeshSize(obj)
            value = obj.OptionsStore.Options.MaxMeshSize;
        end

        function value = get.MeshAccelerator(obj)
            value = obj.OptionsStore.Options.MeshAccelerator;
        end

        function value = get.MeshContraction(obj)
            value = obj.OptionsStore.Options.MeshContraction;
        end

        function value = get.MeshExpansion(obj)
            value = obj.OptionsStore.Options.MeshExpansion;
        end
        
        function value = get.MeshRotate(obj)
            value = obj.OptionsStore.Options.MeshRotate;
        end
        
        function value = get.OutputFcns(obj)
            value = obj.OptionsStore.Options.OutputFcns;
        end
        
        function value = get.PenaltyFactor(obj)
            value = obj.OptionsStore.Options.PenaltyFactor;
        end

        function value = get.PlotFcns(obj)
            value = obj.OptionsStore.Options.PlotFcns;
        end
        
        function value = get.PlotInterval(obj)
            value = obj.OptionsStore.Options.PlotInterval;
        end

        function value = get.PollingOrder(obj)
            value = obj.OptionsStore.Options.PollingOrder;
        end

        function value = get.SearchMethod(obj)
            value = obj.OptionsStore.Options.SearchMethod;
        end

        function value = get.TimeLimit(obj)
            value = obj.OptionsStore.Options.TimeLimit;
        end
        
        function value = get.TolBind(obj)
            value = obj.OptionsStore.Options.TolBind;
        end

        function value = get.TolCon(obj)
            value = obj.OptionsStore.Options.TolCon;
        end

        function value = get.TolFun(obj)
            value = obj.OptionsStore.Options.TolFunValue;
        end

        function value = get.TolMesh(obj)
            value = obj.OptionsStore.Options.TolMesh;
        end

        function value = get.TolX(obj)
            value = obj.OptionsStore.Options.TolX;
        end

        function value = get.Vectorized(obj)
            value = obj.OptionsStore.Options.Vectorized;
        end
        
    end
    
    methods (Hidden)
        function OptionsStruct = extractOptionsStructure(obj)
            %EXTRACTOPTIONSSTRUCTURE Extract options structure from OptionsStore
            %
            %   OPTIONSSTRUCT = EXTRACTOPTIONSSTRUCTURE(OBJ) extracts a plain structure
            %   containing the options from obj.OptionsStore. The solver calls
            %   convertForSolver, which in turn calls this method to obtain a plain
            %   options structure.
            
            % Call superclass to generate structure with options set by the
            % user and all other fields empty
            OptionsStruct = getSetByUserOptionsStruct(obj);
            
            % Add TolFun to the options struct since patternsearch2 doesn't
            % recognize TolFunValue
            OptionsStruct.TolFun = OptionsStruct.TolFunValue;
        end
        
    end   % Hidden methods
    
end

function OS = createOptionsStore
%CREATEOPTIONSSTORE Create the OptionsStore
%
%   OS = createOptionsStore creates the OptionsStore structure. This
%   structure contains the options and meta-data for option display, e.g.
%   data determining whether an option has been set by the user. This
%   function is only called when the class is first instantiated to create
%   the OptionsStore structure in its default state. Subsequent
%   instantiations of this class pick up the default OptionsStore from the
%   MCOS class definition.
%
%   Class authors must create a structure containing all the options in a
%   field of OS called Defaults. This structure must then be passed to the
%   optim.options.generateSingleAlgorithmOptionsStore2 function to create
%   the full OptionsStore. See below for an example for PatternsearchOptions2.

% Define the option defaults for the solver
OS.Defaults.Cache = 'off';
OS.Defaults.CacheSize = 10000;
OS.Defaults.CacheTol = eps;
OS.Defaults.CompletePoll = 'off';
OS.Defaults.CompleteSearch = 'off';
OS.Defaults.Display = 'final';
OS.Defaults.InitialMeshSize = 1.0;
OS.Defaults.InitialPenalty = 10;
OS.Defaults.MaxFunEvals = '2000*numberOfVariables';
OS.Defaults.MaxIter = '100*numberOfVariables';
OS.Defaults.MaxMeshSize = Inf;
OS.Defaults.MeshAccelerator = 'off';
OS.Defaults.MeshContraction = 0.5;
OS.Defaults.MeshExpansion = 2.0;
OS.Defaults.MeshRotate = 'on';
OS.Defaults.OutputFcns = [];
OS.Defaults.PenaltyFactor = 100;
OS.Defaults.PlotFcns = [];
OS.Defaults.PlotInterval = 1;
OS.Defaults.PollMethod = 'GPSPositiveBasis2N';
OS.Defaults.PollingOrder = 'consecutive';
OS.Defaults.ScaleMesh = 'on';
OS.Defaults.SearchMethod = [];
OS.Defaults.TimeLimit = Inf;
OS.Defaults.TolBind = 1e-3;
OS.Defaults.TolCon = 1e-6;
OS.Defaults.TolFunValue = 1e-6;
OS.Defaults.TolMesh = 1e-6;
OS.Defaults.TolX = 1e-6;
OS.Defaults.UseParallel = false;
OS.Defaults.Vectorized = 'off';

% Call the package function to generate the OptionsStore
OS = optim.options.generateSingleAlgorithmOptionsStore2(OS);
end
