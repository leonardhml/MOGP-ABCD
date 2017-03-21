classdef (Sealed) Particleswarm2 < optim.options.SingleAlgorithm2
%Particleswarm2 Options for PARTICLESWARM2
%
%   The OPTIM.OPTIONS.PARTICLESWARM2 class allows the user to create a set
%   of options for the PARTICLESWARM2 solver. For a list of options that can
%   be set, see the documentation for particleswarm2.
%
%   OPTS = OPTIM.OPTIONS.PARTICLESWARM2 creates a set of options for
%   PARTICLESWARM2 with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.PARTICLESWARM2(PARAM, VAL, ...) creates a set of
%   options for PARTICLESWARM2 with the named parameters altered with the
%   specified values.
%
%   OPTS = OPTIM.OPTIONS.PARTICLESWARM2(OLDOPTS, PARAM, VAL, ...) creates a
%   copy of OLDOPTS with the named parameters altered with the specified
%   values.
%
%   See also OPTIM.OPTIONS.SingleAlgorithm2, OPTIM.OPTIONS.SolverOptions2

%   Copyright 2012-2015 The MathWorks, Inc.
    
    properties (Dependent)
%CREATIONFCN Handle to the function that creates the initial particles
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        CreationFcn
        
%DISPLAY Level of display
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        Display
        
%FUNCTIONTOLERANCE Termination tolerance on function value
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        FunctionTolerance
        
%HYBRIDFCN Run HybridFcn at the end of iterations of the solver
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        HybridFcn
        
%INERTIARANGE Two-element vector of the adaptive interia factor range
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        InertiaRange
        
%INITIALSWARMMATRIX Initial positions of swarm particles used to seed the
%algorithm
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        InitialSwarmMatrix
        
%INITIALSWARMSPAN Scalar or vector specifying the range of initial
%positions of the swarm
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        InitialSwarmSpan        
        
%MAXITERATIONS Maximum number of iterations allowed
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        MaxIterations
                
%MAXSTALLITERATIONS Maximum number of stalled iterations
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        MaxStallIterations
        
%MAXSTALLTIME Maximum time of stalled iterations
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        MaxStallTime
        
%MAXTIME The algorithm stop after running for MaxTime seconds
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        MaxTime
        
%MINNEIGHBORSFRACTION Minimum fraction of SwarmSize that will be considered
%as neighbors
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        MinNeighborsFraction
        
%OBJECTIVELIMIT Minimum objective function value desired
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        ObjectiveLimit
        
%OUTPUTFCN Functions that get iterative data and can change options
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        OutputFcn
        
%PLOTFCN Plots various measures of progress while the algorithm executes
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        PlotFcn
        
%SELFADJUSTMENTWEIGHT PARTICLESWARM2 algorithm adjustment amount due to
%particles individual best objective function value.
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        SelfAdjustmentWeight
        
%SOCIALADJUSTMENTWEIGHT PARTICLESWARM2 algorithm adjustment amount due to
%neighbors best objective function values.
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        SocialAdjustmentWeight
        
%SWARMSIZE Number of particles in the swarm
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        SwarmSize                
                        
%USEPARALLEL Compute the objective function of particles in parallel
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        UseParallel
        
%USEVECTORIZED Compute the objective function with a vectorized function call
%
%   For more information, type "doc particleswarm2" and see the "Options"
%   section in the PARTICLESWARM2 documentation page.
        UseVectorized
    end

    properties (Hidden, Access = protected)
        
%OPTIONSSTORE Contains the option values and meta-data for the class
%
        OptionsStore = createOptionsStore;
    end
    
    properties (Hidden)
        
%SOLVERNAME Name of the solver that the options are intended for
%
        SolverName = 'particleswarm2';
    end
    
    properties (Hidden, SetAccess = private, GetAccess = public)
       
        % Version number for the optim.options.ParticleSwarm2 class. 
        % We do not change the Version property (base class).
        ParticleswarmVersion;
    end
    
% -------------------------------------------------------
% Old properties that will be removed in a future release or are removed
% now (i.e. MaxFunEvals).
    properties (Dependent, Hidden)

%DISPLAYINTERVAL Interval for iterative display
%
        DisplayInterval
        
%FUNVALCHECK Check whether the objective function values are finite.
%
        FunValCheck
        
%INITIALSWARM Initial positions of swarm particles used to seed the
%algorithm
%
        InitialSwarm     
        
%MAXITER Maximum number of iterations allowed
%
        MaxIter

%MINFRACTIONNEIGHBORS Minimum fraction of SwarmSize that will be considered
%as neighbors
%
        MinFractionNeighbors
        
%OUTPUTFCNS Functions that get iterative data and can change options
%
        OutputFcns
        
%PLOTFCNS Plots various measures of progress while the algorithm executes
%
        PlotFcns
        
%SELFADJUSTMENT PARTICLESWARM2 algorithm adjustment amount due to particles
%individual best objective function value.
%
        SelfAdjustment
        
%SOCIALADJUSTMENT PARTICLESWARM2 algorithm adjustment amount due to
%neighbors best objective function values.
%
        SocialAdjustment
        
%STALLITERLIMIT Maximum number of stalled iterations
%
        StallIterLimit
        
%STALLTIMELIMIT Maximum time of stalled iterations
%
        StallTimeLimit
        
%TOLFUN Termination tolerance on function value
%
        TolFun
        
%VECTORIZED Compute the objective function with a vectorized function call
%
        Vectorized                
    end    
% -------------------------------------------------------

    methods (Hidden)
        
        function obj = Particleswarm2(varargin)
%Particleswarm2 Options for PARTICLESWARM2
%
%   The OPTIM.OPTIONS.PARTICLESWARM2 class allows the user to create a set
%   of options for the PARTICLESWARM2 solver. For a list of options that can
%   be set, see the documentation for PARTICLESWARM2.
%
%   OPTS = OPTIM.OPTIONS.PARTICLESWARM2 creates a set of options for
%   PARTICLESWARM2 with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.PARTICLESWARM2(PARAM, VAL, ...) creates a set of
%   options for PARTICLESWARM2 with the named parameters altered with the
%   specified values.
%
%   OPTS = OPTIM.OPTIONS.PARTICLESWARM2(OLDOPTS, PARAM, VAL, ...) creates a
%   copy of OLDOPTS with the named parameters altered with the specified
%   values.
%
%   See also OPTIM.OPTIONS.SingleAlgorithm2, OPTIM.OPTIONS.SolverOptions2
            
            % Call the superclass constructor
            obj = obj@optim.options.SingleAlgorithm2(varargin{:});
            obj.Version = 1;
            obj.ParticleswarmVersion = 2;
        end
        
    end
    
    % Set/get methods
    methods
        function obj = set.CreationFcn(obj, value)
            obj = setProperty(obj, 'CreationFcn', value);
        end

        function obj = set.Display(obj, value)
            obj = setProperty(obj, 'Display', value, ...
                {'none','off','iter','final'});
        end
        
        function obj = set.DisplayInterval(obj, value)
            obj = setProperty(obj, 'DisplayInterval', value);
        end

        function obj = set.FunValCheck(obj, value)
            obj = setProperty(obj, 'FunValCheck', value);
        end

        function obj = set.HybridFcn(obj, value)
            obj = setProperty(obj, 'HybridFcn', value, ...
                {'fminsearch','fminunc','fmincon2','patternsearch2'});
        end

        function obj = set.InertiaRange(obj, value)
            obj = setProperty(obj, 'InertiaRange', value);
        end

        function obj = set.InitialSwarm(obj, value)
            obj = setProperty(obj, 'InitialSwarm', value);
        end

        function obj = set.InitialSwarmMatrix(obj, value)
            obj = setAliasProperty(obj, 'InitialSwarmMatrix', 'InitialSwarm', value);
        end        
        
        function obj = set.InitialSwarmSpan(obj, value)
            obj = setProperty(obj, 'InitialSwarmSpan', value);
        end        
        
        function obj = set.MaxIter(obj, value)
            obj = setProperty(obj, 'MaxIter', value);
        end
        
        function obj = set.MaxIterations(obj, value)
            obj = setAliasProperty(obj, 'MaxIterations', 'MaxIter', value);
        end        
        
        function obj = set.MaxTime(obj, value)
            obj = setProperty(obj, 'MaxTime', value);
        end
        
        function obj = set.MinNeighborsFraction(obj, value)
            obj = setAliasProperty(obj, 'MinNeighborsFraction', 'MinFractionNeighbors', value);
        end
        
        function obj = set.MinFractionNeighbors(obj, value)
            obj = setProperty(obj, 'MinFractionNeighbors', value);
        end
       
        function obj = set.ObjectiveLimit(obj, value)
            obj = setProperty(obj, 'ObjectiveLimit', value);
        end

        function obj = set.OutputFcns(obj, value)
            obj = setProperty(obj, 'OutputFcns', value);
        end
        
        function obj = set.OutputFcn(obj, value)
            obj = setAliasProperty(obj, 'OutputFcn', 'OutputFcns', value);
        end        
        
        function obj = set.PlotFcns(obj, value)
            obj = setProperty(obj, 'PlotFcns', value);
        end        

        function obj = set.PlotFcn(obj, value)
            obj = setAliasProperty(obj, 'PlotFcn', 'PlotFcns', value);
        end           
        
        function obj = set.SelfAdjustment(obj, value)
            obj = setProperty(obj, 'SelfAdjustment', value);
        end

        function obj = set.SelfAdjustmentWeight(obj, value)
            obj = setAliasProperty(obj, 'SelfAdjustmentWeight', 'SelfAdjustment', value);
        end        
        
        function obj = set.SocialAdjustment(obj, value)
            obj = setProperty(obj, 'SocialAdjustment', value);
        end

        function obj = set.SocialAdjustmentWeight(obj, value)
            obj = setAliasProperty(obj, 'SocialAdjustmentWeight', 'SocialAdjustment', value);
        end        
        
        function obj = set.StallIterLimit(obj, value)
            obj = setProperty(obj, 'StallIterLimit', value);
        end
        
        function obj = set.MaxStallIterations(obj, value)
            obj = setAliasProperty(obj, 'MaxStallIterations', 'StallIterLimit', value);
        end        

        function obj = set.StallTimeLimit(obj, value)
            obj = setProperty(obj, 'StallTimeLimit', value);
        end

        function obj = set.MaxStallTime(obj, value)
            obj = setAliasProperty(obj, 'MaxStallTime', 'StallTimeLimit', value);
        end        
        
        function obj = set.SwarmSize(obj, value)
            obj = setProperty(obj, 'SwarmSize', value);
        end
        
        function obj = set.TolFun(obj, value)
            obj = setAliasProperty(obj, 'TolFun', 'TolFunValue', value);
        end

        function obj = set.FunctionTolerance(obj, value)
            obj = setAliasProperty(obj, 'FunctionTolerance', 'TolFunValue', value);
        end        
        
        function obj = set.UseParallel(obj, value)
            obj = setProperty(obj, 'UseParallel', value);
        end

        function obj = set.Vectorized(obj, value)
            obj = setProperty(obj, 'Vectorized', value);
        end

        function obj = set.UseVectorized(obj, value)
            obj = setNewProperty(obj, 'UseVectorized', value);
        end        
        
        %---------------------- Get functions -----------------------------
        
        function value = get.CreationFcn(obj)
            value = obj.OptionsStore.Options.CreationFcn;
        end

        function value = get.Display(obj)
            value = obj.OptionsStore.Options.Display;
        end

        function value = get.DisplayInterval(obj)
            value = obj.OptionsStore.Options.DisplayInterval;
        end

        function value = get.FunValCheck(obj)
            value = obj.OptionsStore.Options.FunValCheck;
        end

        function value = get.HybridFcn(obj)
            value = obj.OptionsStore.Options.HybridFcn;
        end

        function value = get.InertiaRange(obj)
            value = obj.OptionsStore.Options.InertiaRange;
        end

        function value = get.InitialSwarm(obj)
            value = obj.OptionsStore.Options.InitialSwarm;
        end
        
        function value = get.InitialSwarmMatrix(obj)
            value = obj.OptionsStore.Options.InitialSwarm;
        end        

        function value = get.InitialSwarmSpan(obj)
            value = obj.OptionsStore.Options.InitialSwarmSpan;
        end
        
        function value = get.MaxIter(obj)
            value = obj.OptionsStore.Options.MaxIter;
        end

        function value = get.MaxIterations(obj)
            value = obj.OptionsStore.Options.MaxIter;
        end        
        
        function value = get.MaxTime(obj)
            value = obj.OptionsStore.Options.MaxTime;
        end
        
        function value = get.MinFractionNeighbors(obj)
            value = obj.OptionsStore.Options.MinFractionNeighbors;
        end

        function value = get.MinNeighborsFraction(obj)
            value = obj.OptionsStore.Options.MinFractionNeighbors;
        end        
        
        function value = get.ObjectiveLimit(obj)
            value = obj.OptionsStore.Options.ObjectiveLimit;
        end

        function value = get.OutputFcns(obj)
            value = obj.OptionsStore.Options.OutputFcns;
        end
        
        function value = get.OutputFcn(obj)
            value = obj.OptionsStore.Options.OutputFcns;
        end        
        
        function value = get.PlotFcns(obj)
            value = obj.OptionsStore.Options.PlotFcns;
        end     
        
        function value = get.PlotFcn(obj)
            value = obj.OptionsStore.Options.PlotFcns;
        end           

        function value = get.SelfAdjustment(obj)
            value = obj.OptionsStore.Options.SelfAdjustment;
        end

        function value = get.SelfAdjustmentWeight(obj)
            value = obj.OptionsStore.Options.SelfAdjustment;
        end        
        
        function value = get.SocialAdjustment(obj)
            value = obj.OptionsStore.Options.SocialAdjustment;
        end

        function value = get.SocialAdjustmentWeight(obj)
            value = obj.OptionsStore.Options.SocialAdjustment;
        end        
        
        function value = get.StallIterLimit(obj)
            value = obj.OptionsStore.Options.StallIterLimit;
        end
        
        function value = get.MaxStallIterations(obj)
            value = obj.OptionsStore.Options.StallIterLimit;
        end        

        function value = get.StallTimeLimit(obj)
            value = obj.OptionsStore.Options.StallTimeLimit;
        end

        function value = get.MaxStallTime(obj)
            value = obj.OptionsStore.Options.StallTimeLimit;
        end        
        
        function value = get.SwarmSize(obj)
            value = obj.OptionsStore.Options.SwarmSize;
        end        

        function value = get.TolFun(obj)
            value = obj.OptionsStore.Options.TolFunValue;
        end

        function value = get.FunctionTolerance(obj)
            value = obj.OptionsStore.Options.TolFunValue;
        end        
        
        function value = get.UseParallel(obj)
            value = obj.OptionsStore.Options.UseParallel;
        end

        function value = get.Vectorized(obj)
            value = obj.OptionsStore.Options.Vectorized;
        end
        
        function value = get.UseVectorized(obj)
            value = strcmp(obj.OptionsStore.Options.Vectorized, 'on');
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
            
            % Call the superclass method
            OptionsStruct = extractOptionsStructure@optim.options.SingleAlgorithm2(obj);
            
            % Add TolFun to the options struct since particleswarm2 doesn't
            % recognize TolFunValue
            OptionsStruct.TolFun = OptionsStruct.TolFunValue;
        end
        
    end   % Hidden methods
    
    % Load old objects
    methods (Static = true)
        function obj = loadobj(obj)    

            % Upgrade from 14b, 15a, 15b. These objects will come in as
            % structures due, in part, to the rename of
            % "ParticeSwarmVersion" to "ParticleswarmVersion"
            if isstruct(obj) && obj.ParticleSwarmVersion == 1           
                % Save the existing structure
                s = obj;
                
                % Create a new object
                obj = optim.options.Particleswarm2;
                
                % Call the superclass method to upgrade the object
                obj = upgradeFrom13a(obj, s); 
                                
                % Add TolFunValue/ remove TolFun
                obj.OptionsStore.Defaults.TolFunValue = 1e-6;
                obj.OptionsStore.SetByUser.TolFunValue = obj.OptionsStore.SetByUser.TolFun;
                obj.OptionsStore.Options.TolFunValue = obj.OptionsStore.Options.TolFun;
                obj.OptionsStore.SetByUser = rmfield(obj.OptionsStore.SetByUser,'TolFun');
                obj.OptionsStore.Options = rmfield(obj.OptionsStore.Options,'TolFun');
                obj.OptionsStore.Defaults = rmfield(obj.OptionsStore.Defaults,'TolFun');
            end
            
            obj.Version = 1;            
            % Set the version number
            obj.ParticleswarmVersion = 2;            
            
        end
    end        
    
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
%   the full OptionsStore. See below for an example for Particleswarm2.

% Define the option defaults for the solver
OS.Defaults.CreationFcn = @pswcreationuniform2;
OS.Defaults.Display = 'final';
OS.Defaults.DisplayInterval = 1;
OS.Defaults.FunValCheck = 'off';
OS.Defaults.HybridFcn = [];
OS.Defaults.InertiaRange = [0.1 1.1];
OS.Defaults.InitialSwarm = [];
OS.Defaults.InitialSwarmSpan = 2e3;
OS.Defaults.MaxIter = '200*numberofvariables';
OS.Defaults.MaxTime = inf;
OS.Defaults.MinFractionNeighbors = 0.25;
OS.Defaults.ObjectiveLimit = -inf;
OS.Defaults.OutputFcns = [];
OS.Defaults.PlotFcns = [];
OS.Defaults.SelfAdjustment = 1.49;
OS.Defaults.SocialAdjustment = 1.49;
OS.Defaults.StallIterLimit = 20;
OS.Defaults.StallTimeLimit = inf;
OS.Defaults.SwarmSize = 'min(100,10*numberofvariables)';
OS.Defaults.TolFunValue = 1e-6;
OS.Defaults.UseParallel = false;
OS.Defaults.Vectorized = 'off';

% Call the package function to generate the OptionsStore
OS = optim.options.generateSingleAlgorithmOptionsStore2(OS);
end
