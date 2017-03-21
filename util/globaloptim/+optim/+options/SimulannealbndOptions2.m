classdef (Sealed) SimulannealbndOptions2 < optim.options.SingleAlgorithm2
%SimulannealbndOptions2 Options for SIMULANNEALBNDOPTIONS2
%
%   The OPTIM.OPTIONS.SIMULANNEALBNDOPTIONS2 class allows the user to create a set
%   of options for the SIMULANNEALBNDOPTIONS2 solver. For a list of options that can
%   be set, see the documentation for simulannealbnd2.
%
%   OPTS = OPTIM.OPTIONS.SIMULANNEALBNDOPTIONS2 creates a set of options for
%   SIMULANNEALBNDOPTIONS2 with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.SIMULANNEALBNDOPTIONS2(PARAM, VAL, ...) creates a set of
%   options for SIMULANNEALBNDOPTIONS2 with the named parameters altered with the
%   specified values.
%
%   OPTS = OPTIM.OPTIONS.SIMULANNEALBNDOPTIONS2(OLDOPTS, PARAM, VAL, ...) creates a
%   copy of OLDOPTS with the named parameters altered with the specified
%   values.
%
%   See also OPTIM.OPTIONS.SingleAlgorithm2, OPTIM.OPTIONS.SolverOptions2

%   Copyright 2015-2015 The MathWorks, Inc.
    
    properties (Dependent)
%ACCEPTANCEFCN Function used to determine if a new point is accepted
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        AcceptanceFcn          
        
%ANNEALINGFCN Function used to generate new points
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        AnnealingFcn        
 
%DATATYPE Type of decision variable
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        DataType        
        
%DISPLAY Level of display
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        Display
        
%FUNCTIONTOLERANCE Tolerance on the function value
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        FunctionTolerance        
        
%HYBRIDFCN Function run during or at the end of iterations of the solver
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        HybridFcn
        
%INITIALTEMPERATURE Initial value of temperature
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        InitialTemperature
        
%MAXFUNCTIONEVALUATIONS Maximum number of objective function evaluations
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        MaxFunctionEvaluations
        
%MAXITERATIONS Maximum number of iterations
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        MaxIterations
        
%MAXSTALLITERATIONS Maximum number of iterations over which the fitness
%function value is allowed to stall
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        MaxStallIterations        
        
%MAXTIME Total time (in seconds) allowed for optimization
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        MaxTime        

%OBJECTIVELIMIT Minimum objective function value desired
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        ObjectiveLimit        
        
%OUTPUTFCN Functions that get iterative data and can change options
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        OutputFcn
        
%PLOTFCN Plots various measures of progress while the algorithm executes
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        PlotFcn       
                      
%REANNEALINTERVAL Reannealing interval
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        ReannealInterval
        
%TEMPERATUREFCN Function used to update temperature schedule
%
%   For more information, type "doc simulannealbnd2" and see the "Options"
%   section in the SIMULANNEALBNDOPTIONS2 documentation page.
        TemperatureFcn
    end
    
        % Hidden properties
    properties (Hidden, Dependent)   
        %DISPLAYINTERVAL Interval for iterative display
        %
        DisplayInterval
        
        %HYBRIDINTERVAL Interval at which hybrid function is called
        %
        HybridInterval        
        
        %MAXFUNEVALS Maximum number of function evaluations allowed
        %
        MaxFunEvals
        
        %MAXITER Maximum number of iterations allowed
        %
        MaxIter

        %OUTPUTFCNS Functions that get iterative data and can change options
        %
        OutputFcns
        
        %PLOTFCNS Plots various measures of progress while the algorithm
        %executes
        %
        PlotFcns
        
        %PLOTINTERVAL Plot functions are called at every interval
        %
        PlotInterval        
        
        %STALLITERLIMIT Maximum number of iterations over which the fitness
        %               function value is allowed to stall
        StallIterLimit       
        
        %TIMELIMIT Total time (in seconds) allowed for optimization
        %
        TimeLimit                
        
        %TOLFUN Termination tolerance on the function value
        %
        TolFun                
    end
    
    properties (Hidden, Access = protected)
        
%OPTIONSSTORE Contains the option values and meta-data for the class
%
        OptionsStore = createOptionsStore;
    end
    
    properties (Hidden)
        
%SOLVERNAME Name of the solver that the options are intended for
%
        SolverName = 'simulannealbnd2';
    end
    
    properties (SetAccess = private, GetAccess = private)
       
        % Version number for the optim.options.SimulannealbndOptions2 class. 
        % We do not change the Version property (base class).
        SimulannealbndVersion = 1
    end
    
% -------------------------------------------------------

    methods (Hidden)
        
        function obj = SimulannealbndOptions2(varargin)
%SimulannealbndOptions2 Options for SIMULANNEALBNDOPTIONS2
%
%   The OPTIM.OPTIONS.SIMULANNEALBNDOPTIONS2 class allows the user to create a set
%   of options for the SIMULANNEALBNDOPTIONS2 solver. For a list of options that can
%   be set, see the documentation for SIMULANNEALBNDOPTIONS2.
%
%   OPTS = OPTIM.OPTIONS.SIMULANNEALBNDOPTIONS2 creates a set of options for
%   SIMULANNEALBNDOPTIONS2 with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.SIMULANNEALBNDOPTIONS2(PARAM, VAL, ...) creates a set of
%   options for SIMULANNEALBNDOPTIONS2 with the named parameters altered with the
%   specified values.
%
%   OPTS = OPTIM.OPTIONS.SIMULANNEALBNDOPTIONS2(OLDOPTS, PARAM, VAL, ...) creates a
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
        function obj = set.AcceptanceFcn(obj, value)
            obj = setProperty(obj, 'AcceptanceFcn', value);
        end
        
        function obj = set.AnnealingFcn(obj, value)
            obj = setProperty(obj, 'AnnealingFcn', value);
        end          
        
        function obj = set.DataType(obj, value)
            obj = setProperty(obj, 'DataType', value);
        end                  
        
        function obj = set.Display(obj, value)
            obj = setProperty(obj, 'Display', value, ...
                {'none','off','iter','diagnose2','final'});
        end
        
        function obj = set.FunctionTolerance(obj, value)
            obj = setAliasProperty(obj, 'FunctionTolerance', 'TolFunValue', value);
        end      

        function obj = set.HybridFcn(obj, value)
            obj = setProperty(obj, 'HybridFcn', value, ...
                 {'fminsearch','fminunc','fmincon2','patternsearch2'});
        end              
        
        function obj = set.InitialTemperature(obj, value)
            obj = setProperty(obj, 'InitialTemperature', value);
        end          
        
        function obj = set.MaxFunctionEvaluations(obj, value)
            obj = setAliasProperty(obj, 'MaxFunctionEvaluations', 'MaxFunEvals', value);
        end         
        
        function obj = set.MaxIterations(obj, value)
            obj = setAliasProperty(obj, 'MaxIterations', 'MaxIter', value);
        end        
        
        function obj = set.MaxStallIterations(obj, value)
            obj = setAliasProperty(obj, 'MaxStallIterations', 'StallIterLimit', value);
        end           
        
        function obj = set.MaxTime(obj, value)
            obj = setAliasProperty(obj, 'MaxTime', 'TimeLimit', value);
        end
        
        function obj = set.ObjectiveLimit(obj, value)
            obj = setProperty(obj, 'ObjectiveLimit', value);
        end        
        
        function obj = set.OutputFcn(obj, value)
            obj = setAliasProperty(obj, 'OutputFcn', 'OutputFcns', value);
        end        
        
        function obj = set.PlotFcn(obj, value)
            obj = setAliasProperty(obj, 'PlotFcn', 'PlotFcns', value);
        end          
        
        function obj = set.ReannealInterval(obj, value)
            obj = setProperty(obj, 'ReannealInterval', value);
        end                 
        
        function obj = set.TemperatureFcn(obj, value)
            obj = setProperty(obj, 'TemperatureFcn', value);
        end                

        %---------------------- Get functions -----------------------------
        
        function value = get.AcceptanceFcn(obj)
            value = obj.OptionsStore.Options.AcceptanceFcn;
        end

        function value = get.AnnealingFcn(obj)
            value = obj.OptionsStore.Options.AnnealingFcn;
        end        
        
        function value = get.DataType(obj)
            value = obj.OptionsStore.Options.DataType;
        end                
        
        function value = get.Display(obj)
            value = obj.OptionsStore.Options.Display;
        end

        function value = get.FunctionTolerance(obj)
            value = obj.OptionsStore.Options.TolFunValue;
        end        
        
        function value = get.HybridFcn(obj)
            value = obj.OptionsStore.Options.HybridFcn;
        end                
        
        function value = get.InitialTemperature(obj)
            value = obj.OptionsStore.Options.InitialTemperature;
        end           
        
        function value = get.MaxFunctionEvaluations(obj)
            value = obj.OptionsStore.Options.MaxFunEvals;
        end                
        
        function value = get.MaxIterations(obj)
            value = obj.OptionsStore.Options.MaxIter;
        end        
        
        function value = get.MaxStallIterations(obj)
            value = obj.OptionsStore.Options.StallIterLimit;
        end         
        
        function value = get.MaxTime(obj)
            value = obj.OptionsStore.Options.TimeLimit;
        end
        
        function value = get.ObjectiveLimit(obj)
            value = obj.OptionsStore.Options.ObjectiveLimit;
        end

        function value = get.OutputFcn(obj)
            value = obj.OptionsStore.Options.OutputFcns;
        end        
               
        function value = get.PlotFcn(obj)
            value = obj.OptionsStore.Options.PlotFcns;
        end           

        function value = get.ReannealInterval(obj)
            value = obj.OptionsStore.Options.ReannealInterval;
        end

        function value = get.TemperatureFcn(obj)
            value = obj.OptionsStore.Options.TemperatureFcn;
        end          
    end % get/set methods
    
    % Set/get methods for hidden options
    methods
        
        %----------------------- Set functions ----------------------------
        function obj = set.DisplayInterval(obj, value)
            obj = setProperty(obj, 'DisplayInterval', value);
        end

        function obj = set.HybridInterval(obj, value)
            if ischar(value)
                if any(strcmpi(value, {'never', 'end'}))
                    obj = setPropertyNoChecks(obj, 'HybridInterval', lower(value));                    
                else
                    error(message('globaloptim:saoptimset2:checkfield2:NotAPosIntEndOrNever',...
                        'HybridInterval', 'end', 'never'));
                end
            else
                obj = setProperty(obj, 'HybridInterval', value);
            end
        end
        
        function obj = set.MaxFunEvals(obj, value)
            obj = setProperty(obj, 'MaxFunEvals', value);
        end
        
        function obj = set.MaxIter(obj, value)
            obj = setProperty(obj, 'MaxIter', value);
        end

        function obj = set.OutputFcns(obj, value)
            obj = setProperty(obj, 'OutputFcns', value);
        end
        
        function obj = set.PlotFcns(obj, value)
            obj = setProperty(obj, 'PlotFcns', value);
        end

        function obj = set.PlotInterval(obj, value)
            obj = setProperty(obj, 'PlotInterval', value);
        end
        
        function obj = set.StallIterLimit(obj, value)
            obj = setProperty(obj, 'StallIterLimit', value);
        end
        
        function obj = set.TimeLimit(obj, value)
            obj = setProperty(obj, 'TimeLimit', value);
        end

        function obj = set.TolFun(obj, value)
            obj = setProperty(obj, 'TolFunValue', value);
        end
        
        %----------------------- Get functions ----------------------------
                
        function value = get.DisplayInterval(obj)
            value = obj.OptionsStore.Options.DisplayInterval;
        end

        function value = get.HybridInterval(obj)
            value = obj.OptionsStore.Options.HybridInterval;
        end
        
        function value = get.MaxFunEvals(obj)
            value = obj.OptionsStore.Options.MaxFunEvals;
        end
        
        function value = get.MaxIter(obj)
            value = obj.OptionsStore.Options.MaxIter;
        end

        function value = get.OutputFcns(obj)
            value = obj.OptionsStore.Options.OutputFcns;
        end
        
        function value = get.PlotFcns(obj)
            value = obj.OptionsStore.Options.PlotFcns;
        end

        function value = get.PlotInterval(obj)
            value = obj.OptionsStore.Options.PlotInterval;
        end
        
        function value = get.StallIterLimit(obj)
            value = obj.OptionsStore.Options.StallIterLimit;
        end
        
        function value = get.TimeLimit(obj)
            value = obj.OptionsStore.Options.TimeLimit;
        end

        function value = get.TolFun(obj)
            value = obj.OptionsStore.Options.TolFunValue;
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
            
            % Add TolFun to the options struct since simulannealbnd2 doesn't
            % recognize TolFunValue
            OptionsStruct.TolFun = OptionsStruct.TolFunValue;
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
%   the full OptionsStore. See below for an example for SimulannealbndOptions2.

% Define the option defaults for the solver
OS.Defaults.AcceptanceFcn = @acceptancesa2;
OS.Defaults.AnnealingFcn = @annealingfast2;
OS.Defaults.DataType = 'double';
OS.Defaults.Display = 'final';
OS.Defaults.DisplayInterval = 10;
OS.Defaults.HybridFcn = [];
OS.Defaults.HybridInterval = 'end';
OS.Defaults.InitialTemperature = 100.0;
OS.Defaults.MaxFunEvals = '3000*numberOfVariables';
OS.Defaults.MaxIter = Inf;
OS.Defaults.ObjectiveLimit = -Inf;
OS.Defaults.OutputFcns = [];
OS.Defaults.PlotFcns = [];
OS.Defaults.PlotInterval = 1;
OS.Defaults.ReannealInterval = 100;
OS.Defaults.StallIterLimit = '500*numberOfVariables';
OS.Defaults.TemperatureFcn = @temperatureexp2;
OS.Defaults.TimeLimit = Inf;
OS.Defaults.TolFunValue = 1e-6;

% Call the package function to generate the OptionsStore
OS = optim.options.generateSingleAlgorithmOptionsStore2(OS);
end
