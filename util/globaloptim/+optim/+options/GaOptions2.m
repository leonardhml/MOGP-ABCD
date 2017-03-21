classdef (Sealed) GaOptions2 < optim.options.GacommonOptions2
%GaOptions2 Options for GAOPTIONS2
%
%   The OPTIM.OPTIONS.GAOPTIONS2 class allows the user to create a set
%   of options for the GAOPTIONS2 solver. For a list of options that can
%   be set, see the documentation for ga2.
%
%   OPTS = OPTIM.OPTIONS.GAOPTIONS2 creates a set of options for
%   GAOPTIONS2 with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.GAOPTIONS2(PARAM, VAL, ...) creates a set of
%   options for GAOPTIONS2 with the named parameters altered with the
%   specified values.
%
%   OPTS = OPTIM.OPTIONS.GAOPTIONS2(OLDOPTS, PARAM, VAL, ...) creates a
%   copy of OLDOPTS with the named parameters altered with the specified
%   values.
%
%   See also OPTIM.OPTIONS.SingleAlgorithm2, OPTIM.OPTIONS.SolverOptions2

%   Copyright 2015-2015 The MathWorks, Inc.
    
    properties (Dependent)
%ELITECOUNT  Count of individuals in current generation guaranteed to 
% survive to next generation
%
%   For more information, type "doc ga2" and see the "Options"
%   section in the GAOPTIONS2 documentation page.
        EliteCount       
        
%FITNESSLIMIT Limit for fitness function value at which the algorithm stops
%
%   For more information, type "doc ga2" and see the "Options"
%   section in the GAOPTIONS2 documentation page.
        FitnessLimit        
        
%FITNESSSCALINGFCN Function that scales the values of the fitness function
%
%   For more information, type "doc ga2" and see the "Options"
%   section in the GAOPTIONS2 documentation page.
        FitnessScalingFcn          
        
%HYBRIDFCN Function that continues the optimization after ga2 finishes
%
%   For more information, type "doc ga2" and see the "Options"
%   section in the GA2 documentation page.
        HybridFcn          
        
%MAXSTALLTIME Maximum time, in seconds, over which the fitness function 
% value can stall before stopping.
%
%   For more information, type "doc ga2" and see the "Options"
%   section in the GAOPTIONS2 documentation page.
        MaxStallTime           
        
%NONLINEARCONSTRAINTALGORITHM Algorithm used for solving nonlinear
% constrained problems
%
%   For more information, type "doc ga2" and see the "Options"
%   section in the GAOPTIONS2 documentation page.
        NonlinearConstraintAlgorithm        
        
%SELECTIONFCN Function that selects parents of crossover and mutation children
%
%   For more information, type "doc ga2" and see the "Options"
%   section in the GA2 documentation page.
        SelectionFcn           
    end

    properties (Hidden, Dependent)
        
%INITIALPENALTY Initial value of penalty parameter
        InitialPenalty
    
%NONLINCONALGORITHM Nonlinear constraint solver algorithm. 
        NonlinConAlgorithm
    
%PENALTYFACTOR Penalty update parameter
        PenaltyFactor
    
%STALLTEST Stopping test used to determine stall
        StallTest
        
%STALLTIMELIMIT Maximum time, in seconds, over which the fitness function 
% value can stall before stopping.
        StallTimeLimit
        
    end % Hidden properties
    
    properties (Hidden, Access = protected)
        
%OPTIONSSTORE Contains the option values and meta-data for the class
%
        OptionsStore = createOptionsStore;
    end
    
    properties (Hidden)
        
%SOLVERNAME Name of the solver that the options are intended for
%
        SolverName = 'ga2';        
        
    end
    
    properties (SetAccess = private, GetAccess = private)
       
        % Version number for the optim.options.GaOptions2 class. 
        % We do not change the Version property (base class).
        GaVersion = 1
    end
    
% -------------------------------------------------------

    methods (Hidden)
        
        function obj = GaOptions2(varargin)
%GaOptions2 Options for GAOPTIONS2
%
%   The OPTIM.OPTIONS.GAOPTIONS2 class allows the user to create a set
%   of options for the GAOPTIONS2 solver. For a list of options that can
%   be set, see the documentation for GAOPTIONS2.
%
%   OPTS = OPTIM.OPTIONS.GAOPTIONS2 creates a set of options for
%   GAOPTIONS2 with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.GAOPTIONS2(PARAM, VAL, ...) creates a set of
%   options for GAOPTIONS2 with the named parameters altered with the
%   specified values.
%
%   OPTS = OPTIM.OPTIONS.GAOPTIONS2(OLDOPTS, PARAM, VAL, ...) creates a
%   copy of OLDOPTS with the named parameters altered with the specified
%   values.
%
%   See also OPTIM.OPTIONS.SingleAlgorithm2, OPTIM.OPTIONS.SolverOptions2
            
            % Call the superclass constructor
            obj = obj@optim.options.GacommonOptions2(varargin{:});
            
        end
        
    end
    
    % Set/get methods
    methods
       
        function obj = set.EliteCount(obj, value)
            obj = setProperty(obj, 'EliteCount', value);
        end 
        
        function obj = set.FitnessLimit(obj, value)
            obj = setProperty(obj, 'FitnessLimit', value);
        end          
        
        function obj = set.FitnessScalingFcn(obj, value)
            obj = setProperty(obj, 'FitnessScalingFcn', value);
        end 
        
        function obj = set.HybridFcn(obj, value)
            obj = setProperty(obj, 'HybridFcn', value, ...
                 {'fminsearch','fminunc','fmincon2','patternsearch2'});
        end  
        
        function obj = set.InitialPenalty(obj, value)
            obj = setProperty(obj, 'InitialPenalty', value);
        end         
        
        function obj = set.MaxStallTime(obj, value)
            obj = setAliasProperty(obj, 'MaxStallTime', 'StallTimeLimit', value);
        end         
        
        function obj = set.NonlinConAlgorithm(obj, value)
            obj = setProperty(obj, 'NonlinConAlgorithm', value);
        end         
                
        function obj = set.NonlinearConstraintAlgorithm(obj, value)
            obj = setAliasProperty(obj, 'NonlinearConstraintAlgorithm', 'NonlinConAlgorithm', value);
        end        
        
        function obj = set.PenaltyFactor(obj, value)
            obj = setProperty(obj, 'PenaltyFactor', value);
        end                 
        
        function obj = set.SelectionFcn(obj, value)
            obj = setProperty(obj, 'SelectionFcn', value);
        end          
        
        function obj = set.StallTest(obj, value)
            obj = setProperty(obj, 'StallTest', value);
        end                 
        
        function obj = set.StallTimeLimit(obj, value)
            obj = setProperty(obj, 'StallTimeLimit', value);
        end                         
        
        %---------------------- Get functions -----------------------------
        
        function value = get.EliteCount(obj)
           value = obj.OptionsStore.Options.EliteCount;
        end 
        
        function value = get.FitnessLimit(obj)
            value = obj.OptionsStore.Options.FitnessLimit;
        end                  
        
        function value = get.FitnessScalingFcn(obj)
            value = obj.OptionsStore.Options.FitnessScalingFcn;
        end         
        
        function value = get.HybridFcn(obj)
            value = obj.OptionsStore.Options.HybridFcn;
        end           
        
        function value = get.InitialPenalty(obj)
            value = obj.OptionsStore.Options.InitialPenalty;
        end         
        
        function value = get.MaxStallTime(obj)
            value = obj.OptionsStore.Options.StallTimeLimit;
        end         

        function value = get.NonlinConAlgorithm(obj)
            value = obj.OptionsStore.Options.NonlinConAlgorithm;
        end          
        
        function value = get.NonlinearConstraintAlgorithm(obj)
            value = obj.OptionsStore.Options.NonlinConAlgorithm;
        end         
        
        function value = get.PenaltyFactor(obj)
            value = obj.OptionsStore.Options.PenaltyFactor;
        end               
        
        function value = get.SelectionFcn(obj)
            value = obj.OptionsStore.Options.SelectionFcn;
        end           
        
        function value = get.StallTest(obj)
            value = obj.OptionsStore.Options.StallTest;
        end          
        
        function value = get.StallTimeLimit(obj)
            value = obj.OptionsStore.Options.StallTimeLimit;
        end          
    end % get/set methods
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
%   the full OptionsStore. See below for an example for GaOptions2.

% Define the option defaults for the solver
OS.Defaults.CreationFcn = @gacreationuniform2;
OS.Defaults.CrossoverFcn = @crossoverscattered2;
OS.Defaults.CrossoverFraction = 0.8;
OS.Defaults.Display = 'final';
OS.Defaults.EliteCount = '0.05*PopulationSize';
OS.Defaults.FitnessLimit = -Inf;
OS.Defaults.FitnessScalingFcn = @fitscalingrank2;
OS.Defaults.Generations = '100*numberOfVariables';
OS.Defaults.HybridFcn = [];
OS.Defaults.InitialPenalty = 10;
OS.Defaults.InitialPopulation = [];
OS.Defaults.InitialScores = [];
OS.Defaults.MigrationDirection = 'forward';
OS.Defaults.MigrationFraction = 0.2;
OS.Defaults.MigrationInterval = 20;
OS.Defaults.MutationFcn = {@mutationgaussian2, 1, 1};
OS.Defaults.NonlinConAlgorithm = 'auglag';
OS.Defaults.OutputFcns = [];
OS.Defaults.PenaltyFactor = 100;
OS.Defaults.PlotFcns = [];
OS.Defaults.PlotInterval = 1;
OS.Defaults.PopInitRange = [];
OS.Defaults.PopulationSize = '50 when numberOfVariables <= 5, else 200';
OS.Defaults.PopulationType = 'doubleVector';
OS.Defaults.SelectionFcn = @selectionstochunif2;
OS.Defaults.StallGenLimit = 50;
OS.Defaults.StallTest = 'averageChange';
OS.Defaults.StallTimeLimit = Inf;
OS.Defaults.TimeLimit = Inf;
OS.Defaults.TolCon = 1e-3;
OS.Defaults.TolFunValue = 1e-6;
OS.Defaults.UseParallel = false;
OS.Defaults.Vectorized = 'off';

% Call the package function to generate the OptionsStore
OS = optim.options.generateSingleAlgorithmOptionsStore2(OS);
end
