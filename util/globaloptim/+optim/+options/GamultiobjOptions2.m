classdef (Sealed) GamultiobjOptions2 < optim.options.GacommonOptions2
%GamultiobjOptions2 Options for GAMULTIOBJOPTIONS2
%
%   The OPTIM.OPTIONS.GAMULTIOBJOPTIONS2 class allows the user to create a set
%   of options for the GAMULTIOBJOPTIONS2 solver. For a list of options that can
%   be set, see the documentation for gamultiobj2.
%
%   OPTS = OPTIM.OPTIONS.GAMULTIOBJOPTIONS2 creates a set of options for
%   GAMULTIOBJOPTIONS2 with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.GAMULTIOBJOPTIONS2(PARAM, VAL, ...) creates a set of
%   options for GAMULTIOBJOPTIONS2 with the named parameters altered with the
%   specified values.
%
%   OPTS = OPTIM.OPTIONS.GAMULTIOBJOPTIONS2(OLDOPTS, PARAM, VAL, ...) creates a
%   copy of OLDOPTS with the named parameters altered with the specified
%   values.
%
%   See also OPTIM.OPTIONS.SingleAlgorithm2, OPTIM.OPTIONS.SolverOptions2

%   Copyright 2015-2015 The MathWorks, Inc.
    
    properties (Dependent)
%DISTANCEMEASUREFCN Function that computes distance measure of individuals
%
%   For more information, type "doc gamultiobj2" and see the "Options"
%   section in the GAMULTIOBJOPTIONS2 documentation page.
        DistanceMeasureFcn        
        
%HYBRIDFCN Function that continues the optimization after ga2 finishes
%
%   For more information, type "doc gamultiobj2" and see the "Options"
%   section in the GAMULTIOBJ2 documentation page.
        HybridFcn          
        
%PARETOFRACTION Fraction of individuals kept on first Pareto front
%
%   For more information, type "doc gamultiobj2" and see the "Options"
%   section in the GAMULTIOBJOPTIONS2 documentation page.
        ParetoFraction     
        
%SELECTIONFCN Function that selects parents of crossover and mutation children
%
%   For more information, type "doc ga2" and see the "Options"
%   section in the GAMULTIOBJ2 documentation page.
        SelectionFcn           
    end

    properties (Hidden, Access = protected)
        
%OPTIONSSTORE Contains the option values and meta-data for the class
%
        OptionsStore = createOptionsStore;
    end
    
    properties (Hidden)
        
%SOLVERNAME Name of the solver that the options are intended for
%
        SolverName = 'gamultiobj2';
    end
    
    properties (SetAccess = private, GetAccess = private)
       
        % Version number for the optim.options.GamultiobjOptions2 class. 
        % We do not change the Version property (base class).
        GamultiobjVersion = 1
    end
    
    properties (Hidden, Constant)
       InvalidSelectionFcns = {'selectionremainder2','selectionroulette2', ...
                    'selectionstochunif2','selectionuniform2'}; 
    end
    
% -------------------------------------------------------

    methods (Hidden)
        
        function obj = GamultiobjOptions2(varargin)
%GamultiobjOptions2 Options for GAMULTIOBJOPTIONS2
%
%   The OPTIM.OPTIONS.GAMULTIOBJOPTIONS2 class allows the user to create a set
%   of options for the GAMULTIOBJOPTIONS2 solver. For a list of options that can
%   be set, see the documentation for GAMULTIOBJOPTIONS2.
%
%   OPTS = OPTIM.OPTIONS.GAMULTIOBJOPTIONS2 creates a set of options for
%   GAMULTIOBJOPTIONS2 with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.GAMULTIOBJOPTIONS2(PARAM, VAL, ...) creates a set of
%   options for GAMULTIOBJOPTIONS2 with the named parameters altered with the
%   specified values.
%
%   OPTS = OPTIM.OPTIONS.GAMULTIOBJOPTIONS2(OLDOPTS, PARAM, VAL, ...) creates a
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
        function obj = set.DistanceMeasureFcn(obj, value)
            obj = setProperty(obj, 'DistanceMeasureFcn', value);
        end 
        
        function obj = set.HybridFcn(obj, value)
            obj = setProperty(obj, 'HybridFcn', value, {'fgoalattain'});
        end  
        
        function obj = set.ParetoFraction(obj, value)
            obj = setProperty(obj, 'ParetoFraction', value);
        end 
        
        function obj = set.SelectionFcn(obj, value)   
            % First check the type, then refine checking
            obj = setProperty(obj, 'SelectionFcn', value);
            if iscell(value)
                value = value{1};
            end
            if any(strcmpi(char(value),optim.options.GamultiobjOptions2.InvalidSelectionFcns))
                msg = getString(message('globaloptim:validate2:InvalidMultiObjSelectionFcn'));
                error('optim:options:GamultiobjOptions2:InvalidMultiObjSelectionFcn',msg);
            end
        end         
        
        %---------------------- Get functions -----------------------------
        
        function value = get.DistanceMeasureFcn(obj)
           value = obj.OptionsStore.Options.DistanceMeasureFcn;
        end 
        
        function value = get.HybridFcn(obj)
            value = obj.OptionsStore.Options.HybridFcn;
        end           
        
        function value = get.ParetoFraction(obj)
            value = obj.OptionsStore.Options.ParetoFraction;
        end         
        
        function value = get.SelectionFcn(obj)
            value = obj.OptionsStore.Options.SelectionFcn;
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
%   the full OptionsStore. See below for an example for GamultiobjOptions2.

% Define the option defaults for the solver
OS.Defaults.CreationFcn = @gacreationuniform2;
OS.Defaults.CrossoverFcn = @crossoverintermediate2;
OS.Defaults.CrossoverFraction = 0.8;
OS.Defaults.Display = 'final';
OS.Defaults.DistanceMeasureFcn = {@distancecrowding2, 'phenotype'};
OS.Defaults.Generations = '200*numberOfVariables';
OS.Defaults.HybridFcn = [];
OS.Defaults.InitialPopulation = [];
OS.Defaults.InitialScores = [];
OS.Defaults.MigrationDirection = 'forward';
OS.Defaults.MigrationFraction = 0.2;
OS.Defaults.MigrationInterval = 20;
OS.Defaults.MutationFcn = @mutationadaptfeasible2;
OS.Defaults.OutputFcns = [];
OS.Defaults.ParetoFraction = 0.35;
OS.Defaults.PlotFcns = [];
OS.Defaults.PlotInterval = 1;
OS.Defaults.PopInitRange = [];
OS.Defaults.PopulationSize = '50 when numberOfVariables <= 5, else 200';
OS.Defaults.PopulationType = 'doubleVector';
OS.Defaults.StallGenLimit = 100;
OS.Defaults.SelectionFcn = {@selectiontournament2, 2};
OS.Defaults.TimeLimit = Inf;
OS.Defaults.TolCon = 1e-3;
OS.Defaults.TolFunValue = 1e-4;
OS.Defaults.UseParallel = false;
OS.Defaults.Vectorized = 'off';

% Call the package function to generate the OptionsStore
OS = optim.options.generateSingleAlgorithmOptionsStore2(OS);
end
