classdef (Abstract) GacommonOptions2 < optim.options.SingleAlgorithm2
%

%GacommonOptions2 Common base class for Ga2 or Gamultiobj2
%
%   The OPTIM.OPTIONS.GACOMMONOPTIONS2 is an abstract class containing the common
%   set of options for the GA2 or GAMULTIOBJ2 solver. For a list of options
%   that can be set, see the documentation for ga2 or gamultiobj2.
%
%   You cannot create an instance of this class directly.  
%   You must create an instance of the derived classes, either
%   OPTIM.OPTIONS.GAOPTIONS2 or OPTIM.OPTIONS.GAMULTIOBJOPTIONS2.
%
%   See also OPTIM.OPTIONS.SingleAlgorithm2, OPTIM.OPTIONS.SolverOptions2

%   Copyright 2015-2015 The MathWorks, Inc.
    
    properties (Dependent)
%CONSTRAINTTOLERANCE Tolerance on the constraints
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        ConstraintTolerance        
        
%CREATIONFCN Function that creates the initial population
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        CreationFcn

%CROSSOVERFCN Function used to create crossover children
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        CrossoverFcn                
        
%CROSSOVERFRACTION Fraction of population created by the crossover function 
%for next generation
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        CrossoverFraction         
        
%DISPLAY Level of display
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        Display
        
%FUNCTIONTOLERANCE Tolerance on the function value
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        FunctionTolerance             
        
%INITIALPOPULATIONMATRIX Matrix containing the initial population for the algorithm
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        InitialPopulationMatrix                
        
%INITIALPOPULATIONRANGE Range of the individuals in the initial population
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        InitialPopulationRange                        
        
%INITIALSCORESMATRIX Initial scores used to determine fitness
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        InitialScoresMatrix         
        
%MAXGENERATIONS Maximum number of generations (iterations) allowed
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        MaxGenerations
        
%MAXSTALLGENERATIONS Maximum number of generations over which the fitness function 
% value can stall before stopping.
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        MaxStallGenerations        
        
%MAXTIME Total time (in seconds) allowed for optimization
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        MaxTime        

%MUTATIONFCN Function that produces mutation children
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        MutationFcn
        
%OUTPUTFCN Functions that get iterative data and can change options
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        OutputFcn
        
%PLOTFCN Plots various measures of progress while the algorithm executes
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        PlotFcn       
                      
%POPULATIONSIZE Number of individuals in the population
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        PopulationSize
        
%POPULATIONTYPE Data type of individuals in the population
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        PopulationType         
        
%USEPARALLEL Compute the objective function of particles in parallel
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        UseParallel
        
%USEVECTORIZED Compute the objective function with a vectorized function call
%
%   For more information, type "doc ga2 or gamultiobj2" and see the "Options"
%   section in the GA2 or GAMULTIOBJ2 documentation page.
        UseVectorized
    end % Public properties
        

    properties (Hidden, Dependent)
%GENERATIONS Maximum number of generations (iterations) allowed
        Generations
        
%INITIALPOPULATION Matrix containing the initial population for the algorithm
        InitialPopulation                
        
%INITIALSCORES Initial scores used to determine fitness
        InitialScores
        
%MIGRATIONDIRECTION Direction of migration
        MigrationDirection
        
%MIGRATIONFRACTION Fraction of each subpopulation that migrates to a
% different subpopulation
        MigrationFraction
        
%MIGRATIONINTERVAL Number of generations between migrations
        MigrationInterval      
        
%OUTPUTFCNS Functions that get iterative data and can change options
        OutputFcns
        
%PLOTFCNS Plots various measures of progress while the algorithm executes
        PlotFcns
        
%PLOTINTERVAL Number of generations between calls to the plot functions
        PlotInterval
       
%POPINITRANGE Range of the individuals in the initial population
        PopInitRange     
        
%STALLGENLIMIT Maximum number of generations over which the fitness function 
% value can stall before stopping.
        StallGenLimit        
                       
%TIMELIMIT Total time (in seconds) allowed for optimization
        TimeLimit
 
%TOLCON Tolerance on the constraints
        TolCon
        
%TOLFUN Tolerance on the function value
        TolFun        
        
%VECTORIZED Compute the objective function with a vectorized function call
        Vectorized        
        
    end % Hidden properties
    
    properties (Abstract, Dependent)
%SELECTIONFCN Function that selects parents of crossover and mutation children
%
%   For more information, type "doc ga2" and see the "Options"
%   section in the GAMULTIOBJ2 documentation page.
        SelectionFcn             
    end
    
    properties (SetAccess = private, GetAccess = private)
       
        % Version number for the optim.options.GacommonOptions2 class. 
        % We do not change the Version property (base class).
        GacommonVersion = 1
    end
    
% -------------------------------------------------------

    methods (Hidden)
        
        function obj = GacommonOptions2(varargin)
%GacommonOptions2 Options for GA2 or GAMULTIOBJ2
%
%   The OPTIM.OPTIONS.GACOMMONOPTIONS2 class allows the user to create a set
%   of options for the GA2 or GAMULTIOBJ2 solver. For a list of options that can
%   be set, see the documentation for GA2 or GAMULTIOBJ2.
%
%   OPTS = OPTIM.OPTIONS.GACOMMONOPTIONS2 creates a set of options for
%   GA2 or GAMULTIOBJ2 with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.GACOMMONOPTIONS2(PARAM, VAL, ...) creates a set of
%   options for GA2 or GAMULTIOBJ2 with the named parameters altered with the
%   specified values.
%
%   OPTS = OPTIM.OPTIONS.GACOMMONOPTIONS2(OLDOPTS, PARAM, VAL, ...) creates a
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
        function obj = set.ConstraintTolerance(obj, value)
            obj = setAliasProperty(obj, 'ConstraintTolerance', 'TolCon', value);
        end          
        
        function obj = set.CreationFcn(obj, value)
            obj = setProperty(obj, 'CreationFcn', value);
        end          
        
        function obj = set.CrossoverFcn(obj, value)
            obj = setProperty(obj, 'CrossoverFcn', value);
        end         
        
        function obj = set.CrossoverFraction(obj, value)
            obj = setProperty(obj, 'CrossoverFraction', value);
        end          
        
        function obj = set.Display(obj, value)
            obj = setProperty(obj, 'Display', value, ...
                {'none','off','iter','diagnose2','final'});
        end       
        
        function obj = set.FunctionTolerance(obj, value)
            obj = setAliasProperty(obj, 'FunctionTolerance', 'TolFunValue', value);
        end      

        function obj = set.Generations(obj, value)
            obj = setProperty(obj, 'Generations', value);
        end                     
        
        function obj = set.InitialPopulation(obj, value)
            obj = setProperty(obj, 'InitialPopulation', value);
        end        
        
        function obj = set.InitialPopulationMatrix(obj, value)
            obj = setAliasProperty(obj, 'InitialPopulationMatrix', 'InitialPopulation', value);
        end       
        
        function obj = set.InitialPopulationRange(obj, value)
            obj = setAliasProperty(obj, 'InitialPopulationRange', 'PopInitRange', value);
        end           
        
        function obj = set.InitialScores(obj, value)
            obj = setProperty(obj, 'InitialScores', value);
        end         
        
        function obj = set.InitialScoresMatrix(obj, value)
            obj = setAliasProperty(obj, 'InitialScoresMatrix', 'InitialScores', value);
        end        
        
        function obj = set.MaxGenerations(obj, value)
            obj = setAliasProperty(obj, 'MaxGenerations', 'Generations', value);
        end        
        
        function obj = set.MaxStallGenerations(obj, value)
            obj = setAliasProperty(obj, 'MaxStallGenerations', 'StallGenLimit', value);
        end         
        
        function obj = set.MaxTime(obj, value)
            obj = setAliasProperty(obj, 'MaxTime', 'TimeLimit', value);
        end
        
        function obj = set.MigrationDirection(obj, value)
            obj = setProperty(obj, 'MigrationDirection', value);
        end  
        
        function obj = set.MigrationFraction(obj, value)
            obj = setProperty(obj, 'MigrationFraction', value);
        end  
        
        function obj = set.MigrationInterval(obj, value)
            obj = setProperty(obj, 'MigrationInterval', value);
        end  
        
        function obj = set.MutationFcn(obj, value)
            obj = setProperty(obj, 'MutationFcn', value);
        end         
        
        function obj = set.OutputFcn(obj, value)
            obj = setAliasProperty(obj, 'OutputFcn', 'OutputFcns', value);
        end        
        
        function obj = set.OutputFcns(obj, value)
            obj = setProperty(obj, 'OutputFcns', value);
        end 
        
        function obj = set.PlotFcn(obj, value)
            obj = setAliasProperty(obj, 'PlotFcn', 'PlotFcns', value);
        end          
        
        function obj = set.PlotFcns(obj, value)
            obj = setProperty(obj, 'PlotFcns', value);
        end         
        
        function obj = set.PlotInterval(obj, value)
            obj = setProperty(obj, 'PlotInterval', value);
        end            
        
        function obj = set.PopInitRange(obj, value)
            obj = setProperty(obj, 'PopInitRange', value);
        end         
        
        function obj = set.PopulationSize(obj, value)
            obj = setProperty(obj, 'PopulationSize', value);
        end  
        
        function obj = set.PopulationType(obj, value)
            obj = setProperty(obj, 'PopulationType', value);
        end          
        
        function obj = set.StallGenLimit(obj, value)
            obj = setProperty(obj, 'StallGenLimit', value);
        end                         
        
        function obj = set.TimeLimit(obj, value)
            obj = setProperty(obj, 'TimeLimit', value);
        end
        
        function obj = set.TolCon(obj, value)
            obj = setProperty(obj, 'TolCon', value);
        end
        
        function obj = set.TolFun(obj, value)
            obj = setProperty(obj, 'TolFunValue', value);
        end        
        
        function obj = set.UseParallel(obj, value)
            obj = setProperty(obj, 'UseParallel', value);
        end

        function obj = set.UseVectorized(obj, value)
            obj = setNewProperty(obj, 'UseVectorized', value);
        end        
        
        function obj = set.Vectorized(obj, value)
            obj = setProperty(obj, 'Vectorized', value);
        end                
        
        %---------------------- Get functions -----------------------------
        
        function value = get.ConstraintTolerance(obj)
            value = obj.OptionsStore.Options.TolCon;
        end        
        
        function value = get.CreationFcn(obj)
            value = obj.OptionsStore.Options.CreationFcn;
        end      
        
        function value = get.CrossoverFcn(obj)
            value = obj.OptionsStore.Options.CrossoverFcn;
        end   
        
        function value = get.CrossoverFraction(obj)
            value = obj.OptionsStore.Options.CrossoverFraction;
        end           
        
        function value = get.Display(obj)
            value = obj.OptionsStore.Options.Display;
        end
        
        function value = get.FunctionTolerance(obj)
            value = obj.OptionsStore.Options.TolFunValue;
        end        
        
        function value = get.Generations(obj)
            value = obj.OptionsStore.Options.Generations;
        end                      
        
        function value = get.InitialPopulation(obj)
            value = obj.OptionsStore.Options.InitialPopulation;
        end         
        
        function value = get.InitialPopulationMatrix(obj)
            value = obj.OptionsStore.Options.InitialPopulation;
        end      
        
        function value = get.InitialPopulationRange(obj)
            value = obj.OptionsStore.Options.PopInitRange;
        end          
        
        function value = get.InitialScores(obj)
            value = obj.OptionsStore.Options.InitialScores;
        end                 
        
        function value = get.InitialScoresMatrix(obj)
            value = obj.OptionsStore.Options.InitialScores;
        end           
        
        function value = get.MaxGenerations(obj)
            value = obj.OptionsStore.Options.Generations;
        end     
        
        function value = get.MaxStallGenerations(obj)
            value = obj.OptionsStore.Options.StallGenLimit;
        end        
        
        function value = get.MaxTime(obj)
            value = obj.OptionsStore.Options.TimeLimit;
        end
        
        function value = get.MigrationDirection(obj)
            value = obj.OptionsStore.Options.MigrationDirection;
        end         
        
        function value = get.MigrationFraction(obj)
            value = obj.OptionsStore.Options.MigrationFraction;
        end 
        
        function value = get.MigrationInterval(obj)
            value = obj.OptionsStore.Options.MigrationInterval;
        end                 
        
        function value = get.MutationFcn(obj)
            value = obj.OptionsStore.Options.MutationFcn;
        end

        function value = get.OutputFcn(obj)
            value = obj.OptionsStore.Options.OutputFcns;
        end
        
        function value = get.OutputFcns(obj)
            value = obj.OptionsStore.Options.OutputFcns;
        end        
               
        function value = get.PlotFcn(obj)
            value = obj.OptionsStore.Options.PlotFcns;
        end

        function value = get.PlotFcns(obj)
            value = obj.OptionsStore.Options.PlotFcns;
        end        
        
        function value = get.PlotInterval(obj)
            value = obj.OptionsStore.Options.PlotInterval;
        end                 
        
        function value = get.PopInitRange(obj)
            value = obj.OptionsStore.Options.PopInitRange;
        end         
        
        function value = get.PopulationSize(obj)
            value = obj.OptionsStore.Options.PopulationSize;
        end
        
        function value = get.PopulationType(obj)
            value = obj.OptionsStore.Options.PopulationType;
        end               
        
        function value = get.StallGenLimit(obj)
            value = obj.OptionsStore.Options.StallGenLimit;
        end         

        function value = get.TimeLimit(obj)
            value = obj.OptionsStore.Options.TimeLimit;
        end        
        
        function value = get.TolCon(obj)
            value = obj.OptionsStore.Options.TolCon;
        end        
        
        function value = get.TolFun(obj)
            value = obj.OptionsStore.Options.TolFunValue;
        end         
        
        function value = get.UseParallel(obj)
            value = obj.OptionsStore.Options.UseParallel;
        end

        function value = get.UseVectorized(obj)
            value = optim.options.OptionAliasStore2.convertToLogical( ...
                        obj.OptionsStore.Options.Vectorized, 'on');            
        end   
        
        function value = get.Vectorized(obj)
            value = obj.OptionsStore.Options.Vectorized;
        end        
    end % get/set methods
    
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
            
            % Add TolFun to the options struct since GA2* solvers don't
            % recognize TolFunValue
            OptionsStruct.TolFun = OptionsStruct.TolFunValue;
        end
        
    end    
end

