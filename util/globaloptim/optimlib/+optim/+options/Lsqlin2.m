classdef (Sealed) lsqlin2 < optim.options.MultiAlgorithm2
%

%lsqlin2 Options for lsqlin2
%
%   The OPTIM.OPTIONS.lsqlin2 class allows the user to create a set of
%   options for the lsqlin2 solver. For a list of options that can be set,
%   see the documentation for lsqlin2.
%
%   OPTS = OPTIM.OPTIONS.lsqlin2 creates a set of options for lsqlin2
%   with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.lsqlin2(PARAM, VAL, ...) creates a set of options
%   for lsqlin2 with the named parameters altered with the specified
%   values.
%
%   OPTS = OPTIM.OPTIONS.lsqlin2(OLDOPTS, PARAM, VAL, ...) creates a copy
%   of OLDOPTS with the named parameters altered with the specified values.
%
%   See also OPTIM.OPTIONS.MultiAlgorithm2, OPTIM.OPTIONS.SolverOptions2   

%   Copyright 2012-2015 The MathWorks, Inc.    
    
    properties (Dependent)
%CONSTAINTTOLERANCE Tolerance on the constraint violation
% 
%   For more information, type "doc lsqlin2" and see the "Options" section
%   in the lsqlin2 documentation page.        
        ConstraintTolerance    
        
%DISPLAY Level of display
%
%   For more information, type "doc lsqlin2" and see the "Options" section
%   in the lsqlin2 documentation page.        
        Display

%FUNCTIONTOLERANCE Termination tolerance on the change in function value
% 
%   For more information, type "doc lsqlin2" and see the "Options" section
%   in the lsqlin2 documentation page.        
        FunctionTolerance        
        
%JACOBIANMULTIPLYFCN Handle to function that performs Jacobian-vector and
%                    Jacobian-matrix products
%
%   For more information, type "doc lsqlin2" and see the "Options" section
%   in the lsqlin2 documentation page.        
        JacobianMultiplyFcn   
        
%MAXITERATIONS Maximum number of iterations allowed 
%
%   For more information, type "doc lsqlin2" and see the "Options" section
%   in the lsqlin2 documentation page.        
        MaxIterations

%OPTIMALITYTOLERANCE Termination tolerance on the first-order optimality
%                    measure
% 
%   For more information, type "doc lsqlin2" and see the "Options" section
%   in the lsqlin2 documentation page.        
        OptimalityTolerance

%SUBPROBLEMALGORITHM Algorithm used to solve a subproblem
% 
%   For more information, type "doc lsqlin2" and see the "Options" section
%   in the lsqlin2 documentation page.        
        SubproblemAlgorithm   
                
%TYPICALX Typical x values        
% 
%   For more information, type "doc lsqlin2" and see the "Options" section
%   in the lsqlin2 documentation page.        
        TypicalX
    end
    
    % Hidden properties
    properties (Hidden, Dependent)
        
%DIAGNOSTICS Display diagnostic information 
%
        Diagnostics
        
%JACOBMULT Function handle for Jacobian multiply function    
%
        JacobMult
        
%MAXITER Maximum number of iterations allowed 
%
        MaxIter
        
%MAXPCGITER Maximum number of PCG (preconditioned conjugate gradient) 
%           iterations        
% 
        MaxPCGIter
        
%PRECONDBANDWIDTH Upper bandwidth of preconditioner for PCG
% 
        PrecondBandWidth
        
%PRESOLVEOPS Operations to perform during presolve.
%
        PresolveOps
        
%TOLCON Tolerance on the constraint violation
% 
        TolCon        
        
%TOLFUN Termination tolerance on the function value
%       
        TolFun
        
%TOLPCG Termination tolerance on the PCG iteration
% 
        TolPCG        
    end
    
    properties (Hidden, Access = protected)
%OPTIONSSTORE Contains the option values and meta-data for the class
%          
        OptionsStore = createOptionsStore;
    end
    
    properties (Hidden)
%SOLVERNAME Name of the solver that the options are intended for
%          
        SolverName = 'lsqlin2';
    end
    
    
    properties (Hidden, SetAccess = private, GetAccess = public)
        
        % New version property added in third version
        LsqlinVersion
    end    
    
    methods (Hidden)
        
        function obj = lsqlin2(varargin)
%lsqlin2 Options for lsqlin2
%
%   The OPTIM.OPTIONS.lsqlin2 class allows the user to create a set of
%   options for the lsqlin2 solver. For a list of options that can be set,
%   see the documentation for lsqlin2.
%
%   OPTS = OPTIM.OPTIONS.lsqlin2 creates a set of options for lsqlin2
%   with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.lsqlin2(PARAM, VAL, ...) creates a set of options
%   for lsqlin2 with the named parameters altered with the specified
%   values.
%
%   OPTS = OPTIM.OPTIONS.lsqlin2(OLDOPTS, PARAM, VAL, ...) creates a copy
%   of OLDOPTS with the named parameters altered with the specified values.
%
%   See also OPTIM.OPTIONS.MultiAlgorithm2, OPTIM.OPTIONS.SolverOptions2   
            
            % Call the superclass constructor
            obj = obj@optim.options.MultiAlgorithm2(varargin{:});
            
            % Record the class version; Update property 'LsqlinVersion'
            % instead of superclass property 'Version'.
            obj.Version = 2;
            obj.LsqlinVersion = 4;    
        end
        
        function optionFeedback = createOptionFeedback2(obj)
%createOptionFeedback2 Create option feedback string 
%
%   optionFeedback = createOptionFeedback2(obj) creates an option feedback
%   strings that are required by the extended exit messages. OPTIONFEEDBACK
%   is a structure containing strings for the options that appear in the
%   extended exit messages. These strings indicate whether the option is at
%   its 'default' value or has been 'selected'.   

            % It is possible for a user to pass in a vector of options to
            % the solver. Silently use the first element in this array.
            obj = obj(1);
            
            % Check if TolFun is the default value
            if obj.OptionsStore.SetByUser.TolFun
                optionFeedback.TolFun = 'selected';
            else
                optionFeedback.TolFun = 'default';
            end
                     
            % Check if TolFunValue is the default value
            if obj.OptionsStore.SetByUser.TolFunValue
                optionFeedback.TolFunValue = 'selected';
            else
                optionFeedback.TolFunValue = 'default';
            end    
            
            % Check if MaxIter is the default value
            if obj.OptionsStore.SetByUser.MaxIter
                optionFeedback.MaxIter = 'selected';
            else
                optionFeedback.MaxIter = 'default';
            end
                        
            % Check if TolCon is the default value
            if obj.OptionsStore.SetByUser.TolCon
                optionFeedback.TolCon = 'selected';
            else
                optionFeedback.TolCon = 'default';
            end
        end
        
    end
    
    % Set/get methods
    methods
        
        function obj = set.MaxPCGIter(obj, value)
            obj = setProperty(obj, 'MaxPCGIter', value);
        end
        
        function obj = set.PrecondBandWidth(obj, value)
            obj = setProperty(obj, 'PrecondBandWidth', value);
        end
        
        function obj = set.TolPCG(obj, value)
            obj = setProperty(obj, 'TolPCG', value);
        end
        
        function obj = set.MaxIter(obj, value)
            obj = setProperty(obj, 'MaxIter', value);
        end
        
        function obj = set.MaxIterations(obj, value)
            obj = setAliasProperty(obj, 'MaxIterations', 'MaxIter', value);
        end        
        
        function obj = set.Display(obj, value)                                                           
            if strcmpi(value, 'testing')
                % Set Display to the undocumented value, 'testing'.
                obj = setPropertyNoChecks(obj, 'Display', 'testing');
            else
                % Pass the possible values that the Display option can take
                % via the fourth input of setProperty.
                obj = setProperty(obj, 'Display', value, ...
                    {'off','none','final', ...
                    'final-detailed','iter','iter-detailed'});
            end

        end
        
        function obj = set.Diagnostics(obj, value)
            obj = setProperty(obj, 'Diagnostics', value);
        end
        
        function obj = set.JacobMult(obj, value)
            obj = setProperty(obj, 'JacobMult', value);
        end
        
        function obj = set.JacobianMultiplyFcn(obj, value)
            obj = setAliasProperty(obj, 'JacobianMultiplyFcn', 'JacobMult', value);
        end        
        
        function obj = set.PresolveOps(obj, value)
            obj = setProperty(obj, 'PresolveOps', value);
        end        
        
        function obj = set.TolCon(obj, value)
            obj = setProperty(obj, 'TolCon', value);
        end        
        
        function obj = set.ConstraintTolerance(obj, value)
            obj = setAliasProperty(obj, 'ConstraintTolerance', 'TolCon', value);
        end 
        
        function obj = set.TolFun(obj, value)
            obj = setNewProperty(obj, 'TolFun', value);
        end
        
        function obj = set.OptimalityTolerance(obj, value)
            obj = setAliasProperty(obj, 'OptimalityTolerance', 'TolFun', value);
        end        
        
        function obj = set.FunctionTolerance(obj, value)
            obj = setAliasProperty(obj, 'FunctionTolerance', 'TolFunValue', value);
        end  
        
        function obj = set.SubproblemAlgorithm(obj, value)
            obj = setNewProperty(obj, 'SubproblemAlgorithm', value);
        end        
        
        function obj = set.TypicalX(obj, value)
            obj = setProperty(obj, 'TypicalX', value);
        end
        
        %----------------- Get functions ----------------------------------
        function value = get.ConstraintTolerance(obj)
            value = obj.OptionsStore.Options.TolCon;
        end            
        
        function value = get.Diagnostics(obj)
            value = obj.OptionsStore.Options.Diagnostics;
        end
        
        function value = get.Display(obj)
            value = obj.OptionsStore.Options.Display;
        end
        
        function value = get.FunctionTolerance(obj)
            value = obj.OptionsStore.Options.TolFunValue;
        end          
        
        function value = get.JacobMult(obj)
            value = obj.OptionsStore.Options.JacobMult;
        end
        
        function value = get.JacobianMultiplyFcn(obj)
            value = obj.OptionsStore.Options.JacobMult;
        end        
        
        function value = get.MaxIter(obj)
            value = obj.OptionsStore.Options.MaxIter;
        end
        
        function value = get.MaxIterations(obj)
            value = obj.OptionsStore.Options.MaxIter;
        end
        
        function value = get.MaxPCGIter(obj)
            value = obj.OptionsStore.Options.MaxPCGIter;
        end
        
        function value = get.OptimalityTolerance(obj)
            value = obj.OptionsStore.Options.TolFun;
        end 
        
        function value = get.PrecondBandWidth(obj)
            value = obj.OptionsStore.Options.PrecondBandWidth;
        end
        
        function obj = get.PresolveOps(obj)
            obj = obj.OptionsStore.Options.PresolveOps;
        end        
        
        function value = get.SubproblemAlgorithm(obj)
            value = optim.options.OptionAliasStore2.mapOptionFromStore('SubproblemAlgorithm', obj.OptionsStore.Options);
        end
        
        function value = get.TolCon(obj)
            value = obj.OptionsStore.Options.TolCon;
        end        
        
        function value = get.TolFun(obj)
            value = obj.OptionsStore.Options.TolFun;
        end
        
        function value = get.TolPCG(obj)
            value = obj.OptionsStore.Options.TolPCG;
        end
        
        function value = get.TypicalX(obj)
            value = obj.OptionsStore.Options.TypicalX;
        end
        
    end
    
    % Hidden utility methods
    methods (Hidden)
        
        function OptionsStruct = mapOptionsForSolver(~, OptionsStruct)
%mapOptionsForSolver Map structure to an optimset one
%
%   OptionsStruct = mapOptionsForSolver(obj, OptionsStruct) maps the
%   specified structure so it can be used in the solver functions and in
%   OPTIMTOOL.            
            % If TolFun is at a default level, set to empty
            if isfield(OptionsStruct, 'TolFun') && ...
                    strcmp(OptionsStruct.TolFun, 'default dependent on problem')
                OptionsStruct.TolFun = [];
                OptionsStruct.TolFunValue = [];
            end
            
            algIsTRR = ~isfield(OptionsStruct, 'Algorithm') || isempty(OptionsStruct.Algorithm) || ...
                        strcmpi(OptionsStruct.Algorithm,'trust2-region-reflective');
                    
            largeScale = ~isfield(OptionsStruct, 'LargeScale') || isempty(OptionsStruct.LargeScale) || ...
                        strcmpi(OptionsStruct.LargeScale,'on');
                    
            if algIsTRR && ~largeScale
                OptionsStruct.Algorithm = 'active-set';
            end
        end
        
        function [obj, OptimsetStruct] = mapOptimsetToOptions(obj, OptimsetStruct)
%mapOptimsetToOptions Map optimset structure to optimoptions2
%
%   obj = mapOptimsetToOptions(obj, OptimsetStruct) maps specified optimset
%   options, OptimsetStruct, to the equivalent options in the specified
%   optimization object, obj.
%
%   [obj, OptionsStruct] = mapOptimsetToOptions(obj, OptimsetStruct)
%   additionally returns an options structure modified with any conversions
%   that were performed on the options object.
            
            if isfield(OptimsetStruct,'LargeScale') && ~isempty(OptimsetStruct.LargeScale)
                if strcmp(OptimsetStruct.LargeScale, 'on')
                    obj.Algorithm = 'trust2-region-reflective';
                else
                    obj.Algorithm = 'active-set';
                end
            end
            % Also modify the incoming structure.
            if nargout > 1
                OptimsetStruct.Algorithm = obj.Algorithm;
                OptimsetStruct = rmfield(OptimsetStruct,'LargeScale');
            end
        end
        
    end
    
    % Load old objects
    methods (Static = true)
        function obj = loadobj(obj)
            
            % Objects saved in R2013a will come in as structures. 
            if isstruct(obj) && obj.Version == 1

                % Save the existing structure
                s = obj;
                
                % Create a new object
                obj = optim.options.lsqlin2;
                
                % Call the superclass method to upgrade the object
                obj = upgradeFrom13a(obj, s); 
                
                % The SolverVersion property was not present in 13a. We
                % clear it here and the remainer of loadobj will set it
                % correctly.
                obj.LsqlinVersion = [];
                
            end
                        
            % Upgrading to 14b
            if obj.Version < 2
                % Update OptionsStore by taking the loaded OptionsStore and
                % add info for the interior-point algorithm
                 os = createOptionsStore();
                 os.SetByUser = obj.OptionsStore.SetByUser;
                 os.SetByUser.TolCon = false; 
				 os.SetByUser.PresolveOps = false;
                 os.Options = obj.OptionsStore.Options;
                 os.Options.TolCon = os.AlgorithmDefaults{3}.TolCon; 
                 os.Options.PresolveOps = os.AlgorithmDefaults{3}.PresolveOps;
                 os.AlgorithmIndex = [obj.OptionsStore.AlgorithmIndex false];
                 obj.OptionsStore = os;
            end
            
            % Upgrade to 16a
            if isempty(obj.LsqlinVersion) || obj.LsqlinVersion < 4
                % Add TolFunValue
                obj.OptionsStore.AlgorithmDefaults{1}.TolFunValue = 1e-6;
                obj.OptionsStore.IsConstantDefault.TolFunValue = true;
                % Set TolFunValue to whatever of TolFun was saved, but only if the selected algorithm has
                % "FunctionTolerance". Otherwise, set to its default value
                % for another algorithm
                if isfield(obj.OptionsStore.AlgorithmDefaults{obj.OptionsStore.AlgorithmIndex},'TolFunValue') && obj.OptionsStore.SetByUser.TolFun
                    obj.OptionsStore.SetByUser.TolFunValue = obj.OptionsStore.SetByUser.TolFun;
                    obj.OptionsStore.Options.TolFunValue = obj.OptionsStore.Options.TolFun;
                else
                    obj.OptionsStore.SetByUser.TolFunValue = false;
                    obj.OptionsStore.Options.TolFunValue = obj.OptionsStore.AlgorithmDefaults{1}.TolFunValue;
                end
                
                % Objects prior to 15b are missing display-related fields
                % in OptionsStore
                obj.OptionsStore = optim.options.generateMultiAlgorithmDisplayOptions2( ...
                                    obj.OptionsStore,'optim.options.lsqlin2'); 
            end            
            
            % Set the version number
            obj.LsqlinVersion = 4;            
            
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
%   Class authors must create a structure containing the following fields:-
%
%   AlgorithmNames   : Cell array of algorithm names for the solver
%   DefaultAlgorithm : String containing the name of the default algorithm
%   AlgorithmDefaults: Cell array of structures. AlgorithmDefaults{i}
%                      holds a structure containing the defaults for 
%                      AlgorithmNames{i}.
%
%   This structure must then be passed to the
%   optim.options.generateMultiAlgorithmOptionsStore2 function to create
%   the full OptionsStore. See below for an example for lsqlin2.

% Define the algorithm names
OS.AlgorithmNames = {'trust2-region-reflective', 'active-set', 'interior-point'};

% Define the default algorithm
OS.DefaultAlgorithm = 'trust2-region-reflective';

% Define the defaults for each algorithm
% trust2-region-reflective
OS.AlgorithmDefaults{1}.Diagnostics = 'off';
OS.AlgorithmDefaults{1}.Display = 'final';
OS.AlgorithmDefaults{1}.JacobMult = [];
OS.AlgorithmDefaults{1}.MaxIter = 200;
OS.AlgorithmDefaults{1}.MaxPCGIter = 'max(1,floor(numberOfVariables/2))';
OS.AlgorithmDefaults{1}.PrecondBandWidth = 0;
OS.AlgorithmDefaults{1}.TolFun = 100*eps;
OS.AlgorithmDefaults{1}.TolFunValue = 100*eps;
OS.AlgorithmDefaults{1}.TolPCG = 0.1;
OS.AlgorithmDefaults{1}.TypicalX = 'ones(numberOfVariables,1)';

% active-set
OS.AlgorithmDefaults{2}.Diagnostics = 'off';
OS.AlgorithmDefaults{2}.Display = 'final';
OS.AlgorithmDefaults{2}.MaxIter = 200;

% interior-point
OS.AlgorithmDefaults{3}.Diagnostics = 'off';
OS.AlgorithmDefaults{3}.Display = 'final';
OS.AlgorithmDefaults{3}.MaxIter = 200;
OS.AlgorithmDefaults{3}.TolFun = 1e-8;
OS.AlgorithmDefaults{3}.TolCon = 1e-8;
OS.AlgorithmDefaults{3}.PresolveOps = [];

% Call the package function to generate the OptionsStore
OS = optim.options.generateMultiAlgorithmOptionsStore2(OS, 'optim.options.lsqlin2');

end
