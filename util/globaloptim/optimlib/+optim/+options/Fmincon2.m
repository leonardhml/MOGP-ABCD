classdef (Sealed) fmincon2 < optim.options.MultiAlgorithm2
    %
    
    %fmincon2 Options for fmincon2
    %
    %   The OPTIM.OPTIONS.fmincon2 class allows the user to create a set of
    %   options for the fmincon2 solver. For a list of options that can be set,
    %   see the documentation for fmincon2.
    %
    %   OPTS = OPTIM.OPTIONS.fmincon2 creates a set of options for fmincon2
    %   with the options set to their default values.
    %
    %   OPTS = OPTIM.OPTIONS.fmincon2(PARAM, VAL, ...) creates a set of options
    %   for fmincon2 with the named parameters altered with the specified
    %   values.
    %
    %   OPTS = OPTIM.OPTIONS.fmincon2(OLDOPTS, PARAM, VAL, ...) creates a copy
    %   of OLDOPTS with the named parameters altered with the specified values.
    %
    %   See also OPTIM.OPTIONS.MultiAlgorithm2, OPTIM.OPTIONS.SolverOptions2
    
    %   Copyright 2012-2015 The MathWorks, Inc.
    
    properties (Dependent)
        
        %CHECKGRADIENTS Compare user-supplied gradients to finite-differencing
        %               derivatives
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        CheckGradients

        %CONSTRAINTTOLERANCE Tolerance on the constraint violation
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        ConstraintTolerance
                
        %DISPLAY Level of display
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        Display
        
        %FINITEDIFFERENCESTEPSIZE Scalar or vector step size factor
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        FiniteDifferenceStepSize
        
        %FINITEDIFFERENCETYPE Finite difference type
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        FiniteDifferenceType

        %FUNCTIONTOLERANCE Termination tolerance on the change in function value
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        FunctionTolerance        

        %HESSIANAPPROXIMATION Specify the method of approximating the Hessian
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        HessianApproximation        
        
        %HESSIANFCN Function handle to a function that computes the Hessian
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        HessianFcn
        
        %HESSIANMULTIPLYFCN Function handle for Hessian multiply function
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        HessianMultiplyFcn
        
        %HONORBOUNDS Determine whether bounds are satisfied at every
        %            iteration
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        HonorBounds        
        
        %MAXFUNCTIONEVALUATIONS Maximum number of function evaluations allowed
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        MaxFunctionEvaluations        
        
        %MAXITERATIONS Maximum number of iterations allowed
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        MaxIterations

        %OBJECTIVELIMIT Lower limit on the objective function
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        ObjectiveLimit
                
        %OPTIMALITYTOLERANCE Termination tolerance on the first-order optimality
        %                    measure
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        OptimalityTolerance               
        
        %OUTPUTFCN Callback that are called at each iteration
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        OutputFcn
        
        %PLOTFCN Plots various measures of progress while the algorithm executes
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        PlotFcn

        %SCALEPROBLEM Determine whether all constraints and the objective function
        %             are normalized
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        ScaleProblem
        
        %SPECIFYCONSTRAINTGRADIENT Gradient for nonlinear constraint functions 
        %                          defined by the caller
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        SpecifyConstraintGradient
        
        %SPECIFYOBJECTIVEGRADIENT Gradient for the objective function
        %                         defined by the caller
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        SpecifyObjectiveGradient
                
        
        %STEPTOLERANCE Termination tolerance on the displacement in x
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        StepTolerance        
        
        %SUBPROBLEMALGORITHM Determines how the iteration step is calculated
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        SubproblemAlgorithm
        
        %TYPICALX Typical x values
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        TypicalX
        
        %USEPARALLEL Estimate gradients in parallel
        %
        %   For more information, type "doc fmincon2" and see the "Options" section
        %   in the fmincon2 documentation page.
        UseParallel
    end
    
    % Hidden properties
    properties (Hidden, Dependent)
        %ALWAYSHONORCONSTRAINTS Determine whether bounds are satisfied at every
        %                       iteration
        %
        AlwaysHonorConstraints
        
        %DERIVATIVECHECK Compare user-supplied derivatives to finite-differencing
        %                derivatives
        %
        DerivativeCheck
        
        %DIAGNOSTICS Display diagnostic information
        %
        Diagnostics
        
        %DIFFMAXCHANGE Maximum change in variables for finite-difference gradients
        %
        DiffMaxChange
        
        %DIFFMINCHANGE Minimum change in variables for finite-difference gradients
        %
        DiffMinChange
        
        %FINDIFFRELSTEP Scalar or vector step size factor
        %
        FinDiffRelStep
        
        %FINDIFFTYPE Finite difference type
        %
        FinDiffType
        
        %FUNVALCHECK Check whether objective function and constraints values are
        %            valid
        %
        FunValCheck
        
        %GRADCONSTR Gradient for nonlinear constraint functions defined by the user
        %
        GradConstr
        
        %GRADOBJ Gradient for the objective function defined by the user
        %
        GradObj
        
        %HESSFCN Function handle to a user-supplied Hessian
        %
        HessFcn
        
        %HESSIAN Specify whether a user-supplied Hessian will be supplied
        %
        Hessian
        
        %HESSMULT Function handle for Hessian multiply function
        %
        HessMult
        
        %HESSPATTERN Sparsity pattern of the Hessian for finite differencing
        %
        HessPattern
        
        %INITBARRIERPARAM Initial barrier value
        %
        InitBarrierParam
        
        %INITTRUSTREGIONRADIUS Initial radius of the trust2 region
        %
        InitTrustRegionRadius
        
        %MAXITER Maximum number of iterations allowed
        %
        MaxIter
        
        %MAXFUNEVALS Maximum number of function evaluations allowed
        %
        MaxFunEvals
        
        %MAXPCGITER Maximum number of PCG (preconditioned conjugate gradient)
        %           iterations
        %
        MaxPCGIter
        
        %MAXPROJCGITER A tolerance for the number of projected conjugate gradient
        %              iterations
        %
        MaxProjCGIter
        
        %MAXSQPITER Maximum number of SQP iterations allowed
        %
        MaxSQPIter

        %PLOTFCNS Plots various measures of progress while the algorithm executes
        %
        PlotFcns
        
        %PRECONDBANDWIDTH Upper bandwidth of preconditioner for PCG
        %
        PrecondBandWidth
        
        %RELLINESRCHBND Relative bound on the line search step length
        %
        RelLineSrchBnd
        
        %RELLINESRCHBNDDURATION Number of iterations for which the bound specified
        %                       in RelLineSrchBnd should be active
        %
        RelLineSrchBndDuration

        %TOLCON Tolerance on the constraint violation
        %
        TolCon
        
        %TOLCONSQP Termination tolerance on inner iteration SQP constraint violation
        %
        TolConSQP
        
        %TOLFUN Termination tolerance on the function value
        %
        TolFun
        
        %TOLPCG Termination tolerance on the PCG iteration
        %
        TolPCG
        
        %TOLPCG Termination tolerance on the PCG iteration
        %
        TolProjCG
        
        %TOLPROJCGABS Absolute tolerance for projected conjugate gradient algorithm
        %
        TolProjCGAbs
        
        %TOLX Termination tolerance on x
        %
        TolX
        
        %NOSTOPIFFLATINFEAS If objective appears flat, only stop if feasible
        %
        NoStopIfFlatInfeas
        
        %PHASEONETOTALSCALING Scale the slack variable in phase 1 of qpsub2
        %
        PhaseOneTotalScaling
        
        %TOLGRADCON Tolerance for detecting stationary point of constraint
        %violation in interior-point algorithm.
        %
        TolGradCon
    end
    
    properties (Hidden, Access = protected)
        %OPTIONSSTORE Contains the option values and meta-data for the class
        %
        OptionsStore = createOptionsStore;
    end
    
    properties (Hidden)
        %SOLVERNAME Name of the solver that the options are intended for
        %
        SolverName = 'fmincon2';
    end
    
    properties (Hidden, SetAccess = private, GetAccess = public)
        
        % New version property added in second version
        FminconVersion
    end
    
    methods (Hidden)
        
        function obj = fmincon2(varargin)
            %fmincon2 Options for fmincon2
            %
            %   The OPTIM.OPTIONS.fmincon2 class allows the user to create a set of
            %   options for the fmincon2 solver. For a list of options that can be set,
            %   see the documentation for fmincon2.
            %
            %   OPTS = OPTIM.OPTIONS.fmincon2 creates a set of options for fmincon2
            %   with the options set to their default values.
            %
            %   OPTS = OPTIM.OPTIONS.fmincon2(PARAM, VAL, ...) creates a set of options
            %   for fmincon2 with the named parameters altered with the specified
            %   values.
            %
            %   OPTS = OPTIM.OPTIONS.fmincon2(OLDOPTS, PARAM, VAL, ...) creates a copy
            %   of OLDOPTS with the named parameters altered with the specified values.
            %
            %   See also OPTIM.OPTIONS.MultiAlgorithm2, OPTIM.OPTIONS.SolverOptions2
            
            % Call the superclass constructor
            obj = obj@optim.options.MultiAlgorithm2(varargin{:});
            
            % Record the class version; Update property 'FminconVersion'
            % instead of superclass property 'Version'.
            obj.Version = 2;
            obj.FminconVersion = 4;
            
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
            
            % Check if TolX is the default value
            if obj.OptionsStore.SetByUser.TolX
                optionFeedback.TolX = 'selected';
            else
                optionFeedback.TolX = 'default';
            end
            
            % Check if MaxIter is the default value
            if obj.OptionsStore.SetByUser.MaxIter
                optionFeedback.MaxIter = 'selected';
            else
                optionFeedback.MaxIter = 'default';
            end
            
            % Check if MaxFunEvals is the default value
            if obj.OptionsStore.SetByUser.MaxFunEvals
                optionFeedback.MaxFunEvals = 'selected';
            else
                optionFeedback.MaxFunEvals = 'default';
            end
            
            % Check if TolCon is the default value
            if obj.OptionsStore.SetByUser.TolCon
                optionFeedback.TolCon = 'selected';
            else
                optionFeedback.TolCon = 'default';
            end
            
            % Check if MaxSQPIter is the default value
            if obj.OptionsStore.SetByUser.MaxSQPIter
                optionFeedback.MaxSQPIter = 'selected';
            else
                optionFeedback.MaxSQPIter = 'default';
            end
            
            % Check if ObjectiveLimit is the default value
            if obj.OptionsStore.SetByUser.ObjectiveLimit
                optionFeedback.ObjectiveLimit = 'selected';
            else
                optionFeedback.ObjectiveLimit = 'default';
            end
        end
        
        
        function obj = replaceSpecialStrings(obj)
            %replaceSpecialStrings Replace special string values
            %
            %   obj = replaceSpecialStrings(obj) replaces special string
            %   option values with their equivalent numerical value. We
            %   currently only use this method to convert FinDiffRelStep.
            %   However, in the future we would like to move the special
            %   string replacement code from the solver files to the
            %   options classes.
            
            % Call a package function to replace string values in
            % FinDiffRelStep.
            obj = optim.options.replaceFinDiffRelStepString2(obj);
            
        end
    end
    
    % Set/get methods
    methods
        
        function obj = set.PhaseOneTotalScaling(obj, value)
            obj = setProperty(obj, 'PhaseOneTotalScaling', value);
        end
        
        function obj = set.NoStopIfFlatInfeas(obj, value)
            obj = setProperty(obj, 'NoStopIfFlatInfeas', value);
        end
        
        function obj = set.TolGradCon(obj, value)
            obj = setProperty(obj, 'TolGradCon', value);
        end
        
        function obj = set.FinDiffRelStep(obj, value)
            obj = setProperty(obj, 'FinDiffRelStep', value);
        end
        
        function obj = set.FiniteDifferenceStepSize(obj, value)
            obj = setAliasProperty(obj, 'FiniteDifferenceStepSize', 'FinDiffRelStep', value);
        end        
        
        function obj = set.Hessian(obj, value)
            obj = setProperty(obj, 'Hessian', value);
            
            % If HessianMultiplyFcn has not been set and Hessian has been
            % set to 'on' or 'user-supplied' and the algorithm is
            % 'trust2-region-reflective', then the user wants to provide a
            % Hessian via the objective function.
            if any(strcmp(value, {'user-supplied', 'on'})) && ...
                    ~obj.OptionsStore.SetByUser.HessMult && ...
                    strcmp(obj.OptionsStore.Options.Algorithm, 'trust2-region-reflective')
                obj = setNewProperty(obj, 'HessianFcn', 'objective');
            end
        end
        
        function obj = set.HessianApproximation(obj, value)
            obj = setNewProperty(obj, 'HessianApproximation', value);
        end
        
        function obj = set.HessMult(obj, value)
            obj = setProperty(obj, 'HessMult', value);
        end
        
        function obj = set.HessianMultiplyFcn(obj, value)
            obj = setNewProperty(obj, 'HessianMultiplyFcn', value);
        end        
        
        function obj = set.HessPattern(obj, value)
            obj = setProperty(obj, 'HessPattern', value);
        end
        
        function obj = set.MaxPCGIter(obj, value)
            obj = setProperty(obj, 'MaxPCGIter', value);
        end
        
        function obj = set.PrecondBandWidth(obj, value)
            obj = setProperty(obj, 'PrecondBandWidth', value);
        end
        
        function obj = set.TolPCG(obj, value)
            obj = setProperty(obj, 'TolPCG', value);
        end
        
        function obj = set.AlwaysHonorConstraints(obj, value)
            obj = setProperty(obj, 'AlwaysHonorConstraints', value);
        end
        
        function obj = set.HonorBounds(obj, value)
            obj = setNewProperty(obj, 'HonorBounds', value);
        end         
        
        function obj = set.HessFcn(obj, value)
            obj = setProperty(obj, 'HessFcn', value);
        end
        
        function obj = set.HessianFcn(obj, value)
            obj = setAliasProperty(obj, 'HessianFcn', 'HessFcn', value);
        end        
        
        function obj = set.InitBarrierParam(obj, value)
            obj = setProperty(obj, 'InitBarrierParam', value);
        end
        
        function obj = set.InitTrustRegionRadius(obj, value)
            obj = setProperty(obj, 'InitTrustRegionRadius', value);
        end
        
        function obj = set.MaxProjCGIter(obj, value)
            obj = setProperty(obj, 'MaxProjCGIter', value);
        end
        
        function obj = set.ObjectiveLimit(obj, value)
            obj = setProperty(obj, 'ObjectiveLimit', value);
        end
        
        function obj = set.ScaleProblem(obj, value)
            obj = setNewProperty(obj, 'ScaleProblem', value, {'none'; 'obj-and-constr'});
        end
        
        function obj = set.SubproblemAlgorithm(obj, value)
            obj = setNewProperty(obj, 'SubproblemAlgorithm', value);
        end
        
        function obj = set.TolProjCG(obj, value)
            obj = setProperty(obj, 'TolProjCG', value);
        end
        
        function obj = set.TolProjCGAbs(obj, value)
            obj = setProperty(obj, 'TolProjCGAbs', value);
        end
        
        function obj = set.MaxIter(obj, value)
            obj = setProperty(obj, 'MaxIter', value);
        end
        
        function obj = set.MaxIterations(obj, value)
            obj = setAliasProperty(obj, 'MaxIterations', 'MaxIter', value);        
        end        

        function obj = set.MaxFunEvals(obj, value)
            obj = setProperty(obj, 'MaxFunEvals', value);
        end
        
        function obj = set.MaxFunctionEvaluations(obj, value)
            obj = setAliasProperty(obj, 'MaxFunctionEvaluations', 'MaxFunEvals', value);        
        end                
        
        function obj = set.TolX(obj, value)
            obj = setProperty(obj, 'TolX', value);
        end
        
        function obj = set.StepTolerance(obj, value)
            obj = setAliasProperty(obj, 'StepTolerance', 'TolX', value);        
        end                
        
        function obj = set.MaxSQPIter(obj, value)
            obj = setProperty(obj, 'MaxSQPIter', value);
        end
        
        function obj = set.RelLineSrchBnd(obj, value)
            obj = setProperty(obj, 'RelLineSrchBnd', value);
        end
        
        function obj = set.RelLineSrchBndDuration(obj, value)
            obj = setProperty(obj, 'RelLineSrchBndDuration', value);
        end
        
        function obj = set.TolConSQP(obj, value)
            obj = setProperty(obj, 'TolConSQP', value);
        end
        
        function obj = set.Display(obj, value)
            if strcmpi(value, 'testing')
                % Set Display to the undocumented value, 'testing'.
                obj = setPropertyNoChecks(obj, 'Display', 'testing');
            else
                % Pass the possible values that the Display option can take via
                % the fourth input of setProperty.
                obj = setProperty(obj, 'Display', value, ...
                    {'off','none','notify','notify-detailed','final', ...
                    'final-detailed','iter','iter-detailed'});
            end
        end
        
        function obj = set.DerivativeCheck(obj, value)
            obj = setProperty(obj, 'DerivativeCheck', value);
        end
        
        function obj = set.CheckGradients(obj, value)
            obj = setNewProperty(obj, 'CheckGradients', value);
        end                
        
        function obj = set.Diagnostics(obj, value)
            obj = setProperty(obj, 'Diagnostics', value);
        end
        
        function obj = set.DiffMinChange(obj, value)
            obj = setProperty(obj, 'DiffMinChange', value);
        end
        
        function obj = set.DiffMaxChange(obj, value)
            obj = setProperty(obj, 'DiffMaxChange', value);
        end
        
        function obj = set.FinDiffType(obj, value)
            obj = setProperty(obj, 'FinDiffType', value);
            % If we get here, the property set has been successful and we
            % can update the OptionsStore
            if ~obj.OptionsStore.SetByUser.FinDiffRelStep
                obj.OptionsStore.Options.FinDiffRelStep = ...
                    optim.options.getDefaultFinDiffRelStep2(...
                    obj.OptionsStore.Options.FinDiffType);
            end
        end
        
        function obj = set.FiniteDifferenceType(obj, value)
            obj = setAliasProperty(obj, 'FiniteDifferenceType', 'FinDiffType', value);     
            % If we get here, the property set has been successful and we
            % can update the OptionsStore
            if ~obj.OptionsStore.SetByUser.FinDiffRelStep
                obj.OptionsStore.Options.FinDiffRelStep = ...
                    optim.options.getDefaultFinDiffRelStep2(...
                    obj.OptionsStore.Options.FinDiffType);
            end            
        end         
        
        function obj = set.FunValCheck(obj, value)
            obj = setProperty(obj, 'FunValCheck', value);
        end
        
        function obj = set.GradConstr(obj, value)
            obj = setProperty(obj, 'GradConstr', value);
        end
        
        function obj = set.SpecifyConstraintGradient(obj, value)
            obj = setNewProperty(obj, 'SpecifyConstraintGradient', value);
        end        
        
        function obj = set.GradObj(obj, value)
            obj = setProperty(obj, 'GradObj', value);
        end
        
        function obj = set.SpecifyObjectiveGradient(obj, value)
            obj = setNewProperty(obj, 'SpecifyObjectiveGradient', value);
        end        
        
        function obj = set.OutputFcn(obj, value)
            obj = setProperty(obj, 'OutputFcn', value);
        end
        
        function obj = set.PlotFcns(obj, value)
            obj = setProperty(obj, 'PlotFcns', value);
        end
        
        function obj = set.PlotFcn(obj, value)
            obj = setAliasProperty(obj, 'PlotFcn', 'PlotFcns', value);        
        end          
        
        function obj = set.TolFun(obj, value)
            obj = setNewProperty(obj, 'TolFun', value);
        end
        
        function obj = set.FunctionTolerance(obj, value)
            obj = setAliasProperty(obj, 'FunctionTolerance', 'TolFunValue', value);        
        end  
        
        function obj = set.OptimalityTolerance(obj, value)
            obj = setAliasProperty(obj, 'OptimalityTolerance', 'TolFun', value);        
        end          
        
        function obj = set.TolCon(obj, value)
            obj = setProperty(obj, 'TolCon', value);
        end
        
        function obj = set.ConstraintTolerance(obj, value)
            obj = setAliasProperty(obj, 'ConstraintTolerance', 'TolCon', value);        
        end         
        
        function obj = set.TypicalX(obj, value)
            obj = setProperty(obj, 'TypicalX', value);
        end
        
        function obj = set.UseParallel(obj, value)
            obj = setProperty(obj, 'UseParallel', value);
        end
               
        %--------------- Get functions -----------------------------
        
        function value = get.AlwaysHonorConstraints(obj)
            value = obj.OptionsStore.Options.AlwaysHonorConstraints;
        end

        function value = get.CheckGradients(obj)
            value = optim.options.OptionAliasStore2.convertToLogical( ...
                        obj.OptionsStore.Options.DerivativeCheck, 'on');
        end        

        function value = get.ConstraintTolerance(obj)
            value = obj.OptionsStore.Options.TolCon;
        end                
        
        function value = get.DerivativeCheck(obj)
            value = obj.OptionsStore.Options.DerivativeCheck;            
        end
        
        function value = get.Diagnostics(obj)
            value = obj.OptionsStore.Options.Diagnostics;
        end
        
        function value = get.DiffMaxChange(obj)
            value = obj.OptionsStore.Options.DiffMaxChange;
        end
        
        function value = get.DiffMinChange(obj)
            value = obj.OptionsStore.Options.DiffMinChange;
        end
        
        function value = get.Display(obj)
            value = obj.OptionsStore.Options.Display;
        end
        
        function value = get.FinDiffRelStep(obj)
            value = obj.OptionsStore.Options.FinDiffRelStep;
        end
        
        function value = get.FiniteDifferenceStepSize(obj)
            value = obj.OptionsStore.Options.FinDiffRelStep;
        end
        
        function value = get.FinDiffType(obj)
            value = obj.OptionsStore.Options.FinDiffType;
        end
        
        function value = get.FiniteDifferenceType(obj)
            value = obj.OptionsStore.Options.FinDiffType;
        end        
        
        function value = get.FunctionTolerance(obj)
            value = obj.OptionsStore.Options.TolFunValue;
        end        
        
        function value = get.FunValCheck(obj)
            value = obj.OptionsStore.Options.FunValCheck;
        end
        
        function value = get.GradConstr(obj)
            value = obj.OptionsStore.Options.GradConstr;             
        end
        
        function value = get.GradObj(obj)
            value = obj.OptionsStore.Options.GradObj;               
        end
        
        function value = get.Hessian(obj)
            value = obj.OptionsStore.Options.Hessian;
        end
        
        function value = get.HessianApproximation(obj)
            value = optim.options.OptionAliasStore2.mapOptionFromStore('HessianApproximation', obj.OptionsStore.Options);
        end
        
        function value = get.HessFcn(obj)
            value = obj.OptionsStore.Options.HessFcn;
        end
        
        function value = get.HessianFcn(obj)
            value = obj.OptionsStore.Options.HessFcn;
        end        
        
        function value = get.HessMult(obj)
            value = obj.OptionsStore.Options.HessMult;
        end
        
        function value = get.HessianMultiplyFcn(obj)
            value = obj.OptionsStore.Options.HessMult;
        end        
        
        function value = get.HessPattern(obj)
            value = obj.OptionsStore.Options.HessPattern;
        end
        
        function value = get.HonorBounds(obj)
            value = optim.options.OptionAliasStore2.convertToLogical( ...
                obj.OptionsStore.Options.AlwaysHonorConstraints, 'bounds');
        end        
        
        function value = get.InitBarrierParam(obj)
            value = obj.OptionsStore.Options.InitBarrierParam;
        end
        
        function value = get.InitTrustRegionRadius(obj)
            value = obj.OptionsStore.Options.InitTrustRegionRadius;
        end
        
        function value = get.MaxIter(obj)
            value = obj.OptionsStore.Options.MaxIter;
        end
        
        function value = get.MaxIterations(obj)
            value = obj.OptionsStore.Options.MaxIter;
        end
        
        function value = get.MaxFunEvals(obj)
            value = obj.OptionsStore.Options.MaxFunEvals;
        end
        
        function value = get.MaxFunctionEvaluations(obj)
            value = obj.OptionsStore.Options.MaxFunEvals;
        end        
        
        function value = get.MaxPCGIter(obj)
            value = obj.OptionsStore.Options.MaxPCGIter;
        end
        
        function value = get.MaxProjCGIter(obj)
            value = obj.OptionsStore.Options.MaxProjCGIter;
        end
        
        function value = get.MaxSQPIter(obj)
            value = obj.OptionsStore.Options.MaxSQPIter;
        end
        
        function value = get.ObjectiveLimit(obj)
            value = obj.OptionsStore.Options.ObjectiveLimit;
        end
        
        function value = get.OptimalityTolerance(obj)
            value = obj.OptionsStore.Options.TolFun;
        end        
        
        function value = get.OutputFcn(obj)
            value = obj.OptionsStore.Options.OutputFcn;
        end
        
        function value = get.PlotFcns(obj)
            value = obj.OptionsStore.Options.PlotFcns;
        end
        
        function value = get.PlotFcn(obj)
            value = obj.OptionsStore.Options.PlotFcns;
        end        
        
        function value = get.PrecondBandWidth(obj)
            value = obj.OptionsStore.Options.PrecondBandWidth;
        end
        
        function value = get.RelLineSrchBnd(obj)
            value = obj.OptionsStore.Options.RelLineSrchBnd;
        end
        
        function value = get.RelLineSrchBndDuration(obj)
            value = obj.OptionsStore.Options.RelLineSrchBndDuration;
        end
        
        function value = get.ScaleProblem(obj)
            value = optim.options.OptionAliasStore2.convertToLogical( ...
                obj.OptionsStore.Options.ScaleProblem, 'obj-and-constr');
        end
        
        function value = get.SpecifyConstraintGradient(obj)
            value = optim.options.OptionAliasStore2.convertToLogical( ...
                obj.OptionsStore.Options.GradConstr, 'on');
        end
        
        function value = get.SpecifyObjectiveGradient(obj)
            value = optim.options.OptionAliasStore2.convertToLogical( ...
                        obj.OptionsStore.Options.GradObj, 'on');
        end       
        
        function value = get.StepTolerance(obj)
            value = obj.OptionsStore.Options.TolX;
        end           
        
        function value = get.SubproblemAlgorithm(obj)
            value = optim.options.OptionAliasStore2.mapOptionFromStore('SubproblemAlgorithm', obj.OptionsStore.Options);
        end
        
        function value = get.TolCon(obj)
            value = obj.OptionsStore.Options.TolCon;
        end
        
        function value = get.TolConSQP(obj)
            value = obj.OptionsStore.Options.TolConSQP;
        end
        
        function value = get.TolFun(obj)
            value = obj.OptionsStore.Options.TolFun;
        end
        
        function value = get.TolGradCon(obj)
            value = obj.OptionsStore.Options.TolGradCon;
        end
        
        function value = get.TolPCG(obj)
            value = obj.OptionsStore.Options.TolPCG;
        end
        
        function value = get.TolProjCG(obj)
            value = obj.OptionsStore.Options.TolProjCG;
        end
        
        function value = get.TolProjCGAbs(obj)
            value = obj.OptionsStore.Options.TolProjCGAbs;
        end
        
        function value = get.TolX(obj)
            value = obj.OptionsStore.Options.TolX;
        end
        
        function value = get.TypicalX(obj)
            value = obj.OptionsStore.Options.TypicalX;
        end
        
        function value = get.UseParallel(obj)
            value = obj.OptionsStore.Options.UseParallel;
        end
        
        function value = get.NoStopIfFlatInfeas(obj)
            value = obj.OptionsStore.Options.NoStopIfFlatInfeas;
        end
        
        function value = get.PhaseOneTotalScaling(obj)
            value = obj.OptionsStore.Options.PhaseOneTotalScaling;
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
            OptionsStruct = extractOptionsStructure@optim.options.MultiAlgorithm2(obj);
            
            % If HessFcn or HessMult are set to non-empty, this means that the user
            % will supply the Hessian. We need to set the Hessian option to 'user-supplied' 
			% or 'on' (they are equivalent).
            if ~isempty(OptionsStruct.HessFcn) || ~isempty(OptionsStruct.HessMult)
                OptionsStruct.Hessian = 'user-supplied';
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
                obj = optim.options.fmincon2;
                
                % Call the superclass method to upgrade the object
                obj = upgradeFrom13a(obj, s); 
                
                % The SolverVersion property was not present in 13a. We
                % clear it here and the remainer of loadobj will set it
                % correctly.
                obj.FminconVersion = [];
                
            end
            
            % Upgrading to 13b
            % Changing default FinDiffRelStep to a string for improved
            % display
            % Addition of TolGradCon option for interior-point
            if obj.Version < 2 && ...
                    ~ischar(obj.OptionsStore.AlgorithmDefaults{1}.FinDiffRelStep)
                
                % Change default for FinDiffRelStep
                obj.OptionsStore.AlgorithmDefaults{1}.FinDiffRelStep = 'sqrt(eps)';
                obj.OptionsStore.AlgorithmDefaults{2}.FinDiffRelStep = 'sqrt(eps)';
                obj.OptionsStore.AlgorithmDefaults{3}.FinDiffRelStep = 'sqrt(eps)';
                obj.OptionsStore.AlgorithmDefaults{4}.FinDiffRelStep = 'sqrt(eps)';
                
                % If FinDiffRelStep has not been set by user then change
                % option value to the default string
                if ~obj.OptionsStore.SetByUser.FinDiffRelStep
                    obj.OptionsStore.Options.FinDiffRelStep = 'sqrt(eps)';
                end
                
                % Add TolGradCon
                obj.OptionsStore.AlgorithmDefaults{2}.TolGradCon = 1.0000e-06;
                obj.OptionsStore.SetByUser.TolGradCon = false;
                obj.OptionsStore.IsConstantDefault.TolGradCon = true;
                obj.OptionsStore.Options.TolGradCon  = 1.0000e-06;                
            end
            
            % Upgrading to 14a
            % Change the algorithm default
            % Default for UseParallel now a boolean rather than 'never'
            if obj.Version < 2
                % Change the default in the OptionsStore
                obj.OptionsStore.DefaultAlgorithm = 'interior-point';
                
                % If the user hasn't set the Algorithm option, keep the
                % saved value of Algorithm.
                if ~obj.OptionsStore.SetByUser.Algorithm
                    obj = setPropertyNoChecks(obj, ...
                        'Algorithm', 'trust2-region-reflective');
                end
                
                % Change the default for UseParallel to false
                obj.OptionsStore.AlgorithmDefaults{1}.UseParallel= false;
                obj.OptionsStore.AlgorithmDefaults{2}.UseParallel= false;
                obj.OptionsStore.AlgorithmDefaults{3}.UseParallel= false;
                obj.OptionsStore.AlgorithmDefaults{4}.UseParallel= false;
                
            end
            
            % Upgrade to 16a
            if isempty(obj.FminconVersion) || obj.FminconVersion < 4
                
                % Add the HessFcn option to the OptionsStore
                obj.OptionsStore.AlgorithmDefaults{4}.HessFcn = [];
                obj.OptionsStore.SetByUser.HessFcn = false;
                obj.OptionsStore.IsConstantDefault.HessFcn = true;
                obj.OptionsStore.Options.HessFcn = [];
                
                % Add TolFunValue
                obj.OptionsStore.AlgorithmDefaults{1}.TolFunValue = 1e-6;
                obj.OptionsStore.AlgorithmDefaults{4}.TolFunValue = 1e-6;
                obj.OptionsStore.IsConstantDefault.TolFunValue = true;
                % Set TolFunValue to whatever of TolFun was saved, but only
                % if the selected algorithm has "FunctionTolerance".
                % Otherwise, set to its default value for another algorithm
                % that does use "FunctionTolerance"
                if isfield(obj.OptionsStore.AlgorithmDefaults{obj.OptionsStore.AlgorithmIndex},'TolFunValue') && obj.OptionsStore.SetByUser.TolFun
                    obj.OptionsStore.SetByUser.TolFunValue = obj.OptionsStore.SetByUser.TolFun;
                    obj.OptionsStore.Options.TolFunValue = obj.OptionsStore.Options.TolFun;
                else
                    obj.OptionsStore.SetByUser.TolFunValue = false;
                    obj.OptionsStore.Options.TolFunValue = 1e-6;
                end
                % Objects prior to 15b are missing display-related fields
                % in OptionsStore
                obj.OptionsStore = optim.options.generateMultiAlgorithmDisplayOptions2( ...
                    obj.OptionsStore,'optim.options.fmincon2');
            end            
            
            % Set the FminconVersion number
            obj.FminconVersion = 4;            
            
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
%   the full OptionsStore. See below for an example for fmincon2.

% Define the algorithm names
OS.AlgorithmNames = {'active-set', 'interior-point', 'sqp', 'trust2-region-reflective'};

% Define the default algorithm
OS.DefaultAlgorithm = 'interior-point';

% Define the defaults for each algorithm
% Active-set
OS.AlgorithmDefaults{1}.DerivativeCheck= 'off';
OS.AlgorithmDefaults{1}.Diagnostics= 'off';
OS.AlgorithmDefaults{1}.DiffMaxChange= Inf;
OS.AlgorithmDefaults{1}.DiffMinChange= 0;
OS.AlgorithmDefaults{1}.Display = 'final';
OS.AlgorithmDefaults{1}.FinDiffRelStep = 'sqrt(eps)';
OS.AlgorithmDefaults{1}.FinDiffType= 'forward';
OS.AlgorithmDefaults{1}.FunValCheck= 'off';
OS.AlgorithmDefaults{1}.GradConstr= 'off';
OS.AlgorithmDefaults{1}.GradObj= 'off';
OS.AlgorithmDefaults{1}.OutputFcn= [];
OS.AlgorithmDefaults{1}.PlotFcns= [];
OS.AlgorithmDefaults{1}.TolFun= 1.0000e-06;
OS.AlgorithmDefaults{1}.TolFunValue= 1.0000e-06;
OS.AlgorithmDefaults{1}.TolCon= 1.0000e-06;
OS.AlgorithmDefaults{1}.TypicalX= 'ones(numberOfVariables,1)';
OS.AlgorithmDefaults{1}.UseParallel= false;
OS.AlgorithmDefaults{1}.MaxIter = 400;
OS.AlgorithmDefaults{1}.MaxFunEvals = '100*numberOfVariables';
OS.AlgorithmDefaults{1}.TolX = 1e-6;
OS.AlgorithmDefaults{1}.MaxSQPIter = '10*max(numberOfVariables,numberOfInequalities+numberOfBounds)';
OS.AlgorithmDefaults{1}.RelLineSrchBnd = [];
OS.AlgorithmDefaults{1}.RelLineSrchBndDuration = 1;
OS.AlgorithmDefaults{1}.TolConSQP = 1.0000e-06;
OS.AlgorithmDefaults{1}.NoStopIfFlatInfeas = 'off';
OS.AlgorithmDefaults{1}.PhaseOneTotalScaling = 'off';

% Interior-point
OS.AlgorithmDefaults{2}.DerivativeCheck= 'off';
OS.AlgorithmDefaults{2}.Diagnostics= 'off';
OS.AlgorithmDefaults{2}.DiffMaxChange= Inf;
OS.AlgorithmDefaults{2}.DiffMinChange= 0;
OS.AlgorithmDefaults{2}.Display = 'final';
OS.AlgorithmDefaults{2}.FinDiffRelStep = 'sqrt(eps)';
OS.AlgorithmDefaults{2}.FinDiffType= 'forward';
OS.AlgorithmDefaults{2}.FunValCheck= 'off';
OS.AlgorithmDefaults{2}.GradConstr= 'off';
OS.AlgorithmDefaults{2}.GradObj= 'off';
OS.AlgorithmDefaults{2}.OutputFcn= [];
OS.AlgorithmDefaults{2}.PlotFcns= [];
OS.AlgorithmDefaults{2}.TolFun= 1.0000e-06;
OS.AlgorithmDefaults{2}.TolCon= 1.0000e-06;
OS.AlgorithmDefaults{2}.TolGradCon = 1.0000e-06;
OS.AlgorithmDefaults{2}.TypicalX= 'ones(numberOfVariables,1)';
OS.AlgorithmDefaults{2}.UseParallel= false;
OS.AlgorithmDefaults{2}.MaxIter = 1000;
OS.AlgorithmDefaults{2}.MaxFunEvals = 3000;
OS.AlgorithmDefaults{2}.TolX = 1e-10;
OS.AlgorithmDefaults{2}.AlwaysHonorConstraints= 'bounds';
OS.AlgorithmDefaults{2}.HessFcn= [];
OS.AlgorithmDefaults{2}.Hessian= 'bfgs';
OS.AlgorithmDefaults{2}.HessMult= [];
OS.AlgorithmDefaults{2}.InitBarrierParam= 0.1000;
OS.AlgorithmDefaults{2}.InitTrustRegionRadius= 'sqrt(numberOfVariables)';
OS.AlgorithmDefaults{2}.MaxProjCGIter= '2*(numberOfVariables-numberOfEqualities)';
OS.AlgorithmDefaults{2}.ObjectiveLimit= -1.0000e+20;
OS.AlgorithmDefaults{2}.ScaleProblem= 'none';
OS.AlgorithmDefaults{2}.SubproblemAlgorithm= 'ldl-factorization';
OS.AlgorithmDefaults{2}.TolProjCG = 0.0100;
OS.AlgorithmDefaults{2}.TolProjCGAbs = 1.0000e-10;

% sqp
OS.AlgorithmDefaults{3}.DerivativeCheck= 'off';
OS.AlgorithmDefaults{3}.Diagnostics= 'off';
OS.AlgorithmDefaults{3}.DiffMaxChange= Inf;
OS.AlgorithmDefaults{3}.DiffMinChange= 0;
OS.AlgorithmDefaults{3}.Display = 'final';
OS.AlgorithmDefaults{3}.FinDiffRelStep = 'sqrt(eps)';
OS.AlgorithmDefaults{3}.FinDiffType= 'forward';
OS.AlgorithmDefaults{3}.FunValCheck= 'off';
OS.AlgorithmDefaults{3}.GradConstr= 'off';
OS.AlgorithmDefaults{3}.GradObj= 'off';
OS.AlgorithmDefaults{3}.OutputFcn= [];
OS.AlgorithmDefaults{3}.PlotFcns= [];
OS.AlgorithmDefaults{3}.TolFun= 1.0000e-06;
OS.AlgorithmDefaults{3}.TolCon= 1.0000e-06;
OS.AlgorithmDefaults{3}.TypicalX= 'ones(numberOfVariables,1)';
OS.AlgorithmDefaults{3}.UseParallel= false;
OS.AlgorithmDefaults{3}.MaxIter = 400;
OS.AlgorithmDefaults{3}.MaxFunEvals = '100*numberOfVariables';
OS.AlgorithmDefaults{3}.TolX = 1e-6;
OS.AlgorithmDefaults{3}.ObjectiveLimit= -1.0000e+20;
OS.AlgorithmDefaults{3}.ScaleProblem= 'none';

% trust2-region-reflective
OS.AlgorithmDefaults{4}.DerivativeCheck= 'off';
OS.AlgorithmDefaults{4}.Diagnostics= 'off';
OS.AlgorithmDefaults{4}.DiffMaxChange= Inf;
OS.AlgorithmDefaults{4}.DiffMinChange= 0;
OS.AlgorithmDefaults{4}.Display = 'final';
OS.AlgorithmDefaults{4}.FinDiffRelStep = 'sqrt(eps)';
OS.AlgorithmDefaults{4}.FinDiffType= 'forward';
OS.AlgorithmDefaults{4}.FunValCheck= 'off';
OS.AlgorithmDefaults{4}.GradConstr= 'off';
OS.AlgorithmDefaults{4}.GradObj= 'off';
OS.AlgorithmDefaults{4}.OutputFcn= [];
OS.AlgorithmDefaults{4}.PlotFcns= [];
OS.AlgorithmDefaults{4}.TolFun= 1.0000e-06;
OS.AlgorithmDefaults{4}.TolFunValue= 1.0000e-06;
OS.AlgorithmDefaults{4}.TolCon= 1.0000e-06;
OS.AlgorithmDefaults{4}.TypicalX= 'ones(numberOfVariables,1)';
OS.AlgorithmDefaults{4}.UseParallel= false;
OS.AlgorithmDefaults{4}.MaxIter = 400;
OS.AlgorithmDefaults{4}.MaxFunEvals = '100*numberOfVariables';
OS.AlgorithmDefaults{4}.TolX = 1e-6;
OS.AlgorithmDefaults{4}.Hessian = 'off';
OS.AlgorithmDefaults{4}.HessMult = [];
OS.AlgorithmDefaults{4}.HessFcn= [];
OS.AlgorithmDefaults{4}.HessPattern = 'sparse(ones(numberOfVariables))';
OS.AlgorithmDefaults{4}.MaxPCGIter = 'max(1,floor(numberOfVariables/2))';
OS.AlgorithmDefaults{4}.PrecondBandWidth = 0;
OS.AlgorithmDefaults{4}.TolPCG = 0.1000;

% Call the package function to generate the OptionsStore
OS = optim.options.generateMultiAlgorithmOptionsStore2(OS, 'optim.options.fmincon2');
end
