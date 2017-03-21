classdef (Abstract) MultiAlgorithm2 < optim.options.SolverOptions2
    %
    
    %MultiAlgorithm2 Options for Optimization Toolbox solvers with multiple
    %               algorithms
    %
    %   MultiAlgorithm2 is an abstract class representing a set of options
    %   for an Optimization Toolbox solver, where the solver has multiple
    %   algorithms. You cannot create instances of this class directly. You
    %   must create an instance of a derived class such
    %   optim.options.fmincon2.
    %
    %   See also OPTIM.OPTIONS.SolverOptions2
    
    %   Copyright 2012-2015 The MathWorks, Inc.
    
    properties (Dependent)
        
        %ALGORITHM Choose the optimization algorithm
        %
        %   For more information, see the "Options" section documentation page of
        %   the solver you are using.
        Algorithm
    end
    
    % Constructor
    methods (Hidden)
        function obj = MultiAlgorithm2(varargin)
            %MultiAlgorithm2 Construct a new MultiAlgorithm2 object
            %
            %   OPTS = OPTIM.OPTIONS.MultiAlgorithm2 creates a set of options with each
            %   option set to its default value.
            %
            %   OPTS = OPTIM.OPTIONS.MultiAlgorithm2(PARAM, VAL, ...) creates a set of
            %   options with the named parameters altered with the specified values.
            %
            %   OPTS = OPTIM.OPTIONS.MultiAlgorithm2(OLDOPTS, PARAM, VAL, ...) creates
            %   a copy of OLDOPTS with the named parameters altered with the specified
            %   value
            
            % Call the superclass constructor
            obj = obj@optim.options.SolverOptions2(varargin{:});
            
        end
    end
    
    
    methods
        function obj = set.Algorithm(obj, value)

            % Map old algorithm names to new values, e.g.
            % 'trust2-region-reflective' for fsolve maps to 'trust2-region'
            value = mapAlgorithmName(obj, value);
            
            % Set Algorithm option
             obj = setProperty(obj, 'Algorithm', value, ...
                obj.OptionsStore.AlgorithmNames);            

            % If we get here, the property set has been successful and we
            % can update the OptionsStore
            
            % Set the algorithm index.
            obj.OptionsStore.AlgorithmIndex = strcmp(...
                obj.OptionsStore.AlgorithmNames, obj.OptionsStore.Options.Algorithm);            
            
            % Set the value of those options that have non-constant
            % defaults and haven't been set by the user
            for i = 1:length(obj.OptionsStore.NonConstantDefaultFields)
                if ~obj.OptionsStore.SetByUser.(obj.OptionsStore.NonConstantDefaultFields{i})
                    obj.OptionsStore.Options.(obj.OptionsStore.NonConstantDefaultFields{i}) = ...
                        obj.OptionsStore.NonConstantDefaults.(obj.OptionsStore.NonConstantDefaultFields{i}){obj.OptionsStore.AlgorithmIndex};
                end
            end
        end
        
        function value = get.Algorithm(obj)
            value = obj.OptionsStore.Options.Algorithm;
        end
        
    end
    
    % Display methods
    methods(Hidden, Access = protected)
        
        function header = addExtraHeader(obj, header)
            %ADDEXTRAHEADER Append text to the header
            %
            %   HEADER = ADDEXTRAHEADER(OBJ, HEADER) appends the following
            %   to the header in the following order:-
            %
            %     * The algorithm header
            %     * "Set by user" section of the display if no options are
            %     set by user.
            
            % Get algorithm header and append it to the main header
            algHeader = getAlgorithmHeader(obj);
            header = sprintf('%s%s', header, algHeader);
            
            % Call superclass addExtraHeader to append the extra header to
            % main header
            header = addExtraHeader@optim.options.SolverOptions2(obj, header);
            
        end
        
        function footer = addExtraFooter(obj, footer)
            %ADDEXTRAFOOTER Add text to the footer
            %
            %   HEADER = ADDEXTRAFOOTER(OBJ, HEADER) prepends the following
            %   to the footer in the following order:-
            %
            %      * "Default" section of the display if all options are
            %      set by user.
            %      * The algorithm footer
            
            % Call superclass addExtraFooter to prepend the extra footer to
            % main footer
            footer = addExtraFooter@optim.options.SolverOptions2(obj, footer);
            
            % Get algorithm footer and prepend it to the main footer
            algFooter = getAlgorithmFooter(obj);
            footer = sprintf('%s%s', footer, algFooter);
            
        end
        
        function algHeader = getAlgorithmHeader(obj)
            %GETALGORITHMHEADER Return the algorithm header for the display
            %
            %   HEADER = GETALGORITHMHEADER(OBJ) returns an algorithm specific string
            %   to be added to the header. An example of this string for fmincon2 is
            %   shown below:
            %
            %   "Options used by current Algorithm ('trust2-region-reflective'):
            %    (Other available algorithms: 'active-set', 'interior-point', 'sqp')"
            
            
            % Now add the algorithm display header
            if isscalar(obj)
                
                %%% Create the first line of the algorithm display header
                currAlgStr = getString(message('optimlib:options:MultiAlgorithm2:CurrentAlgorithmStr', ...
                    generateAlgorithmString, obj.Algorithm));
                algLine1 = sprintf('\n   %s', currAlgStr);
                
                %%% Create the second line of the algorithm display header
                
                % Get a cell array of the other algorithm names
                otherAlgNames = obj.OptionsStore.AlgorithmNames(~obj.OptionsStore.AlgorithmIndex);
                
                % Put the other algorithm names in alphabetical order
                otherAlgNames = sort(otherAlgNames);
                
                % Create a comma-separated list of the other algorithm
                % names. We considered using cellfun to generate this list
                % but it is 3x slower than a for loop we also used at one
                % time. We use strjoin since it's built-in and requires
                % less code. Note: the surrounding quotes add the single
                % quote to the first and last algorithm name.
                otherAlgNamesStr = [ '''', strjoin(otherAlgNames, ''', '''), ''''];
                
                % Get the "other available algorithms line of the header
                otherAvailAlgsStr = addLink2( ...
                    getString(message('optimlib:options:MultiAlgorithm2:OtherAvailableAlgorithmsStr')), ...
                    'optim','helptargets.map','choose_algorithm',false);
                algLine2 = sprintf('   (%s: %s)', otherAvailAlgsStr,otherAlgNamesStr);
                
                %%% Create the full algorithm header
                algHeader = sprintf('%s\n%s\n', algLine1, algLine2);
                
            else
                algHeader = '';
            end
            
        end
        
        function algFooter = getAlgorithmFooter(obj)
            %GETALGORITHMFOOTER Return the algorithm footer for the display
            %
            %   HEADER = GETALGORITHMHEADER(OBJ) returns an algorithm specific string
            %   to be added to the footer. An example of this string for fmincon2 is
            %   shown below:
            %
            %   "Show options not used by current Algorithm('trust2-region-reflective')"
            
            % We're only customizing the display for scalar objects
            if isscalar(obj) && ~algorithmHasAllDisplayOptions(obj)
                
                % Create the algorithm footer
                
                % Below, we use evalc to capture the structure display for
                % the options that aren't used by the algorithm. evalc
                % always sets hotlinks to true (see g591312 for the
                % rationale for this). As such we need to capture the
                % current strong start and end tags and pass them to
                % showOptionsUnusedByAlgorithm.
                currentStrongStartTag = optim.options.SolverOptions2.generateStrongStartTag; %#ok 
                currentStrongEndTag = optim.options.SolverOptions2.generateStrongEndTag; %#ok
                
                % Generate the "Show options" link. Not keen on using evalc
                % here, but structs do not have a "toString" method. Note
                % that we hope to re-implement this in the future.
                showOptsStr = evalc('showOptionsUnusedByAlgorithm(obj, currentStrongStartTag, currentStrongEndTag)');
                if optim.options.SolverOptions2.enableLinks
                    showOptsStr = regexprep(showOptsStr, '''', '''''');
                    showOptsStr = regexprep(showOptsStr, '\n', '\\n');
                    linkCmdStartTag = sprintf('<a href="matlab: fprintf(''%s'')">', ...
                        showOptsStr);
                    linkCmdEndTag = '</a>';
                    algFooterStr = getString(message('optimlib:options:MultiAlgorithm2:ShowUnusedOptionsLinks', ...
                        linkCmdStartTag, linkCmdEndTag, ...
                        generateAlgorithmString, obj.Algorithm));
                    algFooter = sprintf('   %s\n', algFooterStr);
                else
                    algFooterStr = getString(message('optimlib:options:MultiAlgorithm2:ShowUnusedOptionsNoLinks', 'Algorithm', obj.Algorithm));
                    algFooter = sprintf('   %s\n%s', algFooterStr, showOptsStr);
                end
                
                % Add a new line if all options are set by the user
                if needExtraFooter(obj)
                    algFooter = sprintf('\n%s', algFooter);
                end
                
            else
                algFooter = '';
            end
            
        end
        
        function allOptionsAtDefault = needExtraHeader(obj)
            %NEEDEXTRAHEADER Determine whether an extra header is needed
            %
            %   ALLOPTIONSATDEFAULT = NEEDEXTRAHEADER(OBJ) returns whether an extra
            %   header is needed.
            
            allOptionsAtDefault = (numCurrentAlgorithmDisplayOptionsSetByUser(obj) == 0);
        end
        
        function allOptionsSetByUser = needExtraFooter(obj)
            %NEEDEXTRAFOOTER Determine whether an extra footer is needed
            %
            %   ALLOPTIONSSETBYUSER = NEEDEXTRAFOOTER(OBJ) returns whether an extra
            %   footer is needed.
            
            allOptionsSetByUser = (numCurrentAlgorithmDisplayOptionsSetByUser(obj) == numCurrentAlgorithmDisplayOptions(obj));
        end
        
        function currentAlgorithmOptions = getDisplayOptions(obj)
            %GETDISPLAYOPTIONS Get the options to be displayed
            %
            %   CURRENTALGORITHMOPTIONS = GETDISPLAYOPTIONS(OBJ) returns a cell array of
            %   options to be displayed. For solver objects that inherit from
            %   MultiAlgorithm2, options for the current algorithm are displayed.
            
            currentAlgorithmOptions = obj.OptionsStore.DisplayOptions{obj.OptionsStore.AlgorithmIndex};           
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
            
            % Replace any special strings in the options object
            obj = replaceSpecialStrings(obj);
            
            % Extract main options structure
            OptionsStruct = obj.OptionsStore.Options;            
            
            % optimoptions2 handles options with defaults that are algorithm
            % dependent. Some of these options may not be defined for the current
            % algorithm. If this is the case, we set the option to empty in the
            % options structure, which is the optimset standard.
            for i = 1:length(obj.OptionsStore.NonConstantDefaultFields)
                % Get the option name
                thisOption = obj.OptionsStore.NonConstantDefaultFields{i};
                % Is the default for thisOption defined for the current
                % algorithm? If it is not defined and the option has not
                % been set by the user, set the value to empty in the
                % structure.
                if ~isSetByUser(obj, thisOption) && ~isfield(obj.OptionsStore.AlgorithmDefaults{obj.OptionsStore.AlgorithmIndex}, thisOption)
                    OptionsStruct.(thisOption) = [];
                end
            end
            
            % Call method to map optimoptions2 for use in the solver and
            % optimtool
            OptionsStruct = mapOptionsForSolver(obj, OptionsStruct);            
            
        end
    end
    
    methods (Access = private)
        
        function numOptions = numCurrentAlgorithmDisplayOptions(obj)
            %NUMCURRENTALGORITHMDISPLAYOPTIONS Number of current algorithm
            %   options that are displayed
            %
            %   NUMOPTIONS = NUMCURRENTALGORITHMDISPLAYOPTIONS(OBJ) returns
            %   the number of options for the current algorithm that are
            %   displayed.
            
            numOptions = obj.OptionsStore.NumDisplayOptions(obj.OptionsStore.AlgorithmIndex);
            
        end
        
        function numSetByUser = numCurrentAlgorithmDisplayOptionsSetByUser(obj)
            %NUMCURRENTALGORITHMDISPLAYOPTIONSSETBYUSER Number of current
            %   algorithm display options that are set by the user
            %
            %   NUMSETBYUSER = NUMCURRENTALGORITHMDISPLAYOPTIONSSETBYUSER(OBJ) 
            %   returns the number of displayed options for the current 
            %   algorithm that are set by the user
            
            % Get names of all the properties for the current algorithm
            allOptions = getDisplayOptions(obj);
            numSetByUser = 0;
            for i = 1:length(allOptions)
                if obj.OptionsStore.SetByUser.( ...
                        optim.options.OptionAliasStore2.getAlias(allOptions{i}, obj.SolverName, obj.OptionsStore.Options))
                    numSetByUser = numSetByUser + 1;
                end
            end
        end
        
        function hasAllOptions = algorithmHasAllDisplayOptions(obj)
%ALGORITHMHASALLDISPLAYOPTIONS Determine whether the algorithm has
%   all the options that the solver can display
%
%   HASALLOPTIONS = ALGORITHMHASALLDISPLAYOPTIONS(OBJ) returns true if the
%   set of options displayed for the current algorithm is identical to the
%   set of options that can be displayed by the solver.
            
            numAlgorithmOptions = numCurrentAlgorithmDisplayOptions(obj);
            numSolverOptions = length(properties(obj));
            hasAllOptions = (numAlgorithmOptions == numSolverOptions);
        end
    end
    
    methods (Access = private)
        
        function [setByUser, default] = getOptionsUnusedByAlgorithm(obj)
            %GETOPTIONSUNUSEDBYALGORITHM Return options unused by algorithm
            %
            %   [SETBYUSER, DEFAULT] = GETOPTIONSUNUSEDBYALGORITHM(OBJ) returns the
            %   options that are not used by the current algorithm. The cell array
            %   SETBYUSER returns the unused options that have been set by the user.
            %   The cell array, DEFAULT, returns the unused options that are at their
            %   default values.
            
            % Get the public (i.e. non-hidden) options for the solver
            allOptions = properties(obj);
            
            % Get the options for the current algorithm
            currAlgOptions = getDisplayOptions(obj);
            
            % Determine the options not used by the current algorithm
            allUnusedOptions = setdiff(allOptions, currAlgOptions);
            
            % Create the setByUser and default structures
            setByUser = struct([]);
            default = struct([]);
            for i = 1:length(allUnusedOptions)
                if obj.OptionsStore.SetByUser.(optim.options.OptionAliasStore2.getAlias(allUnusedOptions{i}, obj.SolverName, obj.OptionsStore.Options))
                    setByUser(1).(allUnusedOptions{i}) = obj.(allUnusedOptions{i});
                else
                    default(1).(allUnusedOptions{i}) = obj.(allUnusedOptions{i});
                end
            end
            
        end
    end
    
    methods
        function showOptionsUnusedByAlgorithm(obj, StrongStartTag, StrongEndTag)
            %SHOWOPTIONSUNUSEDBYALGORITHM Show options unused by algorithm
            %
            %   SHOWOPTIONSUNUSEDBYALGORITHM(OBJ, STRONGSTARTTAG,
            %   STRONGENDTAG) displays the options that are not used by the
            %   current algorithm.
            %
            %   See also optimoptions2
            
            % Note that this method is called via evalc to capture the
            % structure display for the options that aren't used by the
            % algorithm. evalc always sets hotlinks to true (see g591312
            % for the rationale for this). As such this method requires the
            % tags from the calling workspace to be passed in, as these
            % tags reflect the hotlinks setting in the (actual) calling
            % workspace and not evalc's workspace.
                
            % Get the options unused by the current algorithm
            [setByUserOpts, defaultOpts] = getOptionsUnusedByAlgorithm(obj);
            
            % Display the unused options that have been set by the user
            if ~isempty(setByUserOpts)
                fprintf('   %s\n', ...
                    getString(message('optimlib:options:SolverOptions2:SetByUserHeader', ...
                    StrongStartTag, StrongEndTag)));
                disp(setByUserOpts);
            end
            
            % Display the unused options that take their default value
            if ~isempty(defaultOpts)
                fprintf('   %s\n', getString(message('optimlib:options:SolverOptions2:DefaultHeader', ...
                    StrongStartTag, StrongEndTag)));
                disp(defaultOpts);
            end
            
        end
        
    end
    
    % Public informal interface methods
    methods
        
        function obj = resetoptions(obj, name)
%RESETOPTIONS Reset optimization options
%
%   OPTS = RESETOPTIONS(OPTS,'option') resets the specified option back to
%   its default value.
%
%   OPTS = RESETOPTIONS(OPTS,OPTIONS) resets more than one option at a time
%   when OPTIONS is a cell array of strings.  
%
%   See also optimoptions2.             

            % Handle non-cell inputs
            if ~iscell(name)
                name = {name};                
            end
            
            % Look to see if Algorithm is being reset
            idxAlgorithm = strcmp(name, 'Algorithm');
            
            % Reset Algorithm if required
            if any(idxAlgorithm)
                
                % Remove Algorithm from list of names
                name(idxAlgorithm) = [];
                
                % Set Algorithm back to default 
                obj.Algorithm = obj.OptionsStore.DefaultAlgorithm;
                
                % Algorithm has not been set by user
                obj.OptionsStore.SetByUser.Algorithm = false;
                
            end
            
            % Loop through each option and reset to default
            for i = 1:numel(name)
                
                % Get the default value for the option
                OptionStoredName = optim.options.OptionAliasStore2.getAlias(name{i}, obj.SolverName, obj.OptionsStore.Options);
                optionDefault = optim.options.getOptionDefaultValue2(obj.OptionsStore, OptionStoredName);

                % Reset the option
                obj.OptionsStore.SetByUser.(OptionStoredName) = false;
                obj.OptionsStore.Options.(OptionStoredName) = optionDefault;
            end
            
        end
        
    end        
     
    methods (Access = protected)
        
        function name = mapAlgorithmName(~, name)
%MAPALGORITHMNAME Map old algorithm name to current one
%
%   NAME = MAPALGORITHMNAME(OBJ, OLDNAME) maps a previous name for an
%   algorithm to its current value.
%
%   This method does not perform any mapping. Subclasses can perform their
%   specific mapping.
            
        end
        
    end
end


function algStr = generateAlgorithmString

StartTag = optim.options.SolverOptions2.generateStrongStartTag;
EndTag = optim.options.SolverOptions2.generateStrongEndTag;
algStr = sprintf('%sAlgorithm%s', StartTag, EndTag);
                
end
