function [state,options,optchanged] = gaoutput2(~,options,state,flag)
%GAOUTPUT2 Helper function that manages the output functions for GA2.
%
%   [STATE, OPTIONS, OPTCHANGED] = GAOUTPUT2(~, OPTIONS, STATE, FLAG)
%   runs each of the display functions in the options.OutputFcns cell
%   array.
%
%   This is a helper function called by ga2 between each generation, and is
%   not typically called directly.

%   Copyright 2003-2015 The MathWorks, Inc.


% get the functions and return if there are none
optchanged = false;
functions = options.OutputFcns;
if(isempty(functions))
    return
end

solverName = options.OutputPlotFcnOptions.SolverName;
% call each output function
stopFlag = state.StopFlag;
args = options.OutputFcnsArgs;
for i = 1:length(functions)
    % Always clear state.StopFlag before calling OutputFcns so that each
    % function can independently set a termination message.
    state.StopFlag = '';
    % Add following flags to state. These flags are needed by
    % gatooloutput2, the output function which manages the
    % stop/pause/restart functionality in optimtool for ga2. 
    if isfield(options,'LinearConstr')
        state.LinearConstrType = options.LinearConstr.type;
    else
        state.LinearConstrType = 'unconstrained';
    end
    state.IsMixedInteger = ...
        isfield(options, 'IntegerVars') && ~isempty(options.IntegerVars);
    [state,optnew,changed] = feval(functions{i},options.OutputPlotFcnOptions,state,flag,args{i}{:});
    if ~isempty(state.StopFlag)
        stopFlag = [stopFlag state.StopFlag ';'];
    end
    if changed %If changes are not duplicates, we will get all the changes
       % Keep LinearConstr out of options and accept new options
       LinearConstr = options.LinearConstr;
       type = LinearConstr.type;
       gLength = size(state.Population,2);
       isMultiObjective = options.MultiObjective;
       % Store integer constraints information
       [intcon, UserSpecPopInitRange, UserVectorized] = ...
           i_storeIntegerInfo(options);        
       
       % Merge new options with existing options structure
       optnew2 = prepareOptionsForSolver2(optnew, solverName);
       options = gaoptimset2(options,optnew2);
       options.MultiObjective = isMultiObjective;       
       options = validate2(options,type,gLength,[],[]); 
       
       options.LinearConstr = LinearConstr;
       % Perform integer constraint specific tasks
       options = i_resetIntegerInfo(options, intcon, ...
           UserSpecPopInitRange, UserVectorized, gLength, LinearConstr);
       optchanged = true;
       
       if strcmpi(solverName,'ga2')
           options.OutputPlotFcnOptions  = optimoptions2(@ga2);
       else % gamultiobj2
           options.OutputPlotFcnOptions  = optimoptions2(@gamultiobj2);           
       end
       options.OutputPlotFcnOptions  = copyForOutputAndPlotFcn(options.OutputPlotFcnOptions,options);
    end
end

state.StopFlag = stopFlag;

%--------------------------------------------------------------------------
function [intcon, UserSpecPopInitRange, UserVectorized] = ...
    i_storeIntegerInfo(options)

intcon = options.IntegerVars;
% Keep UserSpecPopInitRange and UserVectorized out of options
UserSpecPopInitRange = [];
UserVectorized = [];
if isfield(options,'UserSpecPopInitRange')
    UserSpecPopInitRange = options.UserSpecPopInitRange;
end
if isfield(options,'UserVectorized')
    UserVectorized = options.UserVectorized;
end

%--------------------------------------------------------------------------
function options = i_resetIntegerInfo(options, intcon, ...
    UserSpecPopInitRange, UserVectorized, gLength, LinearConstr)

% Add IntegerVars field back to options structure
options.IntegerVars = intcon;

% The penalty reformulation overwrites options, restore them if needed.
if ~isempty(UserVectorized)
    % It is possible for users to alter the Vectorized option in an output
    % function. Say the user originally set the Vectorized option to 'off'
    % and they set it to 'on' in their output function. Because we set
    % options.Vectorized = 'on' for GA2 with a penalty reformulation, it
    % will look like this option has not changed. As such, we cannot tell
    % whether or not the user changed options.Vectorized to 'on' in an
    % output function. Note that we can determine if a user sets
    % options.Vectorized to 'off'.
    %
    % So, for penalty GA2, we ignore any change that a user makes to
    % the Vectorized option in an output function. We will use
    % gaminlpoverwriteoptions2 to overwrite the Vectorized option (to 'on')
    % and set the internal UserVectorized option to reflect the Vectorized
    % status when the user called GA2.
    %
    % gapenalty2 sets up UserVectorized from the value in the Vectorized
    % option, so we set options.Vectorized to reflect the value in
    % UserVectorized here.
    options.Vectorized = 'on';
    options.UserVectorized = UserVectorized;
    % Penalty function requires tournament selection.
	options.SelectionFcn = @selectiontournament2;
	options.SelectionFcnArgs = {2};
    
    % If there are any integer constraints, need to perform extra
    % validation checks
    if ~isempty(intcon)
        % Add UserSpecPopInitRange back to options
        options.UserSpecPopInitRange = UserSpecPopInitRange;
        % Validate2 mixed integer GA2 options
        gaminlpvalidateoptions2(options);
        % The call to gaoptimset2 removes the mixed integer GA2 specific
        % options. We need to recreate them here.
        options = gaminlpoverwriteoptions2(options, gLength, ...
            LinearConstr.lb, LinearConstr.ub, intcon);
    end
end