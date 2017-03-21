function [state,optnew,updateoptions]  = gatooloutput2(optlast,state,flag)
%GATOOLOUTPUT2 private to OPTIMTOOL. 

%   Copyright 2004-2015 The MathWorks, Inc.

%Initialize
updateoptions = false;
optnew = optlast;
drawnow;
STOP = 0;
RUN_RESUME = 1;
PAUSE = 2;
modelCheck =  false;

switch flag
    case 'init'
        optimtoolGui = javaMethodEDT('getOptimGUI','com.mathworks.toolbox.optim.OptimGUI'); % Get a handle to the GUI
        [msg,id] = getLastWarnForOptimtool;
        if ~isempty(msg)
            % Make sure that exactly one newline character is at the end of 'msg'
            newLineIndex = regexp(msg,'\n');
            if isempty(newLineIndex) || newLineIndex(end) ~=length(msg)
                msg = sprintf('%s\n',msg);
            end
            javaMethodEDT('appendResults',optimtoolGui,['Warning: ',msg]); % Append to the 'Status and Results' panel
        end
        setappdata(0,'last_warning_id_for_optimtool',id);
        %Set up for the GUI.
        return; %Nothing to do now
    case 'iter'
        optimtoolGui = javaMethodEDT('getOptimGUI','com.mathworks.toolbox.optim.OptimGUI');
        if isempty(optimtoolGui)
            state.StopFlag = 'Abnormal termination';
            return;
        end
        [msg,id] = getLastWarnForOptimtool;
        if ~isempty(msg) && ~strcmp(id,getappdata(0,'last_warning_id_for_optimtool'))
            newLineIndex = regexp(msg,'\n');
            if isempty(newLineIndex) || newLineIndex(end) ~=length(msg)
                msg = sprintf('%s\n',msg);
            end
            javaMethodEDT('appendResults',optimtoolGui,['Warning: ',msg]);
            lastwarn('');
        end
        javaMethodEDT('setIteration',optimtoolGui,value2RHS2(state.Generation));
        %Action based on run mode of GUI
        RunMode = javaMethodEDT('getRunMode',optimtoolGui);
        switch RunMode
            case RUN_RESUME
            % Nothing to do                
            case STOP   %Stop
                state.StopFlag = 'Stop requested';
                return;
            case PAUSE  %Pause
                fprintf('%s\n%s\n','OPTIMTOOL is paused. MATLAB Command prompt will', ...
                    'not be accessible until the genetic algorithm solver is completed.');
                PauseTimeStart = tic; 
                %If in pause state keeping looping here.
                while true
                    drawnow
                    if isempty(javaMethodEDT('getOptimGUI','com.mathworks.toolbox.optim.OptimGUI'))
                        state.StopFlag = 'Abnormal termination';
                        return;
                    end
                    mode = javaMethodEDT('getRunMode',optimtoolGui);
                    if mode == STOP
                        state.StopFlag = 'Stop requested';
                        return;
                    elseif mode == RUN_RESUME
                        modelCheck = true;
                        break;
                    end
                    pause(0.1);
                end % End while
                PauseTimeTotal  = (tic - PauseTimeStart); %Total time the gui is paused
                %Pause state of the GUI is not included in stall tolerances.
                state.LastImprovementTime = state.LastImprovementTime + PauseTimeTotal;
                state.StartTime = state.StartTime + PauseTimeTotal;
            otherwise
                return;
        end
        
        if  javaMethodEDT('getChangedState',optimtoolGui) || modelCheck
            hashOptions = javaMethodEDT('getChangedOptionModelAndClear',optimtoolGui);
            hashProblem = javaObjectEDT('java.util.Hashtable');
            [~,optnew,~,errOpt] = readOptimHashTable(hashProblem, hashOptions);
            try 
                if ~isempty(errOpt)
                    error(message('globaloptim:gatooloutput2:invalidModel', errOpt));
                end
                
                % Attempt to update changed options
                if state.IsMixedInteger
                    extradefault = {'PopInitRange', [-1e4 + 1; 1e4 + 1]};
                else
                    extradefault = {'PopInitRange', [-10;10]};
                end  
                validatefcn = @(optnew2, optnew)i_validatefcn(optnew2, optnew, state);
                optnew = optimtoolUpdateChangedOptions2(...
                    optnew, optlast, 'ga2', extradefault, validatefcn);
                
                % Indicate that the options have been updated
                updateoptions = true;
            catch ME
                errordlg(ME.message,'GA2 run time error');    
                javaMethodEDT('setRunMode',optimtoolGui,2) %Set GUI to pause.
                optnew = optlast;  %Don't change the options if it contains errors
            end
        else
            optnew = optlast;
        end
    case 'interrupt'
        optimtoolGui = javaMethodEDT('getOptimGUI','com.mathworks.toolbox.optim.OptimGUI');
        if isempty(optimtoolGui)
            state.StopFlag = 'Abnormal termination';
            return;
        end
        [msg,id] = getLastWarnForOptimtool;
        if ~isempty(msg) && ~strcmp(id,getappdata(0,'last_warning_id_for_optimtool'))
            newLineIndex = regexp(msg,'\n');
            if isempty(newLineIndex) || newLineIndex(end) ~=length(msg)
                msg = sprintf('%s\n',msg);
            end
            javaMethodEDT('appendResults',optimtoolGui,['Warning: ',msg]);
            lastwarn('');
        end
        %Action based on run mode of GUI
        RunMode = javaMethodEDT('getRunMode',optimtoolGui);
        switch RunMode
            case RUN_RESUME
            % Nothing to do
            case STOP   %Stop
                state.StopFlag = 'Stop requested';
                return;
            case PAUSE  %Pause
                fprintf('%s\n%s\n','OPTIMTOOL is paused. MATLAB Command prompt will', ...
                    'not be accessible until the genetic algorithm solver is completed.');
                PauseTimeStart =cputime;
                %If in pause state keeping looping here.
                while true
                    drawnow
                    if isempty(javaMethodEDT('getOptimGUI','com.mathworks.toolbox.optim.OptimGUI'))
                        state.StopFlag = 'Abnormal termination';
                        return;
                    end
                    mode = javaMethodEDT('getRunMode',optimtoolGui);
                    if mode == STOP
                        state.StopFlag = 'Stop requested';
                        return;
                    elseif mode == RUN_RESUME
                        break;
                    end
                    pause(0.1);
                end % End while
                PauseTimeTotal  = (cputime - PauseTimeStart); %Total time the gui was paused
                %Pause state of the GUI is not included in stall tolerances.
                state.LastImprovementTime = state.LastImprovementTime + PauseTimeTotal;
                state.StartTime = state.StartTime + PauseTimeTotal;
        end
    case 'done'
        optimtoolGui = javaMethodEDT('getOptimGUI','com.mathworks.toolbox.optim.OptimGUI');
        if isempty(optimtoolGui)
            state.StopFlag = 'Abnormal termination';
            return;
        end
         javaMethodEDT('setIteration',optimtoolGui,value2RHS2(state.Generation));
        %If using a hybrid function, give a message.
        if ~isempty(optlast.HybridFcn)
           javaMethodEDT('appendResults',optimtoolGui,sprintf('%s\n','Switching to hybrid function.'));
        end
        [msg,id] = getLastWarnForOptimtool;
        if ~isempty(msg) && ~strcmp(id,getappdata(0,'last_warning_id_for_optimtool'))
            javaMethodEDT('appendResults',optimtoolGui,['Warning: ',lastwarn]);
            lastwarn('');
        end
        if isappdata(0,'last_warning_id_for_optimtool')
            rmappdata(0,'last_warning_id_for_optimtool');
        end
end 


function optnewstructwithdefault = i_validatefcn(...
    optnewstructwithdefault, optnewstruct, state)

type = state.LinearConstrType;
optnewstructwithdefault.MultiObjective = false;
optnewstructwithdefault = validate2(optnewstructwithdefault,type,size(state.Population,2),[],[],optnewstruct);
