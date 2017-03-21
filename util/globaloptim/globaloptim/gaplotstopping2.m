function state = gaplotstopping2(options,state,flag)
%GAPLOTSTOPPING2 Display stopping criteria levels.
%   STATE = GAPLOTSTOPPING2(OPTIONS,STATE,FLAG) plots the current percentage
%   of the various criteria for stopping.  
%
%   Example:
%    Create an options structure that uses GAPLOTSTOPPING2
%    as the plot function
%      options = optimoptions2('ga2','PlotFcn',@gaplotstopping2);
%
%    (Note: If calling gamultiobj2, replace 'ga2' with 'gamultiobj2') 

%   Copyright 2003-2015 The MathWorks, Inc. 

CurrentGen = state.Generation;
stopCriteria(1) = CurrentGen / options.MaxGenerations;
stopString{1} = 'Generation';
stopCriteria(2) = toc(state.StartTime) / options.MaxTime;
stopString{2} = 'Time';

if isfield(state,'LastImprovement')
    stopCriteria(3) = (CurrentGen - state.LastImprovement) / options.MaxStallGenerations;
    stopString{3} = 'Stall (G)';
end

if isfield(state,'LastImprovementTime')
    stopCriteria(4) = toc(state.LastImprovementTime) / options.MaxStallTime;
    stopString{4} = 'Stall (T)';
end

ydata = 100 * stopCriteria;
switch flag
    case 'init'
        barh(ydata,'Tag','gaplotstopping2')
        set(gca,'xlim',[0,100],'yticklabel', ...
            stopString,'climmode','manual')
        xlabel('% of criteria met','interp','none')
        title('Stopping Criteria','interp','none')
    case 'iter'
        ch = findobj(get(gca,'Children'),'Tag','gaplotstopping2');
        set(ch,'YData',ydata);
end

