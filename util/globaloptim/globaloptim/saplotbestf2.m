function stop = saplotbestf2(~,optimvalues,flag)
%SAPLOTBESTF2 PlotFcn to plot best function value.
%   STOP = SAPLOTBESTF2(OPTIONS,OPTIMVALUES,FLAG) where OPTIMVALUES is a
%   structure with the following fields:
%              x: current point
%           fval: function value at x
%          bestx: best point found so far
%       bestfval: function value at bestx
%    temperature: current temperature
%      iteration: current iteration
%      funccount: number of function evaluations
%             t0: start time
%              k: annealing parameter
%
%   OPTIONS: The options object created by using optimoptions2
%
%   FLAG: Current state in which PlotFcn is called. Possible values are:
%           init: initialization state
%           iter: iteration state
%           done: final state
%
%   STOP: A boolean to stop the algorithm.
%
%   Example:
%    Create an options structure that will use SAPLOTBESTF2
%    as the plot function
%     options = optimoptions2('simulannealbnd2','PlotFcn',@saplotbestf2);

%   Copyright 2006-2015 The MathWorks, Inc.

persistent thisTitle

stop = false;
switch flag
    case 'init'
        plotBest = plot(optimvalues.iteration,optimvalues.bestfval, '.b');
        set(plotBest,'Tag','saplotbestf2');
        xlabel('Iteration','interp','none');
        ylabel('Function value','interp','none')
        thisTitle = title(sprintf('Best Function Value: %g',optimvalues.bestfval),'interp','none');
    case 'iter'
        plotBest = findobj(get(gca,'Children'),'Tag','saplotbestf2');
        newX = [get(plotBest,'Xdata') optimvalues.iteration];
        newY = [get(plotBest,'Ydata') optimvalues.bestfval];
        set(plotBest,'Xdata',newX, 'Ydata',newY);
        if isempty(thisTitle)
            set(get(gca,'Title'),'String',sprintf('Best Function Value: %g',optimvalues.bestfval));
        else
            set(thisTitle,'String',sprintf('Best Function Value: %g',optimvalues.bestfval));
        end

end
