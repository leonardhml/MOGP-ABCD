function stop = psplotbestf2(optimvalues,flag)
%PSPLOTBESTF2 PlotFcn to plot best function value.
%   STOP = PSPLOTBESTF2(OPTIMVALUES,FLAG) where OPTIMVALUES is a structure
%   with the following fields:
%              x: current point X
%      iteration: iteration count
%           fval: function value
%       meshsize: current mesh size
%      funccount: number of function evaluations
%         method: method used in last iteration
%         TolFun: tolerance on function value in last iteration
%           TolX: tolerance on X value in last iteration
%
%   FLAG: Current state in which PlotFcn is called. Possible values are:
%           init: initialization state
%           iter: iteration state
%           done: final state
%
%   STOP: A boolean to stop the algorithm.
%
%   See also PATTERNSEARCH2, GA2, optimoptions2.


%   Copyright 2003-2015 The MathWorks, Inc.

stop = false;
switch flag
    case 'init'
        plotBest = plot(optimvalues.iteration,optimvalues.fval, '.b');
        set(plotBest,'Tag','psplotbestf2');
        xlabel('Iteration','interp','none'); 
        ylabel('Function value','interp','none')
        title(sprintf('Best Function Value: %g',optimvalues.fval),'interp','none');
    case 'iter'
        plotBest = findobj(get(gca,'Children'),'Tag','psplotbestf2');
        newX = [get(plotBest,'Xdata') optimvalues.iteration];
        newY = [get(plotBest,'Ydata') optimvalues.fval];
        set(plotBest,'Xdata',newX, 'Ydata',newY);
        set(get(gca,'Title'),'String',sprintf('Best Function Value: %g',optimvalues.fval));
end
