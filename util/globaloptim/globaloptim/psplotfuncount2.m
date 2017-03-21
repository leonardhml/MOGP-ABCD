function stop = psplotfuncount2(optimvalues,flag)
%PSPLOTFUNCOUNT2 Plot function evaluations every iteration.
%   STOP = PSPLOTFUNCOUNT2(OPTIMVALUES,FLAG) where OPTIMVALUES is a
%   structure with the following fields: 
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
    plotFunC = plot(optimvalues.iteration, optimvalues.funccount,'kd', ...
        'MarkerSize',5,'MarkerFaceColor',[1 0 1]);
    set(plotFunC,'Tag','psplotfuncount2');
    xlabel('Iteration','interp','none'); 
    ylabel('Function evaluations per interval','interp','none');
    title(sprintf('Total Function Evaluations: %i',optimvalues.funccount),'interp','none');
    case 'iter'
    plotFunC = findobj(get(gca,'Children'),'Tag','psplotfuncount2');
    totalFunCount = sum(get(plotFunC,'Ydata'));
    newX = [get(plotFunC,'Xdata') optimvalues.iteration];
    newY = [get(plotFunC,'Ydata') (optimvalues.funccount - totalFunCount)];
    set(plotFunC,'Xdata',newX,'Ydata',newY);
    set(get(gca,'Title'),'String',sprintf('Total Function Evaluations: %i',optimvalues.funccount));
end

