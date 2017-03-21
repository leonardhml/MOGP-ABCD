function stop = psplotmeshsize2(optimvalues,flag)
%PSPLOTMESHSIZE2 PlotFcn to plot mesh size.
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
        plotMesh = plot(optimvalues.iteration,optimvalues.meshsize, 'm');
        set(plotMesh,'Tag','psplotmeshsize2');
        xlabel('Iteration','interp','none'); 
        ylabel('Mesh size','interp','none');
        title(sprintf('Current Mesh Size: %g',optimvalues.meshsize),'interp','none')
    case 'iter'
        plotMesh = findobj(get(gca,'Children'),'Tag','psplotmeshsize2');
        newX = [get(plotMesh,'Xdata') optimvalues.iteration];
        newY = [get(plotMesh,'Ydata') optimvalues.meshsize];
        set(plotMesh,'Xdata',newX, 'Ydata',newY);
        set(get(gca,'Title'),'String',sprintf('Current Mesh Size: %g',optimvalues.meshsize)); 
end

