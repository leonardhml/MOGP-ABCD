function stop = psplotbestx2(optimvalues,flag)
%PSPLOTBESTX2 PlotFcn to plot best X value.
%   STOP = PSPLOTBESTX2(OTIMVALUES,FLAG) where OPTIMVALUES is a structure
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
        set(gca,'xlimmode','manual','zlimmode','manual', ...
            'alimmode','manual')
        title('Current Best Point','interp','none')
        Xlength = numel(optimvalues.x);
        xlabel(sprintf('Number of variables (%i)',Xlength),'interp','none');
        ylabel('Current best point','interp','none');
        plotBestX = bar(optimvalues.x(:));
        set(plotBestX,'Tag','psplotbestx2');
        set(plotBestX,'edgecolor','none')
        set(gca,'xlim',[0,1 + Xlength])
    case 'iter'
        plotBestX = findobj(get(gca,'Children'),'Tag','psplotbestx2');
        set(plotBestX,'Ydata',optimvalues.x(:))
end





