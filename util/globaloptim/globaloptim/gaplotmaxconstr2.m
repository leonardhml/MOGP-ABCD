function state = gaplotmaxconstr2(options,state,flag)
%GAPLOTMAXCONSTR2 Plots the maximum nonlinear constraint violation by GA2.
%   STATE = GAPLOTMAXCONSTR2(OPTIONS,STATE,FLAG) plots the maximum nonlinear
%   constraint violation
%
%   Example:
%    Create an options structure that will use GAPLOTMAXCONSTR2
%    as the plot function
%     options = optimoptions2('ga2','PlotFcn',@gaplotmaxconstr2);
%
%    (Note: If calling gamultiobj2, replace 'ga2' with 'gamultiobj2') 

%   Copyright 2005-2015 The MathWorks, Inc.

% If the flag is 'interrupt' simply return
if strcmpi(flag,'interrupt')
    return;
end
if ~(isfield(state,'NonlinIneq') || isfield(state,'NonlinEq'))
    title('zero nonlinear constraint'); 
    return;
end
% Maximum constraint violation
maxConstr = 0;
if ~isempty(state.NonlinIneq)
  maxConstr = max([0;state.NonlinIneq(:)]);
end
if ~isempty(state.NonlinEq)
    maxConstr = max([maxConstr; abs(state.NonlinEq(:))]);
end

switch flag
    case 'init'
        hold on;
        set(gca,'xlim',[0,options.MaxGenerations]);
        xlabel('Generation','interp','none');
        ylabel('Max constraint','interp','none');
        plotConstr = plot(state.Generation,maxConstr,'.k');
        set(plotConstr,'Tag','gaplotconstr');
         title(sprintf('Max constraint: %g',maxConstr),'interp','none');
    case 'iter'
        plotConstr = findobj(get(gca,'Children'),'Tag','gaplotconstr');
        newX = [get(plotConstr,'Xdata') state.Generation];
        newY = [get(plotConstr,'Ydata') maxConstr];
        set(plotConstr,'Xdata',newX, 'Ydata',newY);
         title(sprintf('Max constraint: %g',maxConstr),'interp','none');
    case 'done'
        hold off;
end
