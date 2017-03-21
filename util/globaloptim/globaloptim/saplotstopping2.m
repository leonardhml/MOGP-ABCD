function stop = saplotstopping2(options,optimvalues,flag)
%SAPLOTSTOPPING2 PlotFcn to plot stopping criteria satisfaction.
%   STOP = SAPLOTSTOPPING2(OPTIONS,OPTIMVALUES,FLAG) where OPTIMVALUES is a
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
%    Create an options structure that will use SAPLOTSTOPPING2
%    as the plot function
%     options = optimoptions2('simulannealbnd2','PlotFcn',@saplotstopping2);

%   Copyright 2006-2015 The MathWorks, Inc.

stop = false;
% Calculate fraction of 'doneness' for each criterion
func = optimvalues.funccount / options.MaxFunctionEvaluations;
iter = optimvalues.iteration / options.MaxIterations;
time = (cputime-optimvalues.t0) / options.MaxTime;

% Multiply ratios by 100 to get percentages
ydata = 100 * [time, iter, func];

switch flag
    case 'init'
        barh(ydata,'Tag','saplotstopping2')
        set(gca,'xlim',[0,100],'yticklabel', ...
            {'Time','Iteration', 'f-count'},'climmode','manual')
        xlabel('% of criteria met','interp','none')
        title('Stopping Criteria','interp','none')
    case 'iter'
        ch = findobj(get(gca,'Children'),'Tag','saplotstopping2');
        set(ch,'YData',ydata);
end