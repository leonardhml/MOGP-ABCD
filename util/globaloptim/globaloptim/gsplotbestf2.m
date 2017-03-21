function stop = gsplotbestf2(optimValues, state)
%GSPLOTBESTF2 Plot best function value.
%
%   STOP = GSPLOTBESTF2(OPTIMVALUES, STATE) plots OPTIMVALUES.BESTFVAL
%   against OPTIMVALUES.LOCALRUNINDEX. This function is called from
%   the global solver with the following inputs:
%
%   OPTIMVALUES: Information after the current local solver call.
%    localsolution: solution returned from current call to the local solver
%    localrunindex: index of current call to the local solver
%        funccount: number of function evaluations
%            bestx: best solution found so far               
%         bestfval: function value at bestx
%
%   OPTIMVALUES.LOCALSOLUTION: Solution returned from current call to the
%                              local solver.
%                X: solution returned from local solver
%             Fval: function value at X
%         Exitflag: exit flag from current call to the local solver
% 
%   STATE: Current state in which plot function is called. 
%          Possible values are:
%             init: initialization state 
%             iter: after every call to the local solver 
%             done: final state
%
%   STOP: A boolean to stop the algorithm.
%
%   See also GLOBALSEARCH2, MULTISTART2

%   Copyright 2010 The MathWorks, Inc.

% Initialize stop boolean to false.
stop = false;

% Update the plot
switch state
    case 'init'
        % Create the plot. As the local solver has not yet been called
        % (optimValues.localrunindex = 0), there is no best function value
        % available (optimValues.bestfval = []). We manually set
        % localrunindex to zero and best function value to NaN in this case.
        plot(0, NaN, 'b.', 'tag', 'gsplotbestf2');
        xlabel('Local solver call', 'interp', 'none');
        ylabel('Function value', 'interp', 'none');
        
        % Create a title for the plot.
        title('Best Function Value');

        % Create a text annotation to display a message when there is no
        % best function value.
        text(0.275, 0.5, 'Best function value is not available', ...
            'Units', 'normalized', 'Tag', 'gstextbestf');
        
    case 'iter'
        % It is possible for optimValues.bestfval to be empty.  
        if isempty(optimValues.bestfval)
            % Set bestfval to NaN so that we can plot it for this call.
            % Note that this point will not be visible in the plot.
            bestfval = NaN;
        else
            % We have a best function value. Remove the message stating
            % that optimValues.bestfval is empty and update the title.
            bestfval = optimValues.bestfval;
            textBest = findobj(get(gca,'Children'), 'Tag', 'gstextbestf');
            set(textBest, 'String', '');
            set(get(gca,'Title'), ...
                'String', sprintf('Best Function Value: %g',bestfval));
        end
                
        % Update the plot with the current best function value.
        plotBest = findobj(get(gca,'Children'), 'Tag', 'gsplotbestf2');
        newX = [get(plotBest,'Xdata') optimValues.localrunindex];
        newY = [get(plotBest,'Ydata') bestfval];       
        set(plotBest, 'Xdata', newX, 'Ydata', newY);
               
    case 'done'
        % No clean up tasks required for this plot function.
        
end