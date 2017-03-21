function state = gaplotexpectation2(~,state,flag)
%GAPLOTEXPECTATION2 Plots raw scores against the expected number of offspring.
%   STATE = GAPLOTEXPECTATION2(OPTIONS,STATE,FLAG) plots the raw scores
%   against the expected number of offspring.
%
%   Example:
%    Create an options structure that uses GAPLOTEXPECTATION2
%    as the plot function
%     options = optimoptions2('ga2','PlotFcn',@gaplotexpectation2);
%
%    (Note: If calling gamultiobj2, replace 'ga2' with 'gamultiobj2') 

%   Copyright 2003-2015 The MathWorks, Inc.

% we have to store scores because the expectation in the state are for the
% last generation and the scores are for the next generation.
persistent scores;

if  ~isfield(state,'Expectation') || size(state.Score,2) > 1
    title('Expectation Plot: not available','interp','none');
    return;
end

switch flag
    case 'init'
        scores = state.Score;
        % Expectation is a vector of zeros because there is no expectation
        % in first iteration at all.
        plotExpec = plot(scores,zeros(1,length(scores)),'.');
        set(plotExpec,'Tag','gaplotexpectation2');
        xlabel('Raw scores','interp','none');
        ylabel('Expectation','interp','none');
        title('Fitness Scaling','interp','none')
    case 'iter'
        %This is a safeguard when population size is reduced at run time.
        if length(state.Score) ~= length(state.Expectation)
            scores = state.Score(length(state.Score)-length(state.Expectation)+1:end);
        end
        plotExpec = findobj(get(gca,'Children'),'Tag','gaplotexpectation2');
        set(plotExpec,'Xdata',scores,'Ydata',state.Expectation);
        scores = state.Score;
end
