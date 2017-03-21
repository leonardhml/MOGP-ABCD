function state = gaplotrange2(options,state,flag)
%GAPLOTRANGE2 Plots the mean and the range of the scores.
%   STATE = GAPLOTRANGE2(OPTIONS,STATE,FLAG) plots the mean and the range
%   (best and the worst) of scores.  
%
%   Example:
%   Create an options structure that uses GAPLOTRANGE2
%   as the plot function
%     options = optimoptions2('ga2','PlotFcn',@gaplotrange2);
%
%    (Note: If calling gamultiobj2, replace 'ga2' with 'gamultiobj2') 

%   Copyright 2003-2015 The MathWorks, Inc.

if isinf(options.MaxGenerations) || size(state.Score,2) > 1
    title('Plot Not Available','interp','none');
    return;
end
generation = state.Generation;
score = state.Score;
smean = mean(score);
Y = smean;
L = smean - min(score);
U = max(score) - smean;

switch flag

    case 'init'
        set(gca,'xlim',[1,options.MaxGenerations+1]);
        plotRange = errorbar(generation,Y,L,U);
        set(plotRange,'Tag','gaplotrange_errorbar');
        title('Best, Worst, and Mean Scores','interp','none')
        xlabel('Generation','interp','none')
    case 'iter'
        plotRange = findobj(get(gca,'Children'),'Tag','gaplotrange_errorbar');
        
        oldX = get(plotRange,'Xdata');
        newX = [oldX(:);generation];
        
        oldY = get(plotRange,'Ydata');
        newY = [oldY(:); Y];
        
        oldL = get(plotRange,'Ldata');
        newL = [oldL(:);L];
        
        oldU = get(plotRange,'Udata');
        newU = [oldU(:);U];
        
        set(plotRange, 'Xdata',newX,'Ydata',newY,'Ldata',newL,'Udata',newU);
end

