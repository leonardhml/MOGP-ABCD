function state = gaplotdistance2(options,state,flag)
%GAPLOTDISTANCE2 Averages several samples of distances between individuals.
%   STATE = GAPLOTDISTANCE2(OPTIONS,STATE,FLAG) plots an averaged distance
%   between individuals.
%
%   Example:
%    Create an options structure that uses GAPLOTDISTANCE2
%    as the plot function
%     options = optimoptions2('ga2','PlotFcn',@gaplotdistance2);
%
%    (Note: If calling gamultiobj2, replace 'ga2' with 'gamultiobj2') 

%   Copyright 2003-2015 The MathWorks, Inc.

samples = 20;
choices = ceil(sum(options.PopulationSize) * rand(samples,2));
switch flag
    case 'init'
        population = state.Population;
        distance = 0;
        for i = 1:samples
            d = population(choices(i,1),:) - population(choices(i,2),:);
            distance = distance + sqrt( sum ( d.* d));
        end
        plotDist = plot(state.Generation,distance/samples,'.');
        set(gca,'xlimmode','manual','zlimmode','manual', ...
            'alimmode','manual')
        set(gca,'xlim',[1,options.MaxGenerations]);
        set(plotDist,'Tag','gaplotdistance2');
        xlabel('Generation','interp','none');
        ylabel('Avergae Distance');
        title('Average Distance Between Individuals','interp','none')

    case 'iter'
        population = state.Population;
        distance = 0;
        for i = 1:samples
            d = population(choices(i,1),:) - population(choices(i,2),:);
            distance = distance + sqrt( sum ( d.* d));
        end
        plotDist = findobj(get(gca,'Children'),'Tag','gaplotdistance2');
        newX = [get(plotDist,'Xdata') state.Generation];
        newY = [get(plotDist,'Ydata') distance/samples];
        set(plotDist,'Xdata',newX,'Ydata',newY);
end
