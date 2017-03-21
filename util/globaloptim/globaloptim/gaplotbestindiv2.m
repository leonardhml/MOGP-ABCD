function state = gaplotbestindiv2(~,state,flag)
%GAPLOTBESTINDIV2 Plots the best individual.
%   STATE = GAPLOTBESTINDIV2(OPTIONS,STATE,FLAG) plots the best 
%   individual's genome as a histogram, with the number of bins
%   in the histogram equal to the length of the genome.
%
%   Example:
%    Create an options structure that uses GAPLOTBESTINDIV2
%    as the plot function
%     options = optimoptions2('ga2','PlotFcn',@gaplotbestindiv2);

%   Copyright 2003-2015 The MathWorks, Inc.

if  size(state.Score,2) > 1
    title('Best Individual Plot: not available','interp','none');
    return;
end

switch flag
    case 'init'
        GenomeLength = size(state.Population,2);
        [~,i] = min(state.Score);
        h = bar(double(state.Population(i,:)));
        set(h,'edgecolor','none','Tag','gaplotbestindiv2')
        set(gca,'xlim',[0,1 + GenomeLength])
        title('Current Best Individual','interp','none')
        xlabel(sprintf('Number of variables (%i)',GenomeLength),'interp','none');
        ylabel('Current best individual','interp','none');

    case 'iter'
        [~,i] = min(state.Score);
        h = findobj(get(gca,'Children'),'Tag','gaplotbestindiv2');
        set(h,'Ydata',double(state.Population(i,:)));
end


