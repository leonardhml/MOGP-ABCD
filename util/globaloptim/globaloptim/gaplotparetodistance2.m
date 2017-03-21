function state = gaplotparetodistance2(options,state,flag)
%GAPLOTPARETODISTANCE2 Plots distance value of each individuals.
%
%   Example:
%    Create an options structure that will use GAPLOTPARETODISTANCE2
%    as the plot function
%     options = optimoptions2('ga2','PlotFcn',@gaplotparetodistance2);

%   Copyright 2007-2015 The MathWorks, Inc.

if ~isfield(state,'Distance') 
      title('Pareto Distance Plot: not available','interp','none');
    return;
end

colors = ['g','m','b','k','r','c'];
subPops = populationIndices2(options.PopulationSize);
[~,c] = size(subPops);

switch flag
    case 'init'
         set(gca,'xlimmode','manual','zlimmode','manual', ...
            'climmode','manual','alimmode','manual')
        % make the x axis span exactly the data size
        set(gca,'xlim',[0,1 + sum(options.PopulationSize)])
        title('Distance of individuals','interp','none')
        xlabel('Individuals');
        ylabel('Distance');
        for i = 1:c
            pop = subPops(:,i);
            range = pop(1):pop(2);
            h = bar(range,state.Distance(range),colors(1 + mod(i,5)));
            str = ['gaplotparetodistance2',int2str(i)];
            set(h,'edgecolor','none','facecolor',colors(1 + mod(i,5)),'Tag',str)
            hold on;
        end
        hold off;

    case 'iter'
       for i = 1:c
            pop = subPops(:,i);
            range = pop(1):pop(2);
            str = ['gaplotparetodistance2',int2str(i)];
            h = findobj(get(gca,'Children'),'Tag',str);
            set(h,'Ydata',state.Distance(range));
        end
end
