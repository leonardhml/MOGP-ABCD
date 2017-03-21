function state = gaplotselection2(options,state,flag)
%GAPLOTSELECTION2 A histogram of parents.
%   STATE = GAPLOTSELECTION2(OPTIONS,STATE,FLAG) plots a histogram of the
%   parents. This will let you see which parents are contributing to each 
%   generation. If there is not enough spread (only a few parents are being 
%   used) then you may want to change some parameters to get more diversity.
%
%   Example:
%    Create an options structure that uses GAPLOTSELECTION2
%    as the plot function
%      options = optimoptions2('ga2','PlotFcn',@gaplotselection2);
%
%    (Note: If calling gamultiobj2, replace 'ga2' with 'gamultiobj2') 

%   Copyright 2003-2015 The MathWorks, Inc. 

switch flag
    case 'init'
        title('Selection Function','interp','none')
        xlabel('Individual','interp','none')
        ylabel('Number of children','interp','none')
        addlistener(gca,'XLim','PostSet',@gahistplotupdate2);
        % Set the default xlim (can be changed interactively by user)
        xlim([0 options.PopulationSize(1)+1])
    case 'iter'
        h = histogram(gca,state.Selection,'BinMethod','integers');
        xLimits = get(gca,'xlim');
        % update the xlim if necessary
        if xLimits(end) < h.BinEdges(end)
            xlim([xLimits(1), h.BinEdges(end)+0.5]);
        end
    case 'done'
        % Add tag to the axis so that it can be easily found.
        set(gca,'Tag','gaplotselection2');
end

