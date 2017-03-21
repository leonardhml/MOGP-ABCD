function state = gaplotrankhist2(~,state,flag)
%GAPLOTRANKHIST2 Plots a histogram of all ranks.
%
%   Example:
%    Create an options structure that will use GAPLOTRANK
%    as the plot function
%     options = optimoptions2(@gamultiobj2,'PlotFcn',@gaplotrank);

%   Copyright 2007-2015 The MathWorks, Inc.

if ~isfield(state,'Rank') 
      title('Rank Histogram: not available','interp','none');
    return;
end

% maximum limit for x axis
maxRank = 5;

switch flag
    case 'init'
        title('Rank histogram','interp','none')
        xlabel('Rank','interp','none')
        ylabel('Number of individuals','interp','none')
        addlistener(gca,'XLim','PostSet',@gahistplotupdate2);
        xlim([0 maxRank])
    case 'iter'
        allRank = state.Rank;
        h = histogram(gca, allRank,'BinMethod','integers');
        % update the xlim if necessary
        if maxRank <= h.NumBins
            set(gca,'xlim', [0, h.NumBins+1])
        end
end