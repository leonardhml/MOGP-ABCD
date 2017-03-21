function gahistplotupdate2(~, eventData)
%GAHISTPLOTUPDATE2 recomputes the x-axis tick labels
% 	based on the set x-axis limits to only show integer values.

%   Copyright 2015 The MathWorks, Inc.

% Get the requested x-axes limits
xlimits = get(eventData.AffectedObject,'xlim');
% Get the closest integers around the requested limits
xmin = floor(xlimits(1));
xmax = ceil(xlimits(end));
% Plot at most 10 labels (otherwise the plot gets too crowded)
nskip = max(1,floor((xmax-xmin+1)/10));
% Set the new integer Xticks
set(eventData.AffectedObject,'Xtick',xmin:nskip:xmax);