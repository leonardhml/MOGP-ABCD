function state = gaplotbestfun2(options,state,flag)
%GAPLOTBESTFUN2 Plots the best score and the mean score.

%   Copyright 2005-2015 The MathWorks, Inc. 

if(strcmp(flag,'init'))
    set(gca,'xlim',[1,options.MaxGenerations]);
    xlabel('Generation','interp','none');
    grid on
    ylabel('Fitness value','interp','none');
end

hold on;
generation = state.Generation;
best = min(state.Score);
plot(generation,best, 'v');
title(['Best: ',num2str(best)],'interp','none')
hold off;
        
