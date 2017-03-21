function parents = selectionroulette2(expectation,nParents,options)
%SELECTIONROULETTE2 Choose parents using roulette wheel.
%   PARENTS = SELECTIONROULETTE2(EXPECTATION,NPARENTS,OPTIONS) chooses
%   PARENTS using EXPECTATION and number of parents NPARENTS. On each 
%   of the NPARENTS trials, every parent has a probability of being selected
%   that is proportional to their expectation.
%
%   Example:
%   Create an options structure using SELECTIONROULETTE2 as the selection
%   function
%     options = optimoptions2('ga2','SelectionFcn',@selectionroulette2);
%
%   (Note: If calling gamultiobj2, replace 'ga2' with 'gamultiobj2') 

%   Copyright 2003-2015 The MathWorks, Inc.

expectation = expectation(:,1);
wheel = cumsum(expectation) / nParents;

parents = zeros(1,nParents);
for i = 1:nParents
    r = rand;
    for j = 1:length(wheel)
        if(r < wheel(j))
            parents(i) = j;
            break;
        end
    end
end
    
