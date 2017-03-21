function parents = selectionuniform2(expectation,nParents,options)
%SELECTIONUNIFORM2 Choose parents at random.
%   PARENTS = SELECTIONUNIFORM2(EXPECTATION,NPARENTS,OPTIONS) chooses
%   PARENTS randomly using the EXPECTATION and number of parents NPARENTS. 
%
%   Parent selection is NOT a function of performance. This selection function 
%   is useful for debugging your own custom selection, or for comparison. It is 
%   not useful for actual evolution of high performing individuals. 
%
%   Example:
%   Create an options structure using SELECTIONUNIFORM2 as the selection
%   function.
%     options = optimoptions2('ga2', 'SelectionFcn', @selectionuniform2);
%
%   (Note: If calling gamultiobj2, replace 'ga2' with 'gamultiobj2') 

%   Copyright 2003-2015 The MathWorks, Inc.

expectation = expectation(:,1);
% nParents random numbers
parents = rand(1,nParents);

% integers on the interval [1, populationSize]
parents = ceil(parents * length(expectation));
