function expectation = fitscalingtop2(scores,nParents,quantity)
%FITSCALINGTOP2 Top individuals reproduce equally (single objective only).
%   EXPECTATION = FITSCALINGTOP2(SCORES,NPARENTS,QUANTITY) calculates the
%   EXPECTATION using the SCORES and number of parents NPARENTS as well as 
%   QUANTITY.  QUANTITY represents either the number of expectations or 
%   the number of expectations in terms of the population size. If QUANTITY 
%   is not specified a default value of 0.4 is used. This function chooses the  
%   best QUANTITY scores for parents. QUANTITY can be an integer or a fraction
%   of the populationSize.
%
%   Example:
%   Create an options structure that uses FITSCALINGTOP2 as the fitness
%   scaling function and 0.5 as the QUANTITY
%     quantity = 0.5;
%     options = optimoptions2('ga2','FitnessScalingFcn',{@fitscalingtop2, quantity}); 

% 	 Each of the best n individuals have an equal chance of reproducing.
% 	 The rest have zero expectation. Expectation looks like:
% 	 [ 0 1/n 1/n 0 0 1/n 0 0 1/n ...]
% 	 this is called "Top fitness scaling" or "Truncation fitness scaling."
% 	 If quantity is set to 1, this is called "Best fitness scaling."

%   Copyright 2003-2015 The MathWorks, Inc.

if nargin < 3 || isempty(quantity)
    quantity = 0.4;
end
scores = scores(:);

% if less than one, it's a fraction, not the absolute number.
if quantity < 1
    quantity = round(quantity * length(scores));
end

expectation = zeros(size(scores));
% find the best ones
[~,i] = sort(scores);
expectation(i(1:quantity)) = nParents / quantity;
