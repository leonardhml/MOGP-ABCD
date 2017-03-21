function temperature = satemperaturefcntemplate2(optimValues,options)
%SATEMPERATUREFCNTEMPLATE2 Template to write custom temperature function
%   TEMPERATURE = SATEMPERATUREFCNTEMPLATE2(optimValues,options) updates the
%   current temperature
%
%   OPTIMVALUES is a structure containing the following information:
%              x: current point 
%           fval: function value at x
%          bestx: best point found so far
%       bestfval: function value at bestx
%    temperature: current temperature
%      iteration: current iteration 
%             t0: start time
%              k: annealing parameter 'k'
%
%   OPTIONS: options object created by using optimoptions2.

%   Copyright 2006-2015 The MathWorks, Inc.

% The body of this function should simply return a new temperature vector 
% The new temperature will generally depend on either
% options.InitialTemperature or optimValues.temperature and on
% optimValues.k 

% Make sure temperature is of same length as optimValues.k (nvars)
temperature = options.InitialTemperature./(optimValues.k).^2;
