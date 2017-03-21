function temperature = temperaturefast2(optimValues,options)
%TEMPERATUREFAST2 Updates the temperature vector for annealing process 
%   TEMPERATURE = TEMPERATUREFAST2(optimValues,options) uses fast 
%   annealing by updating the current temperature based on the initial
%   temperature and the current annealing parameter k
%
%   OPTIMVALUES is a structure containing the following information:
%              x: current point 
%           fval: function value at x
%          bestx: best point found so far
%       bestfval: function value at bestx
%    temperature: current temperature
%      iteration: current iteration 
%             t0: start time
%              k: annealing parameter
%
%   OPTIONS: options object created by using optimoptions2.
%
%   Example:
%    Create an options structure using TEMPERATUREFAST2 as the annealing
%    function
%    options = optimoptions2('simulannealbnd2','TemperatureFcn',@temperaturefast2);

%   Copyright 2006-2015 The MathWorks, Inc.

temperature = options.InitialTemperature./optimValues.k;

