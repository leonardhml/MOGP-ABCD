function [stop,options,optchanged] = saoutputfcntemplate2(options,optimvalues,flag)
%SAOUTPUTFCNTEMPLATE2 Template to write custom OutputFcn for SIMULANNEALBND2
%   [STOP,OPTIONS,OPTCHANGED] = SAOUTPUTFCNTEMPLATE2(OPTIONS,OPTIMVALUES,FLAG) 
%   OPTIMVALUES is a structure with the following fields:  
%              x: current point 
%           fval: function value at x
%          bestx: best point found so far
%       bestfval: function value at bestx
%    temperature: current temperature
%      iteration: current iteration
%      funccount: number of function evaluations
%             t0: start time
%              k: annealing parameter 'k'
%
%   OPTIONS: The options object created by using optimoptions2
%
%   FLAG: Current state in which OutputFcn is called. Possible values are:
%           init: initialization state
%           iter: iteration state
%           done: final state
%
%   STOP: A boolean to stop the algorithm.
% 		
%   OPTCHANGED: A boolean indicating if the options have changed.
%

%   Copyright 2006-2015 The MathWorks, Inc.

stop = false;
optchanged = false;
switch flag
   case 'init'
        disp('Initializing output function');
    case 'iter'
        disp('Iterating ...')
    case 'done'
        disp('Performing final task');
end
