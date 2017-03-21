function [state, options,optchanged] = gaoutputfcntemplate2(options,state,flag)
%GAOUTPUTFCNTEMPLATE2 Template to write custom OutputFcn for GA2.
%   [STATE, OPTIONS, OPTCHANGED] = GAOUTPUTFCNTEMPLATE2(OPTIONS,STATE,FLAG)
%   where OPTIONS is an options structure used by GA2. 
%
%   STATE: A structure containing the following information about the state 
%   of the optimization:
%             Population: Population in the current generation
%                  Score: Scores of the current population
%             Generation: Current generation number
%              StartTime: Time when GA2 started 
%               StopFlag: String containing the reason for stopping
%              Selection: Indices of individuals selected for elite,
%                         crossover and mutation
%            Expectation: Expectation for selection of individuals
%                   Best: Vector containing the best score in each generation
%        LastImprovement: Generation at which the last improvement in
%                         fitness value occurred
%    LastImprovementTime: Time at which last improvement occurred
%
%   FLAG: Current state in which OutputFcn is called. Possible values are:
%         init: initialization state 
%         iter: iteration state
%    interrupt: intermediate state
%         done: final state
% 		
%   STATE: Structure containing information about the state of the
%          optimization.
%
%   OPTCHANGED: Boolean indicating if the options have changed.
%
%	See also PATTERNSEARCH2, GA2, optimoptions2

%   Copyright 2004-2015 The MathWorks, Inc.


optchanged = false;

switch flag
 case 'init'
        disp('Starting the algorithm');
    case {'iter','interrupt'}
        disp('Iterating ...')
    case 'done'
        disp('Performing final task');
end
