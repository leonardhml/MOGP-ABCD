function [stop,options,optchanged]  = psoutputfcntemplate2(optimvalues,options,flag)
%PSOUTPUTFCNTEMPLATE2 Template to write custom OutputFcn for PATTERNSEARCH2.
%   [STOP,OPTIONS,OPTCHANGED] = PSOUTPUTFCNTEMPLATE2(OPTIMVALUES,OPTIONS,FLAG) 
%   where OPTIMVALUES is a structure containing information about the state
%   of the optimization:
%            x: current point X 
%    iteration: iteration number
%         fval: function value 
%     meshsize: current mesh size 
%    funccount: number of function evaluations
%       method: method used in last iteration 
%       TolFun: change in fval from previous iteration
%         TolX: norm of change in X from previous iteration
%   nonlinineq: nonlinear inequality constraint values, when nonlinear constraints exist
%     nonlineq: nonlinear equality constraint values, when nonlinear constraints exist
%
%   OPTIONS: Options object used by PATTERNSEARCH2.
%
%   FLAG: Current state in which OutPutFcn is called. Possible values are:
%         'init': initialization state 
%         'iter': iteration state
%    'interrupt': subproblem for nonlinear constraints state
%         'done': final state
% 		
%   STOP: A boolean to stop the algorithm.
%
%   OPTCHANGED: A boolean indicating if the options have changed.
%
%	See also PATTERNSEARCH2, GA2, optimoptions2, SEARCHFCNTEMPLATE2

%   Copyright 2003-2015 The MathWorks, Inc.

stop = false;
optchanged = false;

switch flag
    case 'init'
        disp('Starting the algorithm');
    case {'iter','interrupt'}
        disp('Iterating ...')
    case 'done'
        disp('Performing final task');
end
  