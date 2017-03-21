function [stop,options,optchanged]  = psoutputfile2(optimvalues,options,flag,fileName,interval)
%PSOUTPUTFILE2 Output function writes PATTERNSEARCH2 iterative history to a file.
%   [STOP, OPTIONS, OPTCHANGED] = PSOUTPUTFILE2(OPTIMVALUES,OPTIONS,FLAG, ...
%   FILE,INTERVAL) where OPTIMVALUES is a structure containing information
%   about the state of the optimization:
%            x: current point X
%    iteration: iteration number
%         fval: function value
%     meshsize: current mesh size
%    funccount: number of function evaluations
%       method: method used in last iteration
%       TolFun: tolerance on fval
%         TolX: tolerance on X
%     nonlinineq: nonlinear inequality constraints at X
%       nonlineq: nonlinear equality constraints at X
%
%   OPTIONS: Options object used by PATTERNSEARCH2.
%
%   FLAG: Current state in which OutPutFcn is called. Possible values are:
%         init: initialization state
%         iter: iteration state
%         done: final state
%
%   FILENAME: A file name where iterative output is written.
%   INTERVAL: Every 'INTERVAL' iterations are written; default = 1.
%
%   STOP: A boolean to stop the algorithm.
%   OPTCHANGED: A boolean indicating if the options have changed.
%
%   See also PSOUTPUTFCNTEMPLATE2, PATTERNSEARCH2, optimoptions2.

%   Copyright 2005-2015 The MathWorks, Inc.

stop = false;
optchanged = false;
mode = 'a';
if nargin < 5
    interval = 1;
end
if interval <= 0
    interval = 1;
end

if (rem(optimvalues.iteration,interval) ~=0)
    return;
end
type = optimvalues.problemtype;
% Maximum constraint
maxConstr = 0;
if strcmpi(type,'nonlinearconstr') && ~isempty(optimvalues.nonlinineq)
    maxConstr = max([maxConstr; optimvalues.nonlinineq(:)]);
end
if strcmpi(type,'nonlinearconstr') && ~isempty(optimvalues.nonlineq)
    maxConstr = max([maxConstr; abs(optimvalues.nonlineq(:))]);
end
[fid,theMessage] = fopen(fileName,mode);
if fid==-1
    error(message('globaloptim:psoutputfile2:fileWriteError', fileName, theMessage))
end

switch flag
    case 'init'
        % make sure that the file name is valid.
        msg = sprintf('\nIterative output generated by the PATTERNSEARCH2 solver on %s time', ...
            datestr(datenum(now)));
        fprintf(fid,msg);
        if ~strcmpi(type,'nonlinearconstr')
            fprintf(fid,'\n\nIter     f-count          f(x)      MeshSize\n');
        else
            fprintf(fid,'\n                                  max\n');
            fprintf(fid,'Iter   f-count      f(x)      constraint   MeshSize\n');
        end
    case 'iter'
        Iter = optimvalues.iteration;
        FunEval = optimvalues.funccount;
        MeshSize = optimvalues.meshsize;
        % Write to the file
        if ~strcmpi(type,'nonlinearconstr')
            fprintf(fid,'%5.0f    %5.0f   %12.6g  %12.4g\n',Iter, FunEval, ...
                optimvalues.fval, MeshSize);
        else
            fprintf(fid,'%3.0f   %4.0f  %12.6g  %12.4g %12.4g', ...
                Iter,FunEval,optimvalues.fval,maxConstr,MeshSize);
            fprintf(fid,'\n');
        end
    case 'done'
        Iter = optimvalues.iteration;
        FunEval = optimvalues.funccount;
        MeshSize = optimvalues.meshsize;
        % Write to the file
        if ~strcmpi(type,'nonlinearconstr')
            fprintf(fid,'%5.0f    %5.0f   %12.6g  %12.4g\n',Iter, FunEval, ...
                optimvalues.fval, MeshSize);
        else
            fprintf(fid,'%3.0f   %4.0f  %12.6g  %12.4g %12.4g', ...
                Iter,FunEval,optimvalues.fval,maxConstr,MeshSize);
            fprintf(fid,'\n');
        end
        fprintf(fid,'\nOptimization terminated.\n');
        fprintf(fid,'Best function value : %g\n',optimvalues.fval);
        fprintf(fid,'Best X found by PATTERNSEARCH2 solver\n[');
        fprintf(fid,' %g \n',optimvalues.x(:));
        fprintf(fid,']');
end
st = fclose(fid);
if st ~= 0
    error(message('globaloptim:psoutputfile2:fileCloseError', fileName))
end