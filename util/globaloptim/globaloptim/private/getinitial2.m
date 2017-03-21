function [optimState,nextIterate,MeshSize,EXITFLAG,run] = ...
    getinitial2(Iterate,numberOfVariables,neqcstr,lb,ub,options)
%GETINITIAL2 is private to pfminlcon2, pfminbnd2 and pfminunc2.

%   Copyright 2003-2006 The MathWorks, Inc.

% Initialization
optimState.Iter = 0;
optimState.FunEval = 1;
optimState.infMessage  = '';
optimState.stopPlot = false;
optimState.stopOutput = false;
optimState.deltaF = NaN;
optimState.deltaX = NaN;
optimState.Successdir = [];
optimState.how = ' ';
optimState.MeshCont = options.MeshContraction;
optimState.scale = ones(numberOfVariables,1);
EXITFLAG = -1;
run = true;
MeshSize = options.InitialMeshSize;

% Calculate scale
if any(strcmpi(options.ScaleMesh,{'dynamic','on'})) && ~neqcstr
    meanX = mean([Iterate.x],2);
    optimState.scale = logscale2(lb,ub,meanX);
end

nextIterate = Iterate;
optimState.StartTime = cputime;
