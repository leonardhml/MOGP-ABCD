function [stop,varargout] = psoutput2(OutputFcns, OutputFcnArgs,...
    optimval,optold,flag,varargin)
%PSOUTPUT2 Helper function that manages the output functions.
%
%   STATE = PSOUTPUT2(OPTIMVAL, OPTOLD, FLAG) runs each of
%   the output functions in the options.OutputFcn cell array. This API is
%   used when FLAG is 'interrupt' or 'done'.
%
%   [STATE, OPTIONS] = PSOUTPUT2(OPTIMVAL, OPTOLD, FLAG, OPTIONS,
%   DEFAULTOPT, NUMBEROFVARIABLES) runs each of the output functions in the
%   options.OutputFcn cell array. This API is used when FLAG is 'init' or
%   'iter'.
%
%   Private to PFMINLCON2, PFMINBND2, PFMINUNC2.

%   Copyright 2003-2015 The MathWorks, Inc.

%Initialize
if any(strcmpi(flag, {'init', 'iter'}))
    options = varargin{1};
    defaultopt = varargin{2};
    numberOfVariables = varargin{3};    
end
stop   = false;
optchanged = false;

% get the functions and return if there are none
if(isempty(OutputFcns))
    return
end
% call each output function
stop = false(length(OutputFcns),1);
for i = 1:length(OutputFcns)
    [stop(i) ,optnew , changed ] = feval(OutputFcns{i},optimval,optold,flag,OutputFcnArgs{i}{:});
    if changed  %If changes are not duplicates, we will get all the changes
        optold = optnew;
        optchanged = true;
    end
end

% Update options if user has changed them. We only have to check this in
% the 'iter' and 'init' cases.
if optchanged && any(strcmpi(flag, {'init', 'iter'}))
    optnew = prepareOptionsForSolver2(optold, 'patternsearch2');
    options = psoptimset2(options, optnew);
    options = checkoptions2(options,defaultopt,numberOfVariables);
    % Add TolFunValue back
    options.TolFunValue = options.TolFun;
    options.OutputPlotFcnOptions = copyForOutputAndPlotFcn(...
        optim.options.PatternsearchOptions2, options);
end

% If any stop(i) is true we set the stop to true
stop = any(stop);

% Return options if flag is 'init' or 'iter'
if any(strcmpi(flag, {'init', 'iter'}))
   varargout{1} = options; 
end
