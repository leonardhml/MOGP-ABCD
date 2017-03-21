function [oFun, cFun] = getSetOptimFcns2(oFunIn, cFunIn)
% Hold the values of oFunIn and cFunIn as persistent storage. When these
% two input arguments are provided, the persistent storage is updated with
% these new values.
% Note: Caller must provide two input arguments to update. No error checks
% in this function because:
% - Private use by optimization solvers only.
% - This function should be efficient because it is used in a loop (finite
%   difference) 

%   Copyright 2014 The MathWorks, Inc.

persistent currentOFun currentCFun

if ~mislocked
  % Make sure this file is locked. This will prevent the persistent
  % variables from being reset if 'clear functions' is executed.
  mlock
end

% Return the existing values.
oFun = currentOFun;
cFun = currentCFun;

% Also, update the persistent storage if new values are provided.
if nargin == 2
    currentOFun = oFunIn;
    currentCFun = cFunIn;
end
