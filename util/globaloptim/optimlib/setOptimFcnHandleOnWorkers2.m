function cleanupObj = setOptimFcnHandleOnWorkers2(useparallel,funfcn,confcn)

% Used for parallel finite difference by fmincon2, fminimax, and fgoalattain
% solvers.
% This function sends objective and constraint functions to workers which
% later are used by parfinitedifference function. parfinitedifference is
% used by fmincon2, fminimax, and fgoalattain solvers. These solvers must
% call this function to setup the workers.
%
% The idea is that we want to avoid sending function handles in every
% optimization iteration. This helps when function handle is large.
%
% This function sends to the workers objective and constraint function
% handles. These handles are stored as persistent on the workers.
% parfinitedifference function gets the function handles from the workers. 
%
% This function also creates onCleanup object and returns to the caller.
% This will make sure that when the fmincon2, fminimax, or fgoalattain
% solvers are done, the persistent storage on the workers will be removed. 
%
% If there are any recursive calls to any of the three solvers, the
% persistent storage is limited to the current stack entry. The older
% values are stored as part of cleanup object. When the stack unwinds, the
% old values are restored.
%
% Inputs:
% useparallel: must be true or false indicates if parallel finite
%              difference is to be used or not
% funfcn: (Objective function information) A cell array used by
%         optimization nonlinear solvers. The first element indicates if
%         finite difference is to be used or not. The third element
%         contains the function handle.
% confcn: Similar syntax as funfcn but for constraint function.
%
% Output:
% cleanupObj: The caller must create a variable as the output of this
%             function. This is a onCleanup object (see help for this
%             class) that will ensure that the cleanup work on the workers
%             is completed

%   Copyright 2014 The MathWorks, Inc.


if ~validateopts_UseParallel(useparallel,true,true)
    % When user does not set UseParallel option to true, return now.
    cleanupObj = onCleanup.empty;
    return;
end

if strcmp(funfcn{1},'fun')
    % Solvers will use finite difference for objective; get the function
    % handle.
    objfun = funfcn{3};
else
    % No finite difference; set the function handle to empty.
    objfun = '';
end

if strcmp(confcn{1},'fun')
    % Solvers will use finite difference for objective; get the function
    % handle.
    confun = confcn{3};
else
    % No finite difference; set the function handle to empty.
    confun = '';
end

% Check for both objfun and confun empty: return immediately if no finite
% difference is needed (irrespective of useparallel value).
if isempty(objfun) && isempty(confun)
    cleanupObj = onCleanup.empty;
    return;
end

% At this point, we know that solvers will execute parfor body for finite
% difference. We need to access the workers now in order to send function
% handles. It will start the pool if user preference is set.
try
    % Get a handle to the pool.
    pool = gcp;
catch
    % PCT not installed or pool creation threw an error
    pool = [];
end

if ~isempty(pool)
    % Get the current value of the persistent storage. Unless we have a
    % plan to deal with failures on the workers, do not use parfeval. Use
    % parfor and halt if there is any error.
    old_objfun = cell(1,pool.NumWorkers);
    old_confun = cell(1,pool.NumWorkers);
    
    parfor ii = 1:pool.NumWorkers
      % Get existing values and set new values
      [old_objfun{ii}, old_confun{ii}] = getSetOptimFcns2(objfun, confun)
    end
    % Setup a cleanup task to restore the persistent storage on the workers.
    cleanupObj = onCleanup(@() CleanupWorkers(pool, old_objfun, old_confun));
else
    % The pool is not available to use. The solver is still going to use the
    % parfor loop (because finite difference is 'on' and UseParallel is true)
    % This means, we must set the persistent storage on the local host.        
    [old_objfun, old_confun] = getSetOptimFcns2(objfun, confun);
    cleanupObj = onCleanup(@() CleanupLocal(old_objfun, old_confun));
end

function CleanupWorkers(pool, old_objfun, old_confun)
% Cleanup function to restore persistent storage (in getSetOptimFcns2 file
% on workers) to old values on the workers.  

parfor ii = 1:pool.NumWorkers
  % Restore the persistent storage (could be [] and it's ok)
  getSetOptimFcns2(old_objfun{ii}, old_confun{ii});
end


function CleanupLocal(old_objfun, old_confun)
% Cleanup function to restore persistent storage (in getSetOptimFcns2 file
% on host) to old values.

% Restore the persistent storage (could be [] and it's ok)
getSetOptimFcns2(old_objfun, old_confun);


