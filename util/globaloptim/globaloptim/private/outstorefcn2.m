function stop = outstorefcn2(x, optimValues, state, osc)
%OUTSTOREFCN2 Output function to retrieve failed local solver call
%information.
%
%   STOP = OUTSTOREFCN2(X, OPTIMVALUES, STATE, OSC) is the output function
%   which is used to retrieve algorithm information from calls to the local
%   solver which error. OSC is an output store container, a structure of
%   nested function handles. This structure allows us to store and update
%   algorithm specific information from the local solver as the solver
%   progresses.
%
%   See also OUTPUTSTORECONTAINER2, GLOBALSEARCHNLP2, FMULTISTART2

%   Copyright 2009 The MathWorks, Inc.

% Get the output store from the container
os = osc.getOutputStore();

% Always record updates in the number of function calls as some algorithms
% call user functions between major iterations. For example the line search2
% phase in active-set needs to evaluate the objective and constraint
% functions to calculate its merit function. 
os.funcCount = optimValues.funccount;

% Want to record the remaining information for completed iterations only
if ~strcmp(state, 'interrupt')
    os.iterations = optimValues.iteration;
end

% Set the updated output store back
osc.setOutputStore(os);

% Don't want to terminate the algorithm in this function
stop = false;


