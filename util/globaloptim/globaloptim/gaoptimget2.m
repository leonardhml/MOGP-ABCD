function o = gaoptimget2(options,name,default,flag)
%GAOPTIMGET2 Get GA2 OPTIONS parameter value.
%   VAL = GAOPTIMGET2(OPTIONS,'NAME') extracts the value of the named parameter
%   from optimization options structure OPTIONS, returning an empty matrix if
%   the parameter value is not specified in OPTIONS.  It is sufficient to
%   type only the leading characters that uniquely identify the
%   parameter.  Case is ignored for parameter names.  [] is a valid OPTIONS
%   argument.
%   
%   VAL = GAOPTIMGET2(OPTIONS,'NAME') extracts the parameter NAME 
%   For example:
%     
%     opts = gaoptimset2('Generations',200);
%     val = gaoptimget2(opts,'Generations');
%   
%   returns val = 200.
%   
%   See also GAOPTIMSET2.

%   Copyright 2003-2015 The MathWorks, Inc.

if nargin < 2
  error(message('globaloptim:gaoptimget2:inputarg'));
end
if nargin < 3
  default = [];
end
if nargin < 4
   flag = [];
end

% undocumented usage for fast access with no error checking
if isequal('fast',flag)
   o = gaoptimgetfast(options,name,default);
   return
end

if ~isempty(options) && ~isa(options,'struct')
  error(message('globaloptim:gaoptimget2:firstargerror'));
end

if isempty(options)
  o = default;
  return;
end

optionsstruct = struct('PopulationType', [], ...
    'PopInitRange',[], ...
    'PopulationSize',[], ...
    'EliteCount',[], ...
    'CrossoverFraction',[], ...
    'ParetoFraction',[], ...    
    'MigrationDirection',[], ...
    'MigrationInterval',[], ...
    'MigrationFraction',[], ...
    'Generations',[], ...
    'TimeLimit',[], ...
    'FitnessLimit',[], ...
    'StallGenLimit',[], ...
    'StallTimeLimit',[], ...
    'StallTest',[], ...
    'TolFun', [], ...
    'TolCon',[], ...
    'InitialPopulation',[], ...
    'InitialScores',[], ...
    'NonlinConAlgorithm',[], ...
    'InitialPenalty', [], ...
    'PenaltyFactor', [], ...
    'PlotInterval',[], ...
    'FitnessScalingFcn',[], ...
    'SelectionFcn',[], ...
    'CrossoverFcn',[], ...
    'MutationFcn',[], ...
    'DistanceMeasureFcn',[], ...
    'Display',[], ...
    'PlotFcns',[], ...
    'OutputFcns',[], ...
    'CreationFcn',[], ...
    'HybridFcn',[], ...
    'Vectorized',[], ...
    'UseParallel',[]);      

Names = fieldnames(optionsstruct);
%[m,n] = size(Names);
names = lower(Names);

lowName = lower(name);
j = strmatch(lowName,names);
if isempty(j) % if no matches
  error(message('globaloptim:gaoptimget2:invalidproperty',name));
elseif length(j) > 1            % if more than one match
  % Check for any exact matches (in case any names are subsets of others)
  k = strmatch(lowName,names,'exact');
  if length(k) == 1
    j = k;
  else
    allNames = ['(' Names{j(1),:}];
    for k = j(2:length(j))'
      allNames = [allNames ', ' Names{k,:}];
    end
    allNames = sprintf('%s).', allNames);
    error(message('globaloptim:gaoptimget2:ambiguousproperty',name,allNames));
  end
end

if any(strcmp(Names,Names{j,:}))
   o = options.(Names{j,:});
  if isempty(o)
    o = default;
  end
else
  o = default;
end

%------------------------------------------------------------------
function value = gaoptimgetfast(options,name,defaultopt)
%OPTIMGETFAST Get OPTIM OPTIONS parameter with no error checking so fast.
%   VAL = OPTIMGETFAST(OPTIONS,FIELDNAME,DEFAULTOPTIONS) will get the
%   value of the FIELDNAME from OPTIONS with no error checking or
%   fieldname completion. If the value is [], it gets the value of the
%   FIELDNAME from DEFAULTOPTIONS, another OPTIONS structure which is 
%   probably a subset of the options in OPTIONS.
%

if ~isempty(options)
        value = options.(name);
else
    value = [];
end

if isempty(value)
    value = defaultopt.(name);
end


