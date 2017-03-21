function o = saoptimget2(options,name,default,flag)
%SAOPTIMGET2 Get SIMULANNEALBND2 OPTIONS parameter value.
%   VAL = SAOPTIMGET2(OPTIONS,'NAME') extracts the value of the named
%   parameter from optimization options structure OPTIONS, returning an
%   empty matrix if the parameter value is not specified in OPTIONS.  It is
%   sufficient to type only the leading characters that uniquely identify
%   the parameter.  Case is ignored for parameter names.  [] is a valid
%   OPTIONS argument.
%   
%   VAL = SAOPTIMGET2(OPTIONS,'NAME',DEFAULT) extracts the named parameter
%   as above, but returns DEFAULT if the named parameter is not specified
%   (is []) in OPTIONS.  For example
%   
%     opt = saoptimset2('TolFun',1e-4);
%     val = saoptimget2(opt,'TolFun');
%   
%   returns val = 1e-4.
%   
%   See also SAOPTIMSET2.

%   Copyright 2006-2015 The MathWorks, Inc.

if nargin < 2
  error(message('globaloptim:saoptimget2:inputarg'));
end
if nargin < 3
  default = [];
end
if nargin < 4
   flag = [];
end

% undocumented usage for fast access with no error checking
if isequal('fast',flag)
   o = saoptimgetfast(options,name,default);
   return
end

if ~isempty(options) && ~isa(options,'struct')
  error(message('globaloptim:saoptimget2:firstargerror'));
end

if isempty(options)
  o = default;
  return;
end

optionsstruct = struct('AnnealingFcn', [], ...
'TemperatureFcn', [], ...
'AcceptanceFcn', [], ...
'TolFun', [], ...
'StallIterLimit', [], ...
'MaxFunEvals', [], ...
'TimeLimit', [], ...
'MaxIter', [], ...
'ObjectiveLimit', [], ...
'Display', [], ...
'DisplayInterval', [], ...
'HybridFcn', [], ...
'HybridInterval', [], ...
'PlotFcns', [], ...
'PlotInterval', [], ...
'OutputFcns', [], ...
'InitialTemperature', [], ...
'ReannealInterval', [], ...
'DataType', []);
 
Names = fieldnames(optionsstruct);
%[m,n] = size(Names);
names = lower(Names);

lowName = lower(name);
j = strmatch(lowName,names);
if isempty(j) % if no matches
  error(message('globaloptim:saoptimget2:invalidProperty','name'));
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
    error(message('globaloptim:saoptimget2:ambiguousProperty',name,allNames));
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
function value = saoptimgetfast(options,name,defaultopt)
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


