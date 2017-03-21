function [handle,args] = functionHandleOrCell2(property,value)
%functionHandleOrCell2 A function Handle or a cell array starting with a function
%handle.

%   Copyright 2007-2011 The MathWorks, Inc.

[handle,args] = isFcn2(value);

if ~isempty(handle)
    return
elseif strcmp(property,'NonconFcn')
    error(message('globaloptim:functionHandleOrCell2:needFunctionHandleConstr'));
else
    error(message('globaloptim:functionHandleOrCell2:needFunctionHandle', property));
end

