function positiveInteger2(property,value)
%positiveInteger2 any positive integer

%   Copyright 2007-2009 The MathWorks, Inc.

valid =  isreal(value) && isscalar(value) && (value > 0) && (value == floor(value));
if(~valid)
   error(message('globaloptim:positiveInteger2:notPosInteger', property));
end

