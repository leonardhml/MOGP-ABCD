function realScalar2(property,value)
%realScalar2 Test for real scalar

%   Copyright 2007-2009 The MathWorks, Inc.

valid = isreal(value) && isscalar(value);
if(~valid)
    error(message('globaloptim:realScalar2:notScalar', property));
end
