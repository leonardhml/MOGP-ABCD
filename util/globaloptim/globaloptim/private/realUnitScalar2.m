function realUnitScalar2(property,value)
%realUnitScalar2 A scalar on the interval [0,1]

%   Copyright 2007-2009 The MathWorks, Inc.

valid = isreal(value) && isscalar(value) && (value >= 0) && (value <= 1);
if(~valid)
    error(message('globaloptim:realUnitScalar2:notScalarOnUnitInterval', property));
end
