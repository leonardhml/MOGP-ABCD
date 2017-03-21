function nonNegScalar2(property,value)
%nonNegScalar2 any scalar >= 0    

%   Copyright 2007-2009 The MathWorks, Inc.

valid =  isreal(value) && isscalar(value) && (value >= 0);
if(~valid)
    error(message('globaloptim:nonNegScalar2:notNonNegScalar', property));
end
