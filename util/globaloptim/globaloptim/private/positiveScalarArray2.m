function positiveScalarArray2(property,value)
%positiveScalarArray2 positive scalar array

%   Copyright 2007-2009 The MathWorks, Inc.

allValid = true;
for i = 1:numel(value)
    valid =  isreal(value(i)) && value(i) > 0;
    allValid = allValid && valid;
end

if(~allValid)
    error(message('globaloptim:positiveScalarArray2:notPosScalarArray', property));
end
