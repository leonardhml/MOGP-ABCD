function positiveIntegerArray2(property,value)
%positiveIntegerArray2 positive integer array

%   Copyright 2007-2009 The MathWorks, Inc.

allValid = true;
for i = 1:numel(value)
    valid =  isreal(value(i)) && value(i) == floor(value(i)) && value(i) > 0;
    allValid = allValid && valid;
end
if(~valid)
    error(message('globaloptim:positiveIntegerArray2:notPosIntegerArray', property));
end
