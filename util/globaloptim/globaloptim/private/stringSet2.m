function stringSet2(property,value,validStrs)
%stringSet2 Check that a given string is one of a valid set of strings.

%   Copyright 2007-2015 The MathWorks, Inc.

if ~ischar(value) || ~any(strcmpi(value,validStrs))
    error(message('globaloptim:stringSet2:notCorrectChoice', ...
        'OPTIONS',property,strjoin(validStrs,', ') ));
end