function y = nonzerosign2(x)
%

%nonzerosign2 Signum function excluding zero
%   For each element of X, nonzerosign2(X) returns 1 if the element
%   is greater than or equal to zero and -1 if it is
%   less than zero. nonzerosign2 differs from SIGN in that nonzerosign2(0)
%   returns 1, while SIGN(0) returns 0.
%

%   Copyright 2008 The MathWorks, Inc.

y = ones(size(x));

y(x < 0) = -1;