function isFeas = isNonlinearFeasible2(C,Ceq,tol)
%isNonlinearFeasible2 Test for nonlinear feasibility with respect to a
%tolerance.

%   Copyright 2014 The MathWorks, Inc.

isFeas = all([C abs(Ceq)] <= tol, 2);
