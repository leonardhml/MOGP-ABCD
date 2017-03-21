function LAMBDA = formBoxLambda2(grad, lb, ub)
%

%formBoxLambda2 Helper function that creates the lambda
% structure for the trust2-region-reflective algorithm.

%   Copyright 2015 The MathWorks, Inc.

% The trust2-region-reflective algorithm does not calculate the 
% Lagrange multipliers during the solution process. As such, they have to
% be computed separately.
%
% Compute the Lagrange Multipliers for the box-constrained nonlinear 
% problem:
%
%              min { F(x) :  lb <= x <= ub }
%
% where F:R^n -> R
%
% INPUTS:
%
%     grad - n-dimensional gradient for the function F() is provided in input of the function.
%       lb - n-dimensional lower-bound constraints
%       ub - n-dimensional upper-bound constraints


% The first order optimality conditions can be stated as 
%
% grad - lambda.lower + lambda.upper = 0
% lambda.lower(x - lb) = 0
% lambda.upper(ub - x) = 0
% lambda.lower, lambda.upper >= 0
%
% There are four cases to consider:
%
% Case 1: lb < x* < ub => lambda.lower = lambda.upper = 0
% Case 2: x* = lb < ub => lambda.lower = grad, lambda.upper = 0
% Case 3: x* = ub > lb => lambda.upper = -grad, lambda.lower = 0
% Case 4: x* = lb = ub: there are several options, and we chose
% the convention:
%   Case 4a: grad >= 0 => lambda.lower = grad, lambda.upper = 0
%   Case 4b: grad < 0  => lambda.upper = -grad, lambda.lower = 0
%
% If we were to implement the cases as written above, we would have to
% introduce an arbitrary tolerance to check whether a variable is at
% the bounds. To avoid such a tolerance, we note that we can rewrite
% the conditions as
%
% grad >= 0 => lambda.lower = grad, lambda.upper = 0
% grad < 0 => lambda.lower = 0, lambda.upper = -grad

% Initialize the Lagrange multipliers
LAMBDA.lower = zeros(length(lb),1);
LAMBDA.upper = zeros(length(ub),1);

% Convert gradient to full
grad = full(grad);

% Elements with non-negative gradient and finite lower bound. Note,
% that isfinite is used for the check because elements of lb could be
% -inf or NaN at this point.
idxNonNegGradAndFiniteLB = grad >= 0 & isfinite(lb);

% Set lower bound Lagrange multipliers for elements with positive
% gradient and finite lower bound.
LAMBDA.lower(idxNonNegGradAndFiniteLB) = grad(idxNonNegGradAndFiniteLB);

% Elements with negative gradient and finite upper bound. Note, that
% isfinite is used for the check because elements of u could be inf or
% NaN at this point.
idxNegGradAndFiniteUB = grad < 0 & isfinite(ub);

% Set upper bound Lagrange multipliers for elements with negative
% gradient and finite upper bound
LAMBDA.upper(idxNegGradAndFiniteUB) = -grad(idxNegGradAndFiniteUB);