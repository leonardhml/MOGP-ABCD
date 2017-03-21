function [feasible,XOUT] = getFeasiblePoints2(X,Aineq,bineq,Aeq,beq,lb,ub,tol)
%GETFEASIBLEPOINTS2 returns all feasible points in point matrix 'X' w.r.t. 
%   linear/bound constraints within the tolerance 'tol'.
%
%   Output arguments 'XOUT' is a matrix of all the feasible points and 'feasible'
%   is a logical array of length size(X,2) indicating the point is feasible (TRUE),
%   or infeasible (FALSE).
%
%   Example:
%     If there are 4 points in 2 dimension space, [2;-4],[1;5],[9;0]
%     and [-2;1] then
%
%     X  =   [2  1 9 -2
%            -4  5 0  1 ]
%
%     and if Aineq = diag([-2,2]), bineq = zeros(2,1);
%     tol = 1e-6; we obtain [Xout,feasible] = getFeasiblePoints2(X,Aineq,bineq,[],[],[],[],tol)
%
%     Xout = [-2;1]
%

%   Copyright 2003-2006 The MathWorks, Inc.


feasible = true(size(X,2),1);
XOUT=[];
for i = 1:size(X,2)
    feasible(i) = isTrialFeasible2(X(:,i),Aineq,bineq,Aeq,beq,lb,ub,tol);
    if feasible(i)
        XOUT(:,end+1) = X(:,i);
    end
end
