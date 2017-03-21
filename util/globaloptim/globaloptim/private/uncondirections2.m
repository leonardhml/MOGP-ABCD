function Basis = uncondirections2(AdaptiveMesh,MeshSize,x)
%UNCONDIRECTIONS2 Finds search2 vectors when no constraints are present.
%
%   AdaptiveMesh is a boolean value which is true if MADS is used to poll2
%   and X is the current point.

%   Copyright 2003-2007 The MathWorks, Inc.


vars = length(x);
% N linearly independent vectors
if ~AdaptiveMesh
    Basis  = eye(vars);
else
    pollParam = 1/sqrt(MeshSize);
    lowerT = tril((round((pollParam+1)*rand(vars)-0.5)),-1);
    diagtemp = pollParam*sign(rand(vars,1)-0.5);
    diagtemp(~diagtemp) = pollParam*sign(0.5-rand);
    diagT  = diag(diagtemp);
    Basis = lowerT + diagT;
    order = randperm(vars);
    Basis = Basis(order,order);
end
