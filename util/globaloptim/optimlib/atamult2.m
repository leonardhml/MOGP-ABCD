function[V] = atamult2(A,Y,flag,varargin)
%

%atamult2 Jacobian-matrix multiply
%
%	V = atamult2(A,Y) computes V = (A'*(A*Y)).
%
%	V = atamult2(A,Y,flag) computes V = (A'*(A*Y)) if flag = 0,
%                                  V = A*Y        if flag > 0,
%                                  V = A'*Y       if flag < 0.
%
% Note: varargin is not used but must be provided in case 
% the objective function has additional problem dependent
% parameters (which will be passed to this routine as well).

%   Copyright 1990-2006 The MathWorks, Inc.

if nargin < 3 || flag == 0
   V = (A'*(A*Y));
elseif flag > 0
   V = A*Y;
else
   V = A'*Y;
end


