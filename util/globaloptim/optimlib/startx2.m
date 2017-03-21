function xstart = startx2(ub,lb,xstart,xstartOutOfBounds_idx)
%

%startx2	Box-centered start point.
%
% This is a helper function.

% xstart = startx2(ub,lb,xstart,xstartOutOfBounds_idx) sets the components
% that violate the bounds to a centered value.

%   Copyright 1990-2011 The MathWorks, Inc.

if nargin == 2
    xstartOutOfBounds_idx = [];
    xstart = [];
elseif nargin ~= 4 % neither 2 nor 4 inputs
    error(message('optimlib:startx2:WrongInputsXstart'))
end
    
arg1 = (ub < inf)  & (lb == -inf);
arg2 = (ub == inf) & (lb > -inf);
arg3 = (ub < inf)  & (lb > -inf);
arg4 = (ub == inf) & (lb == -inf);

if ~isempty(xstart)
    arg1 = arg1 & xstartOutOfBounds_idx;
    arg2 = arg2 & xstartOutOfBounds_idx;
    arg3 = arg3 & xstartOutOfBounds_idx;
    arg4 = arg4 & xstartOutOfBounds_idx;    
else
    n = length(ub);
    xstart = zeros(n,1);
end

%
w = max(abs(ub),1);
xstart(arg1) = ub(arg1) - .5*w(arg1);
%
ww = max(abs(lb),1);
xstart(arg2) = lb(arg2) + .5*ww(arg2);
%
xstart(arg3) = (ub(arg3) + lb(arg3))/2;
%
xstart(arg4) = 1;