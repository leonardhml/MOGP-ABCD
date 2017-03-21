function [x,lb,ub,msg] = checkbounds2(xin,lbin,ubin,nvars)
%

%checkbounds2 Verify that the bounds are valid with respect to initial point.
%
% This is a helper function.

%   [X,LB,UB,X,FLAG] = checkbounds2(X0,LB,UB,nvars)
%   checks that the upper and lower
%   bounds are valid (LB <= UB) and the same length as X (pad with -inf/inf
%   if necessary); warn if too short or too long.  Also make LB and UB vectors
%   if not already. Finally, inf in LB or -inf in UB throws an error.

%   Copyright 1990-2012 The MathWorks, Inc.

msg = [];
% Turn into column vectors
lb = lbin(:);
ub = ubin(:);
xin = xin(:);

lenlb = length(lb);
lenub = length(ub);

% Check lb length
if lenlb > nvars
    warning(message('optimlib:checkbounds2:IgnoringExtraLbs'));
    lb = lb(1:nvars);
    lenlb = nvars;
elseif lenlb < nvars % includes empty lb case
    if lenlb > 0
        % lb is non-empty and length(lb) < nvars.
        warning(message('optimlib:checkbounds2:PadLbWithMinusInf'));
    end
    lb = [lb; -inf*ones(nvars-lenlb,1)];
    lenlb = nvars;
end

% Check ub length
if lenub > nvars
    warning(message('optimlib:checkbounds2:IgnoringExtraUbs'));
    ub = ub(1:nvars);
    lenub = nvars;
elseif lenub < nvars % includes empty ub case
    if lenub > 0
        % ub is non-empty and length(ub) < nvars.
        warning(message('optimlib:checkbounds2:PadUbWithInf'));
    end
    ub = [ub; inf*ones(nvars-lenub,1)];
    lenub = nvars;
end

% Check feasibility of bounds
len = min(lenlb,lenub);
if any( lb( (1:len)' ) > ub( (1:len)' ) )
    count = full(sum(lb>ub));
    if count == 1
        msg=sprintf(['Exiting due to infeasibility:  %i lower bound exceeds the' ...
            ' corresponding upper bound.'],count);
    else
        msg=sprintf(['Exiting due to infeasibility:  %i lower bounds exceed the' ...
            ' corresponding upper bounds.'],count);
    end
end

% Check if -inf in ub or inf in lb
if any(eq(ub, -inf))
    error(message('optimlib:checkbounds2:MinusInfUb'));
elseif any(eq(lb,inf))
    error(message('optimlib:checkbounds2:PlusInfLb'));
end

x = xin;
