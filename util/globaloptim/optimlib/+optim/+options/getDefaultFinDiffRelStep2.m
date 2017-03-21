function FinDiffRelStep = getDefaultFinDiffRelStep2(FinDiffType)
%getDefaultFinDiffRelStep2 Get the default for FinDiffRelStep
%
%   FINDIFFRELSTEP = getDefaultFinDiffRelStep2(FINDIFFTYPE) returns the
%   FinDiffRelStep default value given the finite difference type.

%   Copyright 2013 The MathWorks, Inc.

switch FinDiffType
    case 'forward'
        FinDiffRelStep = 'sqrt(eps)';
    case 'central'
        FinDiffRelStep = 'eps^(1/3)';
end
