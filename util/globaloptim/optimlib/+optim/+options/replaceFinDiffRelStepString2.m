function options = replaceFinDiffRelStepString2(options)
%replaceFinDiffRelStepString2 Replace string value in FinDiffRelStep
%
%   OPTIONS = replaceFinDiffRelStepString2(OPTIONS) replaces a string value
%   in FinDiffRelStep with the equivalent numerical value.

%   Copyright 2013 The MathWorks, Inc.

if ischar(options.FinDiffRelStep)
    % At this point we know that FinDiffRelStep will either be
    % 'sqrt(eps)' or 'eps^(1/3)'
    idx = strcmpi(options.FinDiffRelStep, {'sqrt(eps)', 'eps^(1/3)'});
    possValues = [sqrt(eps), eps^(1/3)];
    options.FinDiffRelStep = possValues(idx);
end
