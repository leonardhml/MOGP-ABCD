function optionDefault = getOptionDefaultValue2(OS, OptionStoredName)
%

%getOptionDefaultValue2 Get a default value for an option
%
%   OPTIONDEFAULT = getOptionDefaultValue2(OS, OPTIONSTOREDNAME) gets a
%   default value for an option. This function returns a default value even
%   if the option is not used by the current algorithm.

% Copyright 2015 The MathWorks, Inc.

if isfield(OS.AlgorithmDefaults{OS.AlgorithmIndex}, OptionStoredName)
    optionDefault = OS.AlgorithmDefaults{OS.AlgorithmIndex}.(OptionStoredName);
else
    RemAlgs = find(~OS.AlgorithmIndex);
    if OS.IsConstantDefault.(OptionStoredName)
        for j = RemAlgs
            if isfield(OS.AlgorithmDefaults{j}, OptionStoredName)
                optionDefault = OS.AlgorithmDefaults{j}.(OptionStoredName);
                break
            end
        end
    else
        optionDefault = OS.NonConstantDefaults.(OptionStoredName){OS.AlgorithmIndex};
    end
end
