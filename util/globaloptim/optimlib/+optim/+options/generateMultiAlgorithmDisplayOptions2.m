function OS = generateMultiAlgorithmDisplayOptions2(OS, solverName)

% Set the DisplayOptions
mc = meta.class.fromName(solverName);
numAlgs = length(OS.AlgorithmDefaults);
OS.NumDisplayOptions = zeros(1, numAlgs);
for i = 1:numAlgs
    OS.DisplayOptions{i} = cell(1, 0);
    thisAlgOptions = fieldnames(OS.AlgorithmDefaults{i});
    
    for j = 1:length(thisAlgOptions)
        thisName = optim.options.OptionAliasStore2.getNameFromAlias(thisAlgOptions{j});
        
        for k = 1:length(thisName)
            displayThisOption = true;
            foundOption = false;
            
            for l = 1:length(mc.PropertyList)
                if strcmp(mc.PropertyList(l).Name, thisName{k})
                    foundOption = true;
                    if mc.PropertyList(l).Hidden
                        displayThisOption = false;
                    end
                    break
                end
            end
            
            if foundOption && displayThisOption
                OS.DisplayOptions{i}(end+1) = {thisName{k}};
            end
        end
    end

    % Add the Algorithm option
    OS.DisplayOptions{i}{end+1} = 'Algorithm';
    
    OS.NumDisplayOptions(i) = length(OS.DisplayOptions{i});
end
