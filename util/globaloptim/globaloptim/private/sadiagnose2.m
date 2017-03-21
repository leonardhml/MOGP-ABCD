function sadiagnose2(options,problem)
% SADIAGNOSE2 prints some diagnostic information about the problem

% Copyright 2006-2010 The MathWorks, Inc.

%   This function is private to SIMULANNEAL2.

properties =  optionsList2('sa');
Output_String = sprintf('\nDiagnostic information.');

Output_String = [Output_String sprintf('\n\tobjective function = %s',value2RHS2(problem.objective))];

Output_String = [Output_String sprintf('\n\tX0 = %s',value2RHS2(problem.x0))];

Output_String = [Output_String sprintf('\n%s','Modified options:')];
for i = 1:length(properties)
    prop = properties{i};
    if(~isempty(prop)) % the property list has blank lines, ignore them
        value = options.(prop);
        if ~isempty(value)  % don't generate Output_String for defaults.
            Output_String = [Output_String sprintf('\n\toptions.%s = %s',prop,value2RHS2(value))];
        end
    end
end
Output_String = [Output_String sprintf('\nEnd of diagnostic information.\n\n')];
fprintf('%s',Output_String)
