function psdiagnose2(FUN,nonlcon,Xin,nineqcstr,neqcstr,ncstr,options)
%PSDIAGNOSE2 prints some diagnostic information about the problem
%   private to PFMINCON2, PFMINLCON2, PFMINBND2, PFMINUNC2.

%   Copyright 2003-2006 The MathWorks, Inc.

properties =  optionsList2('ps');
defaultOpt = psoptimset2;
Output_String = sprintf('\nDiagnostic information.');

Output_String = [Output_String sprintf('\n\tobjective function = %s',value2RHS2(FUN))];
%print some information about constraints
if ~isempty(nonlcon)
    Output_String = [Output_String sprintf('\n\tnonlinear constraint function = %s',value2RHS2(nonlcon))];
end
if ~isempty(nineqcstr)
    Output_String = [Output_String sprintf('\n\t%d Inequality constraints',nineqcstr)];
end
if ~isempty(neqcstr)
    Output_String = [Output_String sprintf('\n\t%d Equality constraints',neqcstr)];
end
if ~isempty(ncstr)
    Output_String = [Output_String sprintf('\n\t%d Total number of linear constraints\n',ncstr)];
end

Output_String = [Output_String sprintf('\n\tX0 = %s',value2RHS2(Xin))];
%Output_String = [Output_String sprintf('\n\tf(X0) = %s',value2RHS2(Iterate.f))];

Output_String = [Output_String sprintf('\n%s','Modified options:')];
for i = 1:length(properties)
    prop = properties{i};
    if(~isempty(prop)) % the property list has blank lines, ignore them
        value = options.(prop);
        if ~isempty(value) && ~isequal(value,defaultOpt.(prop)) % don't generate Output_String for defaults.
            Output_String = [Output_String sprintf('\n\toptions.%s = %s',prop,value2RHS2(value))];
        end
    end
end
Output_String = [Output_String sprintf('\nEnd of diagnostic information.')];
fprintf('%s',Output_String)
