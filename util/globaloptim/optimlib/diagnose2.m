function msg = diagnose2(caller,OUTPUT,gradflag,hessflag,constflag,gradconstflag, ...
    XOUT,non_eq,non_ineq,lin_eq,lin_ineq,LB,UB,funfcn,confcn)
%

%diagnose2 prints diagnostic information about the function to be minimized
%    or solved.
%
% This is a helper function.

%   Copyright 1990-2014 The MathWorks, Inc. 

msg = [];

beginStr = getString(message('optimlib:diagnose2:DiagnosticHeader'));
separatorTxt = repmat('_',60,1);
fprintf('\n%s\n   %s\n\n',separatorTxt,beginStr);

if ~isempty(funfcn{1})
    funformula =  getformula(funfcn{3});
    gradformula = getformula(funfcn{4});
    hessformula = getformula(funfcn{5});
else
    funformula =  '';
    gradformula = '';
    hessformula = '';
end

if ~isempty(confcn{1})
    conformula = getformula(confcn{3});
    gradcformula = getformula(confcn{4});
else
    conformula = '';
    gradcformula = '';
end    

fprintf(getString(message('optimlib:diagnose2:NumVars',length(XOUT))))
if ~isempty(funfcn{1})
    fprintf(getString(message('optimlib:diagnose2:FuncHeader')))
    switch funfcn{1}
    case 'fun'
        % display 
        fprintf(getString(message('optimlib:diagnose2:ObjFun',funformula)));
        
    case 'fungrad'
        if gradflag
            fprintf(getString(message('optimlib:diagnose2:ObjGradFun',funformula)));
        else
            fprintf(getString(message('optimlib:diagnose2:ObjFun',funformula)));
            fprintf(getString(message('optimlib:diagnose2:SetGradObj'))) 
        end
        
    case 'fungradhess'
        if gradflag && hessflag
            fprintf(getString(message('optimlib:diagnose2:ObjGradHessFun',funformula)));
        elseif gradflag
            fprintf(getString(message('optimlib:diagnose2:ObjGradFun',funformula)));
            fprintf(getString(message('optimlib:diagnose2:SetHess'))) 
        else
            fprintf(getString(message('optimlib:diagnose2:ObjFun',funformula)));
            fprintf(getString(message('optimlib:diagnose2:SetGradObj')))
            fprintf(getString(message('optimlib:diagnose2:SetHess'))) 
        end
        
        
    case 'fun_then_grad'
        fprintf(getString(message('optimlib:diagnose2:ObjFun',funformula)));
        if gradflag
            fprintf(getString(message('optimlib:diagnose2:GradFun',gradformula)));
        end
        if hessflag
            fprintf(getString(message('optimlib:diagnose2:IgnoreHessOpt')))
        end
        
    case 'fun_then_grad_then_hess'
        fprintf(getString(message('optimlib:diagnose2:ObjFun',funformula)));
        if gradflag && hessflag
            fprintf(getString(message('optimlib:diagnose2:GradFun',gradformula)));
            fprintf(getString(message('optimlib:diagnose2:HessFun',hessformula,'')));
        elseif gradflag
            fprintf(getString(message('optimlib:diagnose2:GradFun',gradformula)));
        end   
    otherwise
        
    end
    
    if ~gradflag
        fprintf(getString(message('optimlib:diagnose2:GradFun',getString(message('optimlib:diagnose2:FinDiff')))))
    end
    % shape of grad
    
    if ~hessflag && (isequal('fmincon2',caller) || isequal('constrsh',caller) || isequal('fminunc',caller))
        fprintf(getString(message('optimlib:diagnose2:HessFun', ...
            getString(message('optimlib:diagnose2:FinDiff')), ...
            getString(message('optimlib:diagnose2:QuasiNewton')) )));
    end
    % shape of hess
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(confcn{1})
    switch confcn{1}
        
    case 'fun'
        fprintf(getString(message('optimlib:diagnose2:NonlinConstr',conformula)));
    case 'fungrad'
        if gradconstflag
            fprintf(getString(message('optimlib:diagnose2:NonlinConstrAndGrad',conformula)));
        else
            fprintf(getString(message('optimlib:diagnose2:NonlinConstr',conformula)));
            fprintf(getString(message('optimlib:diagnose2:SetGradConstr'))) 
        end
        
    case 'fun_then_grad'
        fprintf(getString(message('optimlib:diagnose2:NonlinConstr',conformula)));
        if gradconstflag
            fprintf(getString(message('optimlib:diagnose2:NonlinConstGrad',gradcformula)));
        end
        
    otherwise
        
    end
    
    if ~constflag
        fprintf(getString(message('optimlib:diagnose2:NonlinConstr', ...
            getString(message('optimlib:diagnose2:DoNotExist')) )));
    end
    if ~gradconstflag
        fprintf(getString(message('optimlib:diagnose2:NonlinConstGrad',...
            getString(message('optimlib:diagnose2:FinDiff')) )));
    end
    fprintf(getString(message('optimlib:diagnose2:ConstrHeader')))  
    fprintf(getString(message('optimlib:diagnose2:NumNonlinIneq',non_ineq)))
    fprintf(getString(message('optimlib:diagnose2:NumNonlinEq',non_eq)))
    
elseif isequal(caller,'fmincon2') || isequal(caller,'constrsh') || isequal(caller,'fminimax') || ...
        isequal(caller,'fgoalattain') || isequal(caller,'fseminf')
    fprintf(getString(message('optimlib:diagnose2:ConstrHeader')));
    fprintf(getString(message('optimlib:diagnose2:NonlinConstr', ...
        getString(message('optimlib:diagnose2:DoNotExist')))));
    
end

fprintf(' \n')

switch caller
case {'fmincon2','constrsh','linprog','quadprog','lsqlin2','fminimax','fseminf','fgoalattain'}
    fprintf(getString(message('optimlib:diagnose2:NumLinIneq',lin_ineq)))
    fprintf(getString(message('optimlib:diagnose2:NumLinEq',lin_eq)))
    fprintf(getString(message('optimlib:diagnose2:NumLB',nnz(~isinf(LB)))))
    fprintf(getString(message('optimlib:diagnose2:NumUB',nnz(~isinf(UB)))))
case {'lsqcurvefit2','lsqnonlin'}
    fprintf(getString(message('optimlib:diagnose2:NumLB',nnz(~isinf(LB)))))
    fprintf(getString(message('optimlib:diagnose2:NumUB',nnz(~isinf(UB)))))
case {'fsolve','fminunc','fsolves'}
otherwise
end

if ~isempty(OUTPUT)
    fprintf(getString(message('optimlib:diagnose2:AlgorithmHeader',OUTPUT.algorithm)));
end

endStr = getString(message('optimlib:diagnose2:EndDiagnose'));
fprintf('\n%s\n   %s\n',separatorTxt,endStr);


%--------------------------------------------------------------------------------
function funformula = getformula(fun)
% GETFORMULA Convert FUN to a string.

if isempty(fun)
    funformula = '';
    return;
end

if ischar(fun) % already a string
    funformula = fun;
elseif isa(fun,'function_handle')  % function handle
    funformula = func2str(fun);
elseif isa(fun,'inline')   % inline object
    funformula = formula(fun);
else % something else with a char method
    try
        funformula = char(fun);
    catch
        funformula = '';
    end
end
