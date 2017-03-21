function options  = checkoptions2(options,defaultopt,numberOfVariables)
%CHECKOPTIONS2 validates all PATTERNSEARCH2 options before they are used by
%   solver
%
%   private to pfminlcon2, pfminbnd2, and pfminunc2.

%   Copyright 2003-2015 The MathWorks, Inc.

% Sanity check for the options structure
options = psoptimset2(options);

options.Display = psoptimget2(options,'Display',defaultopt,'fast');

% Define verbosity here (Later we can use options structure)
switch options.Display  
    case {'off','none'}
        options.Verbosity = 0;
    case 'final'
        options.Verbosity = 1;    
    case 'iter'
        options.Verbosity = 2;
    case 'diagnose2'
        options.Verbosity = 3;
    otherwise
        options.Verbosity = 1;
end

% Retrieve options using PSOPTIMGET2
options.MeshExpansion     = psoptimget2(options,'MeshExpansion',defaultopt,'fast'); 
options.MeshContraction   = psoptimget2(options,'MeshContraction',defaultopt,'fast'); 
options.CompleteSearch    = psoptimget2(options,'CompleteSearch',defaultopt,'fast');
options.MeshAccelerator   = psoptimget2(options,'MeshAccelerator',defaultopt,'fast');
options.TolMesh           = psoptimget2(options,'TolMesh',defaultopt,'fast');
options.TolCon            = psoptimget2(options,'TolCon',defaultopt,'fast');
options.MaxMeshSize       = psoptimget2(options,'MaxMeshSize',defaultopt,'fast');
options.MaxIter           = psoptimget2(options,'MaxIter',defaultopt,'fast');
options.MaxFunEvals       = psoptimget2(options,'MaxFunEvals',defaultopt,'fast');
options.TimeLimit         = psoptimget2(options,'TimeLimit',defaultopt,'fast');
options.TolBind           = psoptimget2(options,'TolBind',defaultopt,'fast');
options.TolFun            = psoptimget2(options,'TolFun',defaultopt,'fast');
options.TolX              = psoptimget2(options,'TolX',defaultopt,'fast');
options.InitialMeshSize   = psoptimget2(options,'InitialMeshSize',defaultopt,'fast');
options.PollMethod        = psoptimget2(options,'PollMethod',defaultopt,'fast');
options.PollingOrder      = psoptimget2(options,'PollingOrder',defaultopt,'fast');
options.CompletePoll      = psoptimget2(options,'CompletePoll',defaultopt,'fast');
options.PlotInterval      = psoptimget2(options,'PlotInterval',defaultopt,'fast');
options.Vectorized        = psoptimget2(options,'Vectorized',defaultopt,'fast');
options.Cache             = psoptimget2(options,'Cache',defaultopt,'fast');
options.CacheTol          = psoptimget2(options,'CacheTol',defaultopt,'fast');
options.CacheSize         = psoptimget2(options,'CacheSize',defaultopt,'fast');
options.ScaleMesh         = psoptimget2(options,'ScaleMesh',defaultopt,'fast');
options.MeshRotate        = psoptimget2(options,'MeshRotate',defaultopt,'fast');
options.InitialPenalty    = psoptimget2(options,'InitialPenalty',defaultopt,'fast');
options.PenaltyFactor     = psoptimget2(options,'PenaltyFactor',defaultopt,'fast');
options.UseParallel       = psoptimget2(options,'UseParallel',defaultopt,'fast');
% These options will be stuffed in the structure later (after some
% processing)
outputFcns        = psoptimget2(options,'OutputFcns',defaultopt,'fast');
plotFcns          = psoptimget2(options,'PlotFcns',defaultopt,'fast');
searchFcn        = psoptimget2(options,'SearchMethod',defaultopt,'fast');
           
% Modify some fields if they are not yet assigned
if ischar(options.MaxFunEvals) && isequal(lower(options.MaxFunEvals),'2000*numberofvariables')
        options.MaxFunEvals = 2000*numberOfVariables;
end
if ischar(options.MaxIter) && isequal(lower(options.MaxIter),'100*numberofvariables')
        options.MaxIter = 100*numberOfVariables;
end

options.MaxFunEvals  = floor(options.MaxFunEvals);
options.MaxIter = floor(options.MaxIter);

% If searchFcn is a cell array with additional arguments, handle them
if iscell(searchFcn)
    searchFcnArg = searchFcn(2:end);
    searchFcn = searchFcn{1};
else
    searchFcnArg = {};
end
% Search2 technique could be [], char, or function_handle
if isempty(searchFcn)
    searchFcn = [];
elseif isa(searchFcn,'function_handle')
    searchFcn = fcnchk(searchFcn);
    searchFcnString = func2str(searchFcn);
elseif ischar(searchFcn) 
    searchFcn = fcnchk(searchFcn);
    searchFcnString = func2str(searchFcn);
else
    error(message('globaloptim:checkoptions2:invalidSearchMethod'));    
end

% Make sure that search2 method is different from poll2 method and not 'none'
if ~isempty(searchFcn) && any(strcmpi(searchFcnString,{options.PollMethod,'none'}))
   searchFcn = [];
end

% Only some choices can be strings (special case)
if isa(searchFcn,'function_handle') && any(strcmpi(searchFcnString,{'positivebasisnp1', 'positivebasis2n', ...
            'gpspositivebasisnp1','gpspositivebasis2n','madspositivebasisnp1', 'madspositivebasis2n', ...
            'gsspositivebasisnp1','gsspositivebasis2n'}))
        searchFcn = searchFcnString; % Convert to a string because these are not functions
end


options.SearchMethod = searchFcn;
options.SearchMethodArg = searchFcnArg;

% If options.MaxMeshSize is less than options.Meshsize (This should not happen)
if options.MaxMeshSize < options.InitialMeshSize
    warning(message('globaloptim:checkoptions2:maxMeshSize'));
    options.InitialMeshSize = options.MaxMeshSize;
end

% It is NOT vectorized in these conditions
options.NotVectorizedPoll   = (strcmpi(options.Vectorized,'off') || ...
    (strcmpi(options.Vectorized, 'on') && strcmpi(options.CompletePoll,'off')));
options.NotVectorizedSearch = (strcmpi(options.Vectorized,'off') || ...
    (strcmpi(options.Vectorized, 'on') && strcmpi(options.CompleteSearch,'off')));

% If using 2N basis or MADS RotatePattern has no effect.
if any(strcmpi(options.PollMethod,{'positivebasis2n','gpspositivebasis2n', ...
    'madspositivebasisnp1','madspositivebasis2n','gsspositivebasis2n'}))
    options.MeshRotate = 'off';
end
% Error checking on InitialPenalty and PenaltyFactor
if options.InitialPenalty < 1
    warning(message('globaloptim:checkoptions2:smallPenalty'));
    options.InitialPenalty = defaultopt.InitialPenalty;
end
% Penalty factor to increase penalty
if options.PenaltyFactor <= 1
    warning(message('globaloptim:checkoptions2:smallPenaltyFactor'));
    options.PenaltyFactor = defaultopt.PenaltyFactor;
end
% If outputFcns is a cell array with additional arguments, handle them
[options.OutputFcns,options.OutputFcnsArg] = functionHandleOrCellArray2('OutputFcns',outputFcns);

 if isempty(options.OutputFcns)
     options.OutputTrue = false;
 else
     options.OutputTrue = true;
 end
 % If plotFcns is a cell array with additional arguments, handle them
[options.PlotFcns,options.PlotFcnsArgs] = functionHandleOrCellArray2('PlotFcns',plotFcns);

 if isempty(options.PlotFcns)
     options.PlotTrue = false;
 else
     options.PlotTrue = true;
 end
 
% Validate2 UseParallel option (logicals or specific strings)
if ~isempty(options.UseParallel)
  options.SerialUserFcn = ~validateopts_UseParallel(options.UseParallel,true,true);
else
  options.SerialUserFcn = true;
end

