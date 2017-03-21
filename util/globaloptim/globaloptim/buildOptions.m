function [options] = buildOptions(di, displayInterval, hybridFcn, outputFcns, plotFcns,  useParallel)
% di: 'off', 'final', 'iter'
% displayInterval: positive integer
% hybridFcn: [], @fmincon, @fminsearch, @patternsearch, @fminunc
% outputFcns: [], {@outfun1, @outfun2,...}
% plotFcns: [], {@plotfun1,...}
% useParallel: 1, 0

options.CreationFcn = @pswcreationuniform2;
options.Display = di;
options.DisplayInterval = displayInterval;
options.FunValCheck = 'off';
options.HybridFcn = hybridFcn;
options.InertiaRange = [0.1000 1.1000];
options.InitialSwarm = [];
options.InitialSwarmSpan = 2000;
options.MaxIter = '200*numberofvariables';
options.MaxTime = Inf;
options.MinFractionNeighbors = 0.2500;
options.ObjectiveLimit = -Inf;
options.OutputFcns = outputFcns;
options.PlotFcns = plotFcns;
options.SelfAdjustment = 1.4900;
options.SocialAdjustment = 1.4900;
options.StallIterLimit = 20;
options.StallTimeLimit = Inf;
options.SwarmSize = 'min(100,10*numberofvariables)';
options.TolFunValue = 1.0000e-06;
options.UseParallel = useParallel;
options.Vectorized = 'off';
options.TolFun = 1.0000e-06;
end
