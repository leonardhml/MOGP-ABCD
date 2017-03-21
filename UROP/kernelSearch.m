% Script to initialise a set of kernels and compute the respective model
% evidence(BIC) for each, in order to select the choice of kernel that best
% descrive the dataset.

%%%
% Sample training and test inputs
%%%
[X,Y,xpred,ypred,outputpred] = generateData();

% Set up common (fixed) variables across all models
g1 = @(x, std1) gauss(x, std1);
g2 = @(x, std2) gauss(x, std2);
n = 6;

% a and b represent the range of input space
a = 0;
b = 1;

cov_options.g1 = g1;
cov_options.g2 = g2;
cov_options.n = n;
cov_options.a = a;
cov_options.b = b;

% Create set of kernels using grammar (how many?)
% Create kernel descriptions for above set
% What we need as output: 
%   Model evidence of each model
%   Optimised hyp of each model
%   Kernel struct
% disp('Building kernel set...');
% kernel_set = buildKernelSet(4);
% 
% disp('Initiating search...');
% results = {};
% for i = 1:length(kernel_set)
%     fprintf('Kernel: %d\n', i);
%     cov_options.k = kernel_set(i);
%     model = MOGP(cov_options);
% 
%     disp('Optimising...');
%     hyp_opt = model.optimise(X,Y);
%     hyp.cov = hyp_opt(1:end-4);
%     hyp.smoothing = hyp_opt(end-3:end-2);
%     hyp.noise = hyp_opt(end-1:end);
%     model.fit(X,Y, hyp);
%     bic = model.modelEvidence();
%     results(i).bic = bic;
%     results(i).k = kernel_set(i);
%     results(i).hyp = hyp_opt;
% end

% Greedy kernel search procedure
depth = 1;
maxDepth = 2;
bestForDepth = struct('bic', inf, 'k', []);
for i = 1:maxDepth
    bestForDepth(i).bic = inf;
    bestForDepth(i).k = [];
end
results = {};

kernel_set = expandKernel();
from = 1;
to = from + length(kernel_set) - 1;
while (depth <= maxDepth)
    fprintf('Exploring depth: %d\n', depth);
    parfor i = from:to
        pause(2);
        fprintf('Kernel: %d\n', i);
        
        cov_options = [];
        cov_options.g1 = g1;
        cov_options.g2 = g2;
        cov_options.n = n;
        cov_options.a = a;
        cov_options.b = b;
        cov_options.k = kernel_set(i - from + 1);
        
        model = MOGP(cov_options);

        hyp = [];
        [hyp_opt, fval] = model.optimise(X,Y);
        hyp.cov = hyp_opt(1:end-4);
        hyp.smoothing = hyp_opt(end-3:end-2);
        hyp.noise = hyp_opt(end-1:end);

        model.fit(X,Y, hyp);
        bic = model.modelEvidence();
        results(i).bic = bic;
        results(i).fval = fval;
        results(i).k = kernel_set(i - from + 1);
        results(i).hyp = hyp_opt;
        results(i).depth = depth;
        
        pause(2);
        fprintf('Optimisation done for Kernel: %d\n', i);
    end
    
    % Find best kernel for a particular depth
    for i = from:to
        if results(i).bic < bestForDepth(depth).bic
            bestForDepth(depth).k = results(i).k;
            bestForDepth(depth).bic = results(i).bic;
        end
    end
    
    fprintf('Done exploring depth: %d\n', depth);
    
    kernel_set = expandKernel(bestForDepth(depth).k);
    from = to + 1;
    to = from + length(kernel_set) - 1;
    depth = depth + 1;
    
   
    fprintf('Saving results\n'); 
    save('results.mat', 'results');
end
fprintf('Search complete\n');