% Script to initialise a set of kernels and compute the respective model
% evidence(BIC) for each, in order to select the choice of kernel that best
% descrive the dataset.

% Create set of kernels using grammar (how many?)
% Create kernel descriptions for above set
% What we need as output: 
%   Model evidence of each model
%   Optimised hyp of each model
%   Kernel struct
disp('Building kernel set...');
kernel_set = buildKernelSet(4);

%%%
% Sample training and test inputs
% Data pre-processing goes here
%%%
% x1 = rand(100,1);                               % 100 training inputs
% y1 = sin(10*x1) + 0.1*gpml_randn(0.9, 100, 1);  % 100 noisy training targets
% x2 = rand(100,1);                               % 100 training inputs
% y2 = sin(9*x2) + 0.1*gpml_randn(0.9, 100, 1);   % 100 noisy training targets
% xs1 = linspace(0, 1, 20)';                      % 20 test inputs for output 1
% ys1 = sin(10*xs1) + 0.1*gpml_randn(0.9,20,1);
% X.x1 = x1;
% X.x2 = x2;
% Y.y1 = y1;
% Y.y2 = y2;

disp('Processing data...');
a = preprocessTrafficData('data/nus_traffic47/e103021891.csv', 1, 1, 12);
X.x1 = a(1:2:end, 1);
X.x2 = a(2:2:end, 1);
Y.y1 = a(1:2:end, 2);
Y.y2 = a(2:2:end, 3);

% Set up common (fixed) variables across all models
g1 = @(x, std1) gauss(x, std1);
g2 = @(x, std2, g2offset) gauss(x, std2, g2offset);
n = 6;
a = -2;
b = 2;

cov_options.g1 = g1;
cov_options.g2 = g2;
cov_options.n = n;
cov_options.a = a;
cov_options.b = b;

disp('Initiating search...');
for i = 1:length(kernel_set)
    fprintf('Kernel: %d\n', i);
    cov_options.k = kernel_set(i);
    model = MOGP(cov_options);

    disp('Optimising...');
    hyp_opt = model.optimise(X,Y);
    hyp.cov = hyp_opt(1:end-4);
    % hyp.cov = [0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1];
    hyp.smoothing = hyp_opt(end-3:end-2);
    % hyp.smoothing = [0.1;0.1;0.1];
    hyp.noise = hyp_opt(end-1:end);
    % hyp.noise = [0.1;0.1];
    
    model.fit(X,Y, hyp);
    bic = model.modelEvidence();
    results(i).bic = bic;
    results(i).k = k;
    results(i).hyp = hyp;
end