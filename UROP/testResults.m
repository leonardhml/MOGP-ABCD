function [ output_args ] = testResults( result )

%%%
% Sample training and test inputs
%%%
[X,Y,xpred,ypred,outputpred] = generateData();

% A sample k
k.kernel = result.k.kernel;
k.components = result.k.components;

% sample g1, g2
g1 = @(x, std1) gauss(x, std1);
g2 = @(x, std2) gauss(x, std2);

n = 6;
a = 0;
b = 1;

% cov_options is a struct that define the options for the kernel,
% such as the latent kernel k, hyperparameters, smoothing kernels, and
% sampling parameters
cov_options.k = k;
cov_options.g1 = g1;
cov_options.g2 = g2;
cov_options.n = n;
cov_options.a = a;
cov_options.b = b;
model = MOGP(cov_options);

% Find MLE estimate for hyp
hyp_opt = result.hyp;
hyp.cov = hyp_opt(1:end-4);
hyp.smoothing = hyp_opt(end-3:end-2);
hyp.noise = hyp_opt(end-1:end);

% Fit and predict
model.fit(X,Y,hyp);
[mu, s2] = model.predict(xpred, 1);
fprintf('Model Evidence is %d.\n', model.modelEvidence);

% Visualise
f = [mu+2*sqrt(diag(s2)); flipdim(mu-2*sqrt(diag(s2)),1)];
fill([xpred; flipdim(xpred,1)], f, [7 7 7]/8)
hold on; plot(xpred, mu, 'b'); plot(X.x1, Y.y1, 'r+'); title(buildKernelDescription(result.k)); xlabel(['bic = ' num2str(result.bic)]); figure
end

function [s] = buildKernelDescription(k)
components = k.components;
operators = k.operators;
s = char(components(end));
for i = length(operators):-1:1
    if strcmp(components(i),'covSum')
        operator = '+';
    else
        operator = 'x';
    end
    s = ['(' char(components(i)) operator s ')'];
end
end

