%1. Build k (with covSum, covProd)
%   a. Range: [0, a]^2 where a is the normalised range of input points
%2. Choose a g1, g2 (smoothing kernels) (hyps?)
%   a. Range: [0, a] for both
%2a. Hyperparameter optimisation: Apply BO algorithm, sampling from
%    (assumed) hyperparameter priors
%    For smoothing kernels, signal variance can be negative
%    Length scale depends on normalised input (degree of
%    interconnectedness) e.g. 1/n, 2/n,...
%3. Build X-cov via summation approximation
%   a. Finding C_ij(x,y) = int(g_i(x-z)g_j(y-z')k(z,z')dzdz'
%                        = conv2(G(x,y),k)
%                        = h_xh_y/9 sum(1...N_x)sum(1...N_y)c_ijgi(x-xi')gj(y-y_j')k(x_i',y_j')
%   b. For 2 outputs, need to find  C_11, C_12, C_21, C_22
%   c. Take note: For discretisation,
%       1. Get samples in ranges N_x = [1,N = 2^n - 1], N_y = [1, M = 2^m - 1]
%          Choices of n and m will decide the accuracy of the model
%       2. Transform k to kbar: Introduce constant c_ij (based on which
%          sample of points x_i, y_j)
%       3. Find grid sizes h_x and h_y
%4. Use GP and X-cov to predict

%%%
% Sample training and test inputs
%%%
[X,Y,xpred,ypred,outputpred] = generateData();

% Base kernels
% See Prior.m for more details
kernel_set = {'covSEiso', 'covLINscaleshift', 'covPeriodic', 'covRQiso'} ;
base_kernels = {'SE', 'LIN', 'PER', 'RQ'}; 

% A sample k
k.kernel = {'covProd', {'covSEiso', 'covPeriodic'}} ;
k.components = {'SE', 'PER'};

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

tic 

% Find MLE estimate for hyp
hyp_opt = model.optimise(X,Y);
hyp.cov = hyp_opt(1:end-4);
hyp.smoothing = hyp_opt(end-3:end-2);
hyp.noise = hyp_opt(end-1:end);

% Fit and predict
model.fit(X,Y,hyp);
[mu, s2] = model.predict(xpred, 1);
fprintf('Model Evidence is %d.\n', model.modelEvidence);

toc

% Visualise
f = [mu+2*sqrt(diag(s2)); flipdim(mu-2*sqrt(diag(s2)),1)];
fill([xpred; flipdim(xpred,1)], f, [7 7 7]/8)
hold on; plot(xpred, mu, 'b'); plot(X.x1, Y.y1, 'r+');

% fileID = fopen('hyp.txt','a');
% fprintf(fileID,'%f\n', hyp_opt);
% fprintf(fileID,'\n');
% fclose(fileID);