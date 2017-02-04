cov = {'covSum',{cgi,cpe}}; hyp = [hypgi; hyppe];      % sum
cov = {@covProd,{cgi,cpe}};   hyp = [hypgi; hyppe];    % product
% See usageCov.m for list of covariance matrices

%1. Build k (with covSum, covProd) (what about hyps?)
%   a. Range: [0, a], [0, b] (how to decide?)
%2. Choose a g1, g2 (smoothing kernels) (hyps?)
%   a. Range: [0, a], [0, b] (how to decide?)
%3. Build X-cov via 2d convolution
%   a. Finding C_ij(x,y) = int(g_i(x-z)g_j(y-z')k(z,z')dzdz'
%                         = conv2(G(x,y),k)
%   b. For 2 outputs, need to find  C_11, C_12, C_21, C_22
%   c. Take note: For discretisation,
%       1. Get samples in ranges N_x = [1,N = 2^n - 1], N_y = [1, M = 2^m - 1] (How do
%          we decide?)
%          Introduce N+1th and M+1th point such that k(x_N+1, y_M+1) = 0
%       2. Transform k to kbar: Introduce constant c_ij (based on which
%          sample of points x_i, y_j)
%       3. Find grid sizes h_x and h_x'
%       4. Discrete convolution summation is then:
%          I(x_i,y_j) = (h_x)(h_x')/9 * sumi'(1...N+1)sumj'(1...M+1) g1(x_i - x_i')g2(y_i -y_i')kbar(x_i',y_i')
%       5. Replacement with cyclic summation: i = [1, 2N_x+2], j = [1, 2N_y+2]
%       6. Apply FFT
%          a. The resulting matrix must be restricted to the domain [1,N+1]
%             and [1,M+1]
%4. Use GP and X-cov to predict
%   a. Prediction points must be within sample range (?) and must be part
%      of grid points (?)

% Base kernels
% covSEiso: k(x,z) = sf^2 * exp(-(x-z)^2/(2*ell^2)) hyp = [log(ell); log(sf)]
% covLINiso: k(x,z) = xz/ell^2 hyp = [log(ell)] (missing l - see lloyd)
base_kernels = {'covSEiso', 'covLIN'
% Sample k
  k = {'covSum', {'covLIN', 'covPER'}}
