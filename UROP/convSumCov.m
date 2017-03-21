function [c] = convSumCov(n, a, b)
% Calculate the (cross) covariance between two outputs at a certain point
% by solving a convolution summation. The formula is given as:
% C_ij(x,y) = h_xh_y/9 sum(1...N_x)sum(1...N_y)c_ijgi(x-x_i')gj(y-y_j')k(x_i',y_j')
% where:
% N_x = N_y = 2^n - 1 is the number of samples
% h_x = h_y = (b - a)/(N_x - 1) is the grid size (separation between samples)
% c_ij = 
%       16 if both i and j are even
%       4 if both i and j are odd
%       8 if i is odd and j is even, or i is even and j is odd
%       4 for i = 1, len and even j's
%       4 for j = 1, len and even i's
%       2 for i = 1, len and odd j's
%       2 for j = 1, len and odd i's
%       1 for (i,j) = (1,1),(len,1),(1,len),(len,len)
%
% We select 2^n - 1 samples from k, g1 and g2 each. The higher the n, the
% more accurate the resulting value of c. These samples are picked from the
% range [a,b] so that the i-th sample is at position a+(i-1)h.
%
% TODO: Need to figure out how to select [a, b]
%
% Expanding beyond 1D inputs: introduce more summations. Treat each
% additional dimension as an additional variable.
% 
% ARGs:
% n: accuracy term. The higher the value, the more samples we pick and the
% more accurate the covariance function.
% [a,b]: range over the convolution
%
% Returns a covariance function
%
% Build C
n_samples = 2^n - 1;
C = ones(n_samples, n_samples);
C(1:2:n_samples, 1:2:n_samples) = 4;
C(2:2:n_samples, 2:2:n_samples) = 16;
C(1:2:n_samples, 2:2:n_samples) = 8;
C(2:2:n_samples, 1:2:n_samples) = 8; 
C(1, 2:2:n_samples) = 4; 
C(n_samples, 2:2:n_samples) = 4;
C(2:2:n_samples, 1) = 4; 
C(2:2:n_samples, n_samples) = 4; 
C(1, 1:2:n_samples) = 2; 
C(n_samples, 1:2:n_samples) = 2;
C(1:2:n_samples, 1) = 2;
C(1:2:n_samples, n_samples) = 2; 
C(1,1) = 1;
C(1,n_samples) = 1;
C(n_samples,1) = 1;
C(n_samples,n_samples) = 1;

    % Function to build the entire covariance matrix between x and y
    % hyp is a struct containing hyperparameters for cov, smoothing, and
    % noise
    % g1 and g2 are structs containing g, the kernel itself, and hyp, the
    % signal variance of g
    function [cov] = evalCov(k, hyp, g1, g2, x, y)
        
    h1 = (b-a + 4*g1.hyp)/(n_samples-1);
    samples1 = [a - 2*g1.hyp : h1 : b + 2*g1.hyp]';
    h2 = (b-a + 4*g2.hyp)/(n_samples-1);
    samples2 = [a - 2*g2.hyp : h2 : b + 2*g2.hyp]';
        % Build k
        K = feval(k{:}, hyp.cov, samples1, samples2);
%         % Evaluate covariance at two points
%         function [v] = eval(x,y)  
%             % Build g's
%             xs = ones(n_samples, 1) * x;
%             ys = ones(n_samples, 1) * y;
%             G1 = g1(xs - samples);
%             G2 = g2(ys - samples);
% 
%             v = (h^2 / 9) * G1' * (K.*C) * G2;
%         end
% 
%         % Compute covariance matrix
%         cov = ones(length(x), length(y));
%         for io=[1:1:length(x)]
%             for jo=[1:1:length(y)]
%                 cov(io,jo) = eval(x(io),y(jo));
%             end
%         end
        
        xs = repmat(x, 1, n_samples) - repmat(samples1', length(x), 1);
        ys = repmat(y, 1, n_samples) - repmat(samples2', length(y), 1);
        G1 = g1.g(xs, g1.hyp);
        G2 = g2.g(ys, g2.hyp);
        cov = (h1*h2 / 9).* G1 * (K.*C) * G2';
    end

c = @evalCov;
end

