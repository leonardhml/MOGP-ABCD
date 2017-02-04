function [ output_args ] = kTransform(cov, hyp, x)

% Transforms the given covariance function k(x,x') to a cyclic definition.
% In the context of discrete convolution summations, this cyclic
% convariance function can then be used in an FFT based method for more
% efficient convolutions.
%
% Note that len = length(x) = 2^n - 1 for some n.
%
% The pipeline is as follows:
% 1. kbar = c_ij * k(x_i, x_j), where:
%       c_ij = 16 if both i and j are even
%       c_ij = 4 if both i and j are odd
%       c_ij = 8 if i is odd and j is even, or i is even and j is odd
%       c_ij = 4 for i = 1, len and even j's
%       c_ij = 4 for j = 1, len and even i's
%       c_ij = 2 for i = 1, len and odd j's
%       c_ij = 2 for j = 1, len and odd i's
%       c_ij = 1 for (i,j) = (1,1),(len,1),(1,len),(len,len)
% 2. Let K be the (len * len) matrix from evaluating kbar. Append a zero
%    row and column to K. K is now a (len + 1 * len + 1) matrix.
% 3. Pad matrix with zeroes for [len + 1, 2*len]

end

