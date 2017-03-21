function [ Cov, cov_eval] = buildCrossCovarianceMatrix(X,Y,cov_options)
% Create the covariance matrix for the combined output.
% For the 2 output case, the resulting covariance matrix will have size
% [(nx+ny) *(nx+ny)] where nx and ny are the sizes of the observations for
% the two outputs respectively.
%
% ARGs:
% X is a struct containing training inputs x1 and x2
% Y is similarly defined
% cov_options is a struct containing the options for the covariance matrix.
%   k:      The kernel of the latent process
%   hyp:    hyp.cov is a vector of hyperparameters for k.
%           hyp.smoothing is a vector of hyperparameters for smoothing
%           kernels.
%           hyp.noise is the hyperparameters for the independent noise
%           term.
%   g1, g2:  The smoothing kernels applied for output 1 and 2 respectively
%   n:      The accuracy term for the convolution summation function
%   a, b:      The range of the (normalised) inputs
%
    [y1l, y1w] = size(Y.y1);
    [y2l, y2w] = size(Y.y2);
    [x1l, x1w] = size(X.x1);
    [x2l, x2w] = size(X.x2);
    assert(y1l == x1l && y2l == x2l && y1w == x1w && y2w == x2w,'Dimensions of x and y must be the same for both outputs');
    
    cov_eval = convSumCov(cov_options.n,cov_options.a,cov_options.b);
    % Create covariance matrix and add independent noise
    % TODO: Add independent Gaussian term
    k = cov_options.k;
    hyp = cov_options.hyp;
    g1.g = cov_options.g1;
    g1.hyp = cov_options.hyp.smoothing(1);
    g2.g = cov_options.g2;
    g2.hyp = cov_options.hyp.smoothing(2);
    x1 = X.x1;
    x2 = X.x2;
    C11 = cov_eval(k, hyp, g1, g1, x1, x1) + (exp(2*hyp.noise(1)) * eye(length(x1)));
    C12 = cov_eval(k, hyp, g1, g2, x1, x2);
    C21 = cov_eval(k, hyp, g2, g1, x2, x1);
    C22 = cov_eval(k, hyp, g2, g2, x2, x2) + (exp(2*hyp.noise(2)) * eye(length(x2)));
    C = [C11 C12; C21 C22];
    Cov.C = C/2 + C'/2;
    Cov.C11 = C11;
    Cov.C12 = C12;
    Cov.C21 = C21;
    Cov.C22 = C22;
    
    [R, p] = chol(Cov.C);
    if p ==0 
        max(max(Cov.C))
        min(min(Cov.C))
        cov_options.hyp
    else 
        max(max(Cov.C))
        min(min(Cov.C))
        cov_options.hyp
    end
%     if (isPSD(Cov.C))
%         
%     else
%         error('Covariance matrix is not PSD!!!');
%     end
%     
end

