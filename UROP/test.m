  x = gpml_randn(0.8, 19, 1);                 % 20 training inputs
  y = sin(3*x) + 0.1*gpml_randn(0.9, 19, 1);  % 20 noisy training targets
  xs = linspace(-3, 3, 40)';                  % 61 test inputs 
  
  % Test: add extra training points
  xextra = [0];
  yextra = sin(3*xextra) + 0.1*gpml_randn(0.9, length(xextra), 1);
  x = [x;xextra];
  y = [y; yextra];
  meanfunc = [];
  covfunc = @covPeriodic;
  likfunc = @likGauss;
  
  hyp = struct('mean', [], 'cov', [0.5 0.5 0.5], 'lik', 0);    % struct encapsulating all possibly hyperparams
  
  hyp2 = minimize(hyp, @gp, -500, @infGaussLik, meanfunc, covfunc, likfunc, x,y);   % optimize hyperparams
  hyp2
  % Find best values for hyperparams - minimise negative log marginal likelihood
  % For ABCD - selection between models is done using BIC, an extension of
  % nlml: BIC(M) = -2*nlml + |M|logn where |M| is the number of kernel
  % parameters (hyperparams) and n is the number of data points
  
  
  
  [mu s2] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, x, y, xs);   % make predictions
  
  ys = sin(3*xs) + 0.1*gpml_randn(0.9,61,1)
  
  f = [mu+2*sqrt(s2); flipdim(mu-2*sqrt(s2),1)];
  fill([xs; flipdim(xs,1)], f, [7 7 7]/8)
  hold on; plot(xs, mu, 'b'); plot(x, y, 'r+'); plot (xs, ys, 'gx')