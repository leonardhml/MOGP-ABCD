function [ hypout ] = expandHyp( hyp )
hypout = hyp.cov;
hypout = [hypout;hyp.smoothing];
hypout = [hypout;hyp.noise];

end

