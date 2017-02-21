classdef Prior < handle
    %Hyperparameter priors. Every hyperparameter across all base kernels
    %will have an associated prior, which is based on section 4.2 of
    %Malkomes, 2016.
    %    The prior model is built from the relevant base kernel
    %    hyperparameter priors, each of which is a univariate Gaussian. The
    %    final prior is a multivariate Gaussian. It is important that the
    %    order of the hyperparameters follows exactly according to the
    %    ordering of the hyperparameters for the composite latent kernel.
    %    After the prior is built, the logPrior can be calculated.
    
    properties
        k_components
        mu
        Sigma
    end
    
    methods
        function obj = Prior(k_components)
            obj.k_components = k_components;
            obj.buildPrior();
            obj.addNoise();
        end
        
        function buildPrior(obj)
            mu = [];
            Sigma_diag = [];
            
            for i = 1: length(obj.k_components)
                [mu_new, var_new] = Prior.getHypPrior(obj.k_components(i));
                mu = [mu; mu_new];
                Sigma_diag = [Sigma_diag; var_new]; 
            end
            
            obj.mu = mu;
            obj.Sigma = diag(Sigma_diag);
        end
        
        function [lp] = logPrior(obj, hyp)
            lp = -(1/2)*(hyp - obj.mu)'*inv(obj.Sigma)*(hyp - obj.mu) - (1/2)*log(det(obj.Sigma)) - (length(obj.mu)/2)*log(2*pi);
        end
        
        % Noise priors for each of the two outputs
        % Independent noise hyperparameter on a log scale for two outputs
        function addNoise(obj)
            obj.mu = [obj.mu; -0.1; -0.1];
            obj.Sigma = [obj.Sigma, zeros(length(obj.Sigma),2); zeros(2,length(obj.Sigma)), diag([0.3^2; 0.3^2])];
        end
    end
    
    methods (Static)
        function [mu, var] = getHypPrior(component)
            
            % Define hyp priors here
            % Base kernels
            % covSEiso: k(x,z) = sf^2 * exp(-(x-z)^2/(2*ell^2))                         hyp = [log(ell); log(sf)]
            % covLINscaleshift: k(x,z) = xz/ell^2                                       hyp = [log(ell); shift] 
            % covPeriodic: k(x,z) = sf^2 * exp( -2*sin^2( pi*(x-z)/p )/ell^2 )          hyp = [log(ell); log(p) ; log(sf)]
            % cpvRQiso: k(x,z) = sf^2 * [1 + (x-z)'*inv(P)*(x-z)/(2*alpha)]^(-alpha)    hyp = [log(ell); log(sf); log(alpha)]

            % (NOT FOUND) covPRiso: k(x,z) = sf^2 * [1 + (x-z)^2/(2*alpha*ell^2)]^(-alpha)          hyp = [log(ell), log(sf), log(alpha)]
            % covConst: k(x,z) = sf^2                                                   hyp = [log(sf)]
            % covWhiteNoise?
            persistent hyp_covSEiso;
            hyp_covSEiso = struct('mu', [0.1; 0.4], 'var', [0.7^2; 0.7^2]);

            persistent hyp_covPeriodic;
            hyp_covPeriodic = struct('mu', [2; 0.1 + log(pi); 0.4], 'var', [0.7^2; 0.7^2; 0.7^2]);

            persistent hyp_covLin;
            hyp_covLin = struct('mu', [-0.8; 0], 'var',  [0.7^2; 2^2]);
            
            
            persistent hyp_covRQiso;
            hyp_covRQiso = struct('mu', [0.1; 0.4; 0.05], 'var', [0.7^2; 0.7^2; 0.7^2]);
            
            if (strcmp(component, 'SE'))
                mu = hyp_covSEiso.mu;
                var = hyp_covSEiso.var;
            elseif (strcmp(component, 'LIN'))
                mu = hyp_covLin.mu;
                var = hyp_covLin.var;
            elseif (strcmp(component, 'PER'))
                mu = hyp_covPeriodic.mu;
                var = hyp_covPeriodic.var;
            elseif (strcmp(component, 'RQ'))
                mu = hyp_covRQiso.mu;
                var = hyp_covRQiso.var;
            else
                Error('Base Kernel not supported');
            end
        end
    end
    
end

