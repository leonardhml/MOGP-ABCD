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
        lb
        ub
    end
    
    methods
        function obj = Prior(k_components)
            obj.k_components = k_components;
            obj.buildPrior(1, 1);
        end
        
        function buildPrior(obj, addNoise, addSmoothing)
            mu = [];
            Sigma_diag = [];
            lb = [];
            ub = [];
            
            for i = 1: length(obj.k_components)
                [mu_new, var_new, lb_new, ub_new] = BaseKernels.getPriorFor(obj.k_components(i));
                mu = [mu; mu_new];
                Sigma_diag = [Sigma_diag; var_new]; 
                lb = [lb; lb_new];
                ub = [ub; ub_new];
            end
            
            if addNoise == 1
                [mu_new, var_new, lb_new, ub_new] = BaseKernels.getPriorFor('noise');
                mu = [mu; mu_new];
                Sigma_diag = [Sigma_diag; var_new]; 
                lb = [lb; lb_new];
                ub = [ub; ub_new];
            end
            
            if addSmoothing == 1
                [mu_new, var_new, lb_new, ub_new] = BaseKernels.getPriorFor('smoothing');
                mu = [mu; mu_new];
                Sigma_diag = [Sigma_diag; var_new]; 
                lb = [lb; lb_new];
                ub = [ub; ub_new];
            end
            
            obj.mu = mu;
            obj.Sigma = diag(Sigma_diag);
            obj.lb = lb;
            obj.ub = ub;
        end
        
        function [lp] = logPrior(obj, hyp)
            lp = -(1/2)*(hyp - obj.mu)'*((obj.Sigma)\(hyp - obj.mu)) - (1/2)*log(det(obj.Sigma)) - (length(obj.mu)/2)*log(2*pi);
        end
    end
end

