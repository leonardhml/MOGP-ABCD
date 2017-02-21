classdef MOGP < handle
    %Class representing the MOGP model. Right now, this model is only able
    %to learn for two outputs and for one dimensional data only. The MOGP
    %uses a convolutional framework to derive a cross covariance function
    %describing the relationship between two outputs.
    %   An MOGP model is defined by the following properties:
    %       k:      the kernel of the latent function.
    %       k_components: the individual base kernels making up k IN ORDER.
    %       g1, g2: the smothing kernels
    %       n:      the number of sampling points for covariance function
    %               construction
    %       a,b:    the range of points at which to sample
    
    properties
        k
        k_components
        g1
        g2
        n
        a
        b
        
        % Properties for fitted model
        Cov
        alpha
        X
        Y
        hyp
        cov_eval
    end
    
    methods
        % argin: 
        %   k: struct containing the kernel itself and names of its
        %   components
        %   g1,g2,n,a,b: see description above
        function obj = MOGP(cov_options)
            obj.k = cov_options.k.kernel;
            obj.k_components = cov_options.k.components;
            obj.n = cov_options.n;
            obj.a = cov_options.a;
            obj.b = cov_options.b;
            obj.g1 = cov_options.g1;
            obj.g2 = cov_options.g2;
        end
        
        % hyp should be a struct containing two fields: cov and noise
        function fit(obj, X, Y, hyp)
            % Set up options
            cov_options.k = obj.k;
            cov_options.hyp = hyp;
            cov_options.g1 = obj.g1;
            cov_options.g2 = obj.g2;
            cov_options.n = obj.n;
            cov_options.a = obj.a;
            cov_options.b = obj.b;
            
            % Build cross covariance matrix
            [Cov,cov_eval] = buildCrossCovarianceMatrix(X,Y, cov_options);
            
            % Set object properties
            obj.Cov = Cov;

            [L, p] = chol(Cov.C);    
            if p == 0                   % If PD
                alpha = L\(L'\[Y.y1;Y.y2]);
            else                        % Not PD
                alpha = Cov.C\[Y.y1;Y.y2];
            end
            obj.alpha = alpha;

            obj.X = X;
            obj.Y = Y;
            obj.hyp = hyp;
            obj.cov_eval = cov_eval;
        end
        
        % Optimise for the hyperparameters
        % Finds the MAP estimate of the hyperparameter posteriors. This
        % assumes that the posterior is a Gaussian distribution, so the MAP
        % is simply the mean of each hyperparameter
        function [hyp] = optimise(obj, X, Y)
            [models, ~] = obj.mcmcPosterior(X,Y);
            hyp = mean(models(:,:), 2);
        end
        
        function [ymeanpred, yvarpred] = predict(obj, xpred, output)
            if output == 1
                g = obj.g1;
            elseif output == 2
                g = obj.g2;
            else
                error('Model currently only supports two outputs');
            end
            k1 = obj.cov_eval(obj.k, obj.hyp.cov, g, obj.g1, xpred, obj.X.x1);     % Covariance between test points and training points of output 1
            k2 = obj.cov_eval(obj.k, obj.hyp.cov, g, obj.g2, xpred, obj.X.x2);     % Covariance between test points and training points of output 2
            kstar = obj.cov_eval(obj.k, obj.hyp.cov, g, g, xpred, xpred);
            ks = [k1 k2]';

            ymeanpred = ks'*obj.alpha;
            [L, p] = chol(obj.Cov.C);
            if p == 0                                   % If PD
                v = L'\ks;
                yvarpred = kstar - v'*v;
            else                                        % Not PD
                yvarpred = kstar - ks*(obj.C\ks');
            end
        end
        
        function [models, logP] = mcmcPosterior(obj,X,Y)
            log_like = @(hyp) obj.logLikelihood(hyp,X,Y);
            prior = Prior(obj.k_components);
            log_prior = @(hyp) prior.logPrior(hyp);
            
            % Apply MCMC hammer
            prior_mu = prior.mu;
            prior_sigma = prior.Sigma;
            minit = mvnrnd(prior_mu', prior_sigma, 2*(length(prior_mu)));
            [models,logP] = gwmcmc(minit', {log_prior log_like}, 5000);
            
            % Bin results?
        end
        
        function [ll] = logLikelihood(obj, hyp, X, Y)
            if isstruct(hyp)
                hyp_struct = hyp;
            else
                % hyp is a vector passed in for MAP calculation
                hyp_struct.cov = hyp(1:end-2);
                hyp_struct.noise = hyp(end-1:end);
            end
            obj.fit(X,Y, hyp_struct);
            y = [Y.y1; Y.y2];
            ll = -((1/2)*y'*obj.alpha) - ((1/2)*log(det(obj.Cov.C))) - ((1/2)*log(2*pi)*length(y));
        end
    end
    
end

