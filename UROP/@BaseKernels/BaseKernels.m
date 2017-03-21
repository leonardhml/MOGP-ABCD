classdef BaseKernels
    %BASEKERNELS Describes all the base kernels used in this UROP project

    %   This class contains information about all the base kernels used,
    %   including the kernels themselves, their hyperparameter priors,
    %   component names, and also operations used to compose kernels.
    %   Also includes lower and upper bounds for the values of each hyp.
    %   This is usually used in optimisation by SA or PSO.
    
    properties (Constant)
        base_kernels = {'covSEiso', 'covLINscaleshift', 'covPeriodic', 'covRQiso'} ;
        components = {'SE', 'LIN', 'PER', 'RQ'}; 
        operators = {'covProd', 'covSum'};

        % Define hyp priors here
        % Base kernels
        % covSEiso: k(x,z) = sf^2 * exp(-(x-z)^2/(2*ell^2))                         hyp = [log(ell); log(sf)]
        % covLINscaleshift: k(x,z) = xz/ell^2                                       hyp = [log(ell); shift] 
        % covPeriodic: k(x,z) = sf^2 * exp( -2*sin^2( pi*(x-z)/p )/ell^2 )          hyp = [log(ell); log(p) ; log(sf)]
        % cpvRQiso: k(x,z) = sf^2 * [1 + (x-z)'*inv(P)*(x-z)/(2*alpha)]^(-alpha)    hyp = [log(ell); log(sf); log(alpha)]

        % (NOT FOUND) covPRiso: k(x,z) = sf^2 * [1 + (x-z)^2/(2*alpha*ell^2)]^(-alpha)          hyp = [log(ell), log(sf), log(alpha)]
        % covConst: k(x,z) = sf^2                                                   hyp = [log(sf)]
        % covWhiteNoise?

        hyp_covSEiso = struct('mu', [0.1; 0.4], 'var', [0.7^2; 0.7^2], 'lb', [-5;-5], 'ub', [5;5]);
        hyp_covPeriodic = struct('mu', [2; 0.1 + log(pi); 0.4], 'var', [0.7^2; 0.7^2; 0.7^2], 'lb', [-5;-5;-5;], 'ub', [5;5;5]);
        hyp_covLin = struct('mu', [0.1; 0], 'var',  [0.7^2; 2^2], 'lb', [-5;-5], 'ub', [5;5]);
        hyp_covRQiso = struct('mu', [0.1; 0.4; 0.05], 'var', [0.7^2; 0.7^2; 0.7^2], 'lb', [-5;-5;-5], 'ub', [5;5;5]);
        hyp_noise = struct('mu', [0.1; 0.1], 'var', [1^2; 1^2], 'lb', [-5;-5], 'ub', [5;5]);
        hyp_smoothing = struct('mu', [0.1; 0.1], 'var', [0.01; 0.01], 'lb', [-5;-5], 'ub', [5;5]);
    end
    
    methods (Static)
        % ARGIN: component is the kernel component name e.g. 'SE'
        %        Also allows for non-kernel components, like noise or
        %        smoothing
        function [mu, var, lb, ub] = getPriorFor(component)
            kernel = [];
            if (strcmp(component, 'SE'))
                kernel = BaseKernels.hyp_covSEiso;
            elseif (strcmp(component, 'LIN'))
                kernel = BaseKernels.hyp_covLin;
            elseif (strcmp(component, 'PER'))
                kernel = BaseKernels.hyp_covPeriodic;
            elseif (strcmp(component, 'RQ'))
                kernel = BaseKernels.hyp_covRQiso;
            elseif (strcmp(component, 'noise'))
                kernel = BaseKernels.hyp_noise;
            elseif (strcmp(component, 'smoothing'))
                kernel = BaseKernels.hyp_smoothing;
            else
                Error('Base Kernel not supported');
            end
            
            mu = kernel.mu;
            var = kernel.var;
            lb = kernel.lb;
            ub = kernel.ub;
        end
    end
    
end

