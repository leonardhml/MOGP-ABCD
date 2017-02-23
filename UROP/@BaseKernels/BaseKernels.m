classdef BaseKernels
    %BASEKERNELS Describes all the base kernels used in this UROP project

    %   This class contains information about all the base kernels used,
    %   including the kernels themselves, their hyperparameter priors,
    %   component names, and also operations used to compose kernels.
    
    properties

    end
    
    methods (Static)
        % ARGIN: component is the kernel component name e.g. 'SE'
        %        Also allows for non-kernel components, like noise or
        %        smoothing
        function [mu, var] = getPriorFor(component)
            if (strcmp(component, 'SE'))
                mu = BaseKernels.hyp_covSEiso.mu;
                var = BaseKernels.hyp_covSEiso.var;
            elseif (strcmp(component, 'LIN'))
                mu = BaseKernels.hyp_covLin.mu;
                var = BaseKernels.hyp_covLin.var;
            elseif (strcmp(component, 'PER'))
                mu = BaseKernels.hyp_covPeriodic.mu;
                var = BaseKernels.hyp_covPeriodic.var;
            elseif (strcmp(component, 'RQ'))
                mu = BaseKernels.hyp_covRQiso.mu;
                var = BaseKernels.hyp_covRQiso.var;
            elseif (strcmp(component, 'noise'))
                mu = BaseKernels.hyp_noise.mu;
                var = BaseKernels.hyp_noise.var;
            elseif (strcmp(component, 'smoothing'))
                mu = BaseKernels.hyp_smoothing.mu;
                var = BaseKernels.hyp_smoothing.var;
            else
                Error('Base Kernel not supported');
            end
        end
        
        function [s] = base_kernels()
            s = {'covSEiso', 'covLINscaleshift', 'covPeriodic', 'covRQiso'} ;
        end
        
        function [s] = components()
            s = {'SE', 'LIN', 'PER', 'RQ'};
        end
        
        function [s] = operators()
            s = {'covProd', 'covSum'};
        end

        % Define hyp priors here
        % Base kernels
        % covSEiso: k(x,z) = sf^2 * exp(-(x-z)^2/(2*ell^2))                         hyp = [log(ell); log(sf)]
        % covLINscaleshift: k(x,z) = xz/ell^2                                       hyp = [log(ell); shift] 
        % covPeriodic: k(x,z) = sf^2 * exp( -2*sin^2( pi*(x-z)/p )/ell^2 )          hyp = [log(ell); log(p) ; log(sf)]
        % cpvRQiso: k(x,z) = sf^2 * [1 + (x-z)'*inv(P)*(x-z)/(2*alpha)]^(-alpha)    hyp = [log(ell); log(sf); log(alpha)]

        % (NOT FOUND) covPRiso: k(x,z) = sf^2 * [1 + (x-z)^2/(2*alpha*ell^2)]^(-alpha)          hyp = [log(ell), log(sf), log(alpha)]
        % covConst: k(x,z) = sf^2                                                   hyp = [log(sf)]
        % covWhiteNoise?
        function [s] = hyp_covSEiso()
            s = struct('mu', [0.1; 0.4], 'var', [0.7^2; 0.7^2]);
        end
        function [s] = hyp_covPeriodic()
            s = struct('mu', [2; 0.1 + log(pi); 0.4], 'var', [0.7^2; 0.7^2; 0.7^2]);
        end
        function [s] = hyp_covLin()
            s = struct('mu', [-0.8; 0], 'var',  [0.7^2; 2^2]);
        end
        function [s] = hyp_covRQiso()
            s = struct('mu', [0.1; 0.4; 0.05], 'var', [0.7^2; 0.7^2; 0.7^2]);
        end
        function [s] = hyp_noise()
            s = struct('mu', [-0.1; -0.1], 'var', [0.3^2; 0.3^2]);
        end
        function [s] = hyp_smoothing()
            s = struct('mu', [0.1; 0.1; 0], 'var', [0.1; 0.1; 0.5^2]);
        end

    end
    
end

