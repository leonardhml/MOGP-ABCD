function [ model ] = buildMOGPModel(X,Y, cov_options)
% Builds and learns a MOGP model for two outputs.
% Contains:
%   C - the output covariance matrix.
    
    % Fitting and Training
    [C,C11,C12,C21,C22,cov_eval] = buildCrossCovarianceMatrix(X,Y, cov_options);
    model.C = C;
    model.C11 = C11;
    model.C12 = C12;
    model.C21 = C21;
    model.C22 = C22;
    
    [L, p] = chol(C);    
    if p == 0                   % If PD
        alpha = L\(L'\[Y.y1;Y.y2]);
    else                        % Not PD
        alpha = C\[Y.y1;Y.y2];
    end
    
    model.alpha = alpha;
    
    model.X = X;
    model.Y = Y;
    model.cov_options = cov_options;
    
    function [ymeanpred, yvarpred] = predict(xpred, output)
        if output == 1
            g = cov_options.g1;
        elseif output == 2
            g = cov_options.g2;
        else
            error('Model currently only supports two outputs');
        end
        k1 = cov_eval(cov_options.k, cov_options.hyp.cov, g, cov_options.g1, xpred, X.x1);     % Covariance between test points and training points of output 1
        k2 = cov_eval(cov_options.k, cov_options.hyp.cov, g, cov_options.g2, xpred, X.x2);     % Covariance between test points and training points of output 2
        kstar = cov_eval(cov_options.k, cov_options.hyp.cov, g, g, xpred, xpred);
        ks = [k1 k2]';
        
        ymeanpred = ks'*alpha;
        if p == 0                                   % If PD
            v = L'\ks;
            yvarpred = kstar - v'*v;
        else                                        % Not PD
            yvarpred = kstar - ks*(C\ks');
        end
        
    end

    model.predict = @predict;
    
        


end

