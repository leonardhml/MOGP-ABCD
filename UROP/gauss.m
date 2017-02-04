function [y] = gauss(x,std,offset)
    if nargin < 1
        error('x,std is a required input');
    end
    
    if nargin < 2
        std = 1.0;
    end
    
    if nargin < 3
        offset = zeros(length(x), 1);
    end
    x = x - offset;
    y = exp(-x.^2/(2*std^2)) / (std*sqrt(2*pi));
end