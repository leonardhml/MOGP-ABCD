function [bool] = isPSD(M)

bool = 0;
        
% Returns true if the input matrix is PSD
if isequal(size(M), size(M')) && all(all(M - M' < 0.00001))
    [~,S,~] = svd(M);
    if (S>=0)
        bool = 1;
        return
    end
end

end

