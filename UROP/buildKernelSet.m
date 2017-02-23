function [kernel_set] = buildKernelSet(depth)
% Builds a set of all kernels over a fixed grammar and fixed set of base
% kernels.
%
%   ARGIN: depth is the maximum depth of the kernel search tree. All
%   kernels up to and including kernels of this depth will be added into
%   the kernel set.
%
%   ARGOUT: 
%       kernel_set is the set of kernel structs k. Each struct contains
%       three components: k.kernel is the kernel itself. k.components is a 
%       vector describing in order the base kernels that make up k.kernel.
%       k.operators describe the order of operators (+ and x) that compose
%       the kernel.
%   Example:
%       k.kernel = {'covProd', {'covSEiso', 'covPeriodic'}} ;
%       k.components = {'SE'; 'PER'};
%       k.operators = {'x'};


base_kernels = BaseKernels.base_kernels;
components = BaseKernels.components;
operators = BaseKernels.operators;

for i = 1 : length(base_kernels)
    kernel_set(i).kernel = base_kernels(i);
    kernel_set(i).components = components(i);
    kernel_set(i).operators = {};
    frontier(i).kernel = base_kernels{i};
    frontier(i).components = components(i);
    frontier(i).operators = {};
end

depth = depth - 1;

while depth > 0
    frontier2 = struct('kernel',{},'components',{},'operators',{});
    for i = 1 : length(frontier)
        for j = 1 : length(base_kernels)
            for k = 1 : length(operators)
                k_new.kernel = {operators{k}, {base_kernels{j}, frontier(i).kernel}};
                k_new.components = [components(j), frontier(i).components];
                k_new.operators = [operators{k}, frontier(i).operators];
                kernel_set(end+1) = k_new;
                frontier2(end+1) = k_new;
            end
        end
    end
    frontier = frontier2;
    depth = depth - 1;
end
end