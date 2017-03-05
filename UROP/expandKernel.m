function [ kernel_set ] = expandKernel( kernel )
%EXPANDKERNEL Expands the given kernel by composing it with all possible
%combinations of base kernels and operators
% 

base_kernels = BaseKernels.base_kernels;
operators = BaseKernels.operators;
components = BaseKernels.components;

if nargin < 1
    for i = 1 : length(base_kernels)
        kernel_set(i).kernel = base_kernels(i);
        kernel_set(i).components = components(i);
        kernel_set(i).operators = {};
    end
else
    kernel_set = struct('kernel',{},'components',{},'operators',{});
    for j = 1 : length(base_kernels)
        for k = 1 : length(operators)
            k_new.kernel = {operators{k}, {base_kernels{j}, kernel.kernel}};
            k_new.components = [components(j), kernel.components];
            k_new.operators = [operators{k}, kernel.operators];
            kernel_set(end+1) = k_new;
        end
    end
end

end

