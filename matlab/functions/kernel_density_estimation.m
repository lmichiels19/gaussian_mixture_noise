function [p_density,p_weights] = kernel_density_estimation(evalf_samples, samples, weights, sig)
%KERNEL_DENSITY_ESTIMATION Implementation of a kernel density estimator

    Np = size(evalf_samples,2); % number of samples in proposed density
    Ns = size(samples, 2); % number of samples in density representation

    % Kernel density
    p_density = evalf_samples;
    p_weights = zeros(1,Np);    
    
    % parallel computing    
    %p_density = parallel.pool.Constant(p_density);
    %samples = parallel.pool.Constant(samples);
    %sig = parallel.pool.Constant(sig);
    
    for i=1:Ns
        p_weights = p_weights + weights(i).*mvnpdf(p_density',samples(:,i)',sig)';
    end
    
    p_weights = p_weights ./ sum(p_weights);
end

