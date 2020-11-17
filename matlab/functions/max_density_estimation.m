function [max_density] = max_density_estimation(samples, weights, bandwidth)
%DENSITYESTIMATION Estimate density of given Dirac-Mixtures

    Ns = size(samples, 2); % number of samples in density representation

    % Kernel density
    density = zeros(1,Ns);    
    
    % parallel computing    
    %p_density = parallel.pool.Constant(p_density);
    %samples = parallel.pool.Constant(samples);
    %sig = parallel.pool.Constant(sig);
    
    for i=1:Ns
        density =  density + weights(i).*mvnpdf(samples',samples(:,i)',diag(bandwidth))';
    end
    
    max_density = max(density);%/mean(density);

end

