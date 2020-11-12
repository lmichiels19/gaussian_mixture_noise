function [cost] = kernel_sig_cost_function(g_samples,gt_weights, g_fused_weights, g_pf_weights, fuse_sig)
%KERNEL_SIG_COST_FUNCTION Summary of this function goes here
%   Detailed explanation goes here


    if size(fuse_sig,1) > 1
        Ns = size(g_samples,2); % number of samples in proposed density

        % Kernel density
        p_weights = zeros(1,Ns);   
        for i=1:Ns
            p_weights = p_weights + g_fused_weights(i).*mvnpdf(g_samples',g_samples(:,i)',fuse_sig(i))';
        end

        g_fused_weights = p_weights ./ sum(p_weights);
    else
    	[~,g_fused_weights] = kernel_density_estimation(g_samples, g_samples, g_fused_weights, fuse_sig);
    end

    g_fused_weights = (g_fused_weights.*g_fused_weights);
    g_fused_weights = g_fused_weights./g_pf_weights;    
    g_fused_weights = g_fused_weights ./ sum(g_fused_weights);

    cost = mean((g_fused_weights-gt_weights).^2);
end

