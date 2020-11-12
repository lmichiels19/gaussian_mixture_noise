%*********************************************
% do monte-carlo runs of tracker
%*********************************************

addpath(genpath('../'));

doPlot = true;

Nc = 501;
limits = [-5 5];
max_run_count = 20;

gauss_mean = [-1.5,1.5];
gauss_sig = diag(0.25);

noise_sig = diag(0.1);

meas_mean = [-2.5 2.5];
meas_sig = diag(0.1);

gauss_count = size(gauss_mean, 2);
Ds = size(gauss_mean, 1);

sig_over_runs = zeros(max_run_count,1);

for run_count=1:max_run_count

% samples
g_samples = linspace(limits(1,1), limits(1,2), Nc);
g_weights = zeros(1,Nc);
for i=1:gauss_count
   g_weights = g_weights + mvnpdf(g_samples', gauss_mean(:,i)', gauss_sig')';
end
g_weights = g_weights ./ sum(g_weights);
g_pf_weights = g_weights;

% noise
n_weights = mvnpdf(g_samples', zeros(Ds,1)', noise_sig');
n_weights = n_weights ./ sum(n_weights);

gt_weights = g_weights;

for runs=1:run_count

    % noise
    g_pf_weights = conv(g_pf_weights,n_weights,'same');   
    
    %likelihood
    g_likelihood = zeros(size(g_samples));
    li_count = size(meas_mean,2);
    for i=1:li_count
       g_likelihood = g_likelihood + (1/li_count).*mvnpdf(g_samples', meas_mean(:,i)', meas_sig')';
    end
    
    % noise
    g_weights = conv(g_weights,n_weights,'same');        
    % update
    g_weights = g_weights .* g_likelihood;
    g_weights = g_weights ./ sum(g_weights);
    
    % noise
    gt_weights = conv(gt_weights,n_weights,'same');        
    % update
    gt_weights = gt_weights .* g_likelihood.^2;
    gt_weights = gt_weights ./ sum(gt_weights);
end

fuse_sig = diag(0.01);
fuse_sig_pf = diag( noise_sig.^4.*run_count.^4 );

[~,g_pf_weights] = kernel_density_estimation(g_samples, g_samples, g_pf_weights, fuse_sig_pf);

g_fused_weights = g_weights;
f_goal = @(x)(1000*kernel_sig_cost_function(g_samples,gt_weights, g_fused_weights, g_pf_weights, x));

[fuse_sig,fval,exitflag,output] = fminbnd(f_goal,0,1);
fuse_sig = fuse_sig
fval = fval

[~,g_fused_weights] = kernel_density_estimation(g_samples, g_samples, g_fused_weights, fuse_sig);

g_fused_weights = (g_fused_weights.*g_fused_weights);
g_fused_weights = g_fused_weights./g_pf_weights;    
g_fused_weights = g_fused_weights ./ sum(g_fused_weights);

sig_over_runs(run_count) = fuse_sig;

end

% plot
if doPlot
       


    
    figure(1); clf; hold on;
    plot(g_samples, g_weights, 'Color', [0 0 1], 'LineWidth', 1.2);
    plot(g_samples, n_weights, 'Color', [0 1 0], 'LineWidth', 1.2);
    plot(g_samples, gt_weights, 'Color', [1 0 1], 'LineWidth', 1.2);
    plot(g_samples, g_fused_weights, 'Color', [0 0 0], 'LineWidth', 1.2);
    plot(g_samples, g_pf_weights, 'Color', [0 1 1], 'LineWidth', 1.2);
    
    legend("Original", "Noise", "GT", "Fused");
    
end
