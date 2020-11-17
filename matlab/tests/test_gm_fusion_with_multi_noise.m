%*********************************************
% do monte-carlo runs of tracker
%*********************************************

addpath(genpath('../'));

doPlot = true;

Nc = 101;
limits = [-5 5; -5 5];
max_run_count = 10;

gauss_mean = [-1.5,2.5; -1.5,2.5];
gauss_sig = diag([0.25 0.25]);

noise_sig = diag([0.4 0.4]);

meas_mean_0 = 2.*randn(2,2);
%meas_mean_0 = [-2.5 2.5; -2.5 2.5];
meas_sig = diag([0.2 0.2]);

gauss_count = size(gauss_mean, 2);
Ds = size(gauss_mean, 1);

sig_over_runs = zeros(max_run_count,1);

for run_count=max_run_count

% samples
g_samples = create_sample_grid(limits, repmat(Nc,Ds,1));%linspace(limits(1,1), limits(1,2), Nc);
Np = size(g_samples,2);
g_weights_org = zeros(1,Np);
for i=1:gauss_count
   g_weights_org = g_weights_org + mvnpdf(g_samples', gauss_mean(:,i)', gauss_sig')';
end
g_weights_org = g_weights_org ./ sum(g_weights_org);
g_weights_1 = g_weights_org;
g_weights_2 = g_weights_org;
g_pf_weights_1 = g_weights_1;
%g_pf_weights_1 = ones(size(g_weights_1));
g_pf_weights_2 = g_weights_2;
%g_pf_weights_2 = ones(size(g_weights_1));

% noise
n_weights = mvnpdf(g_samples', zeros(Ds,1)', noise_sig')';
n_weights = n_weights ./ sum(n_weights);

gt_weights = g_weights_1;

g_likelihood_1 = ones(size(g_weights_1));
g_likelihood_2 = ones(size(g_weights_2));

meas_mean = meas_mean_0;

for runs=1:run_count

    gamma = 1;
    % noise
    if( runs > 1 )
        g_pf_weights_1 = g_pf_weights_1.*conv(g_likelihood_1,n_weights.^gamma,'same');
    end
    g_pf_weights_1 = conv(g_pf_weights_1,n_weights,'same');

    % noise
    if( runs > 1 )
        g_pf_weights_2 = g_pf_weights_2.*conv(g_likelihood_2,n_weights.^gamma,'same');
    end
    g_pf_weights_2 = conv(g_pf_weights_2,n_weights,'same');

    
    %likelihood
    g_likelihood_1 = zeros(size(g_weights_1));
    g_likelihood_2 = zeros(size(g_weights_2));
    li_count = size(meas_mean,2);
    mean_1 = mvnrnd(meas_mean', meas_sig);
    mean_2 = mvnrnd(meas_mean', meas_sig);
    for i=1:li_count %(1/li_count).*
       g_likelihood_1 = g_likelihood_1 + mvnpdf(g_samples', mean_1(i,:), meas_sig')'.*sqrt(2*pi*det(meas_sig));
       g_likelihood_2 = g_likelihood_2 + mvnpdf(g_samples', mean_2(i,:), meas_sig')'.*sqrt(2*pi*det(meas_sig));
       %meas_mean(:,i) = meas_mean(:,i) + 0.1;
    end
        
    % noise
    g_weights_1 = conv(g_weights_1,n_weights,'same');        
    % update
    g_weights_1 = g_weights_1 .* g_likelihood_1;
    g_weights_1 = g_weights_1 ./ sum(g_weights_1);
    
    % noise
    g_weights_2 = conv(g_weights_2,n_weights,'same');        
    % update
    g_weights_2 = g_weights_2 .* g_likelihood_2;
    g_weights_2 = g_weights_2 ./ sum(g_weights_2);
    
    
    % noise
    gt_weights = conv(gt_weights,n_weights,'same');        
    % update
    gt_weights = gt_weights .* g_likelihood_1 .* g_likelihood_2;
    gt_weights = gt_weights ./ sum(gt_weights);

end

g_pf_weights = (g_pf_weights_1.*g_pf_weights_2);
g_pf_weights = g_pf_weights ./ sum(g_pf_weights);

g_fused_weights = (g_weights_1.*g_weights_2);
g_fused_weights = exp( log(g_fused_weights)-log(g_pf_weights) );    
g_fused_weights = g_fused_weights ./ sum(g_fused_weights);

g_gmd = sqrt(g_weights_1.*g_weights_2);
g_gmd = g_gmd ./ sum(g_gmd);

kl_cor = calculate_kl_divergence(g_fused_weights,gt_weights)
kl_gmd = calculate_kl_divergence(g_gmd,gt_weights)

md_gt = max_density_estimation(g_samples,gt_weights,([0.1 0.1]))
md_cor = max_density_estimation(g_samples,g_fused_weights,([0.1 0.1]))
md_gmd = max_density_estimation(g_samples,g_gmd,([0.1 0.1]))

% plot
if doPlot
       
%     figure(1); clf; hold on;
%        
%     plot3(g_samples(1,:),g_samples(2,:), g_weights_org, 'Color', [0.5 0 0.6], 'LineWidth', 1.2);
%     plot3(g_samples(1,:),g_samples(2,:), g_weights_1, 'Color', [0 0 1], 'LineWidth', 1.2);
%     plot3(g_samples(1,:),g_samples(2,:), g_weights_2, 'Color', [0 0 1], 'LineWidth', 1.2, 'LineStyle', '-.');
%     plot3(g_samples(1,:),g_samples(2,:), n_weights, 'Color', [0 1 0], 'LineWidth', 1.2);
%     plot3(g_samples(1,:),g_samples(2,:), gt_weights, 'Color', [1 0 1], 'LineWidth', 1.2);
%     plot3(g_samples(1,:),g_samples(2,:), g_fused_weights, 'Color', [0 0 0], 'LineWidth', 1.2);
%     plot3(g_samples(1,:),g_samples(2,:), g_naive, 'Color', [0 0 0], 'LineWidth', 1.2, 'LineStyle', '-.');
%     plot3(g_samples(1,:),g_samples(2,:), g_pf_weights, 'Color', [0 1 1], 'LineWidth', 1.2);
%     
%     legend("Original", "W1", "W2", "Noise", "GT", "Fused", "Naive", "Common");

    fig1 = figure(1); clf; 
    subplot(3,2,1); hold on; view(3);
    subtitle("Central Filter, MD: "+md_gt);
    plot_grid_data(fig1,g_samples,gt_weights);
    subplot(3,2,3); hold on; view(3);
    subtitle("Geometric Mean, MD: "+md_gmd);
    plot_grid_data(fig1,g_samples,g_gmd);
    subplot(3,2,5); hold on; view(3);
    subtitle("Information Filter, MD: "+md_cor);
    plot_grid_data(fig1,g_samples,g_fused_weights)
    
    subplot(3,2,2); hold on; view(3);
    subtitle("Common Information");
    plot_grid_data(fig1,g_samples,g_pf_weights);
    subplot(3,2,4); hold on; view(3);
    subtitle("Error GMD, KL: "+kl_gmd);
    plot_grid_data(fig1,g_samples,g_gmd-gt_weights);
    subplot(3,2,6); hold on; view(3);
    subtitle("Error IF, KL: "+kl_cor);
    plot_grid_data(fig1,g_samples,g_fused_weights-gt_weights);
    
%     figure(2); clf; hold on;
%     stem3(g_samples(1,:),g_samples(2,:), g_weights_1, 'Color', [0 0 1], 'LineWidth', 1.2);
%     stem3(g_samples(1,:),g_samples(2,:), g_weights_2, 'Color', [1 0 0], 'LineWidth', 1.2, 'LineStyle', '-.');
%     
%     figure(3); clf; hold on;
%     plot3(g_samples(1,:),g_samples(2,:), g_likelihood_1.*g_weights_1, 'Color', [0 0 1], 'LineWidth', 1.2);
%     plot3(g_samples(1,:),g_samples(2,:), g_likelihood_2.*g_weights_1, 'Color', [1 0 0], 'LineWidth', 1.2, 'LineStyle', '-.');
end

end



