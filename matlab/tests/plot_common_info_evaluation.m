%*********************************************
% do monte-carlo runs of tracker
%*********************************************

addpath(genpath('../'));

doPlot = true;

Nc = 501;
limits = [-10 10];
max_run_count = 10;

gauss_mean = [-2.5,1.5];
gauss_sig = diag([0.25]);

noise_sig = diag([0.5]);

%meas_mean_0 = 2.*randn(1,2);
%meas_mean_0 = [-1 2.5];
meas_sig = diag([0.2]);

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
g_pf_weights_prior = g_weights_1;
g_pf_weights_noise = g_weights_2;%ones(size(g_weights_1))./Np;

% noise
n_weights = mvnpdf(g_samples', zeros(Ds,1)', noise_sig')';
n_weights = n_weights ./ sum(n_weights);

gt_weights = g_weights_1;

meas_mean = meas_mean_0;

for runs=1:run_count

    gamma = 1;
    % noise
    g_pf_weights_prior = conv(g_pf_weights_prior,n_weights,'same');

    % noise
    if( runs > 1 )
        old_likelihood = g_likelihood;% .* g_weights_1.^0.1;
        g_pf_weights_noise = g_pf_weights_noise.*conv(old_likelihood,n_weights.^gamma,'same');
        
        figure(12); clf; hold on;
        plot(g_samples,old_likelihood, 'Color', [0 0.5 0], 'LineWidth', 1.2);
        pause(0);
    end
    g_pf_weights_noise = conv(g_pf_weights_noise,n_weights,'same');

    
    %likelihood
    g_likelihood = zeros(size(g_weights_1));
    li_count = size(meas_mean,2);
    mean_1 = mvnrnd(meas_mean', meas_sig);
    for i=1:li_count %(1/li_count).*
       g_likelihood = g_likelihood + mvnpdf(g_samples', mean_1(i,:), meas_sig')'.*sqrt(2*pi*det(meas_sig));
       meas_mean(:,i) = meas_mean(:,i) + 0.1;
    end
        
    % noise
    g_weights_1 = conv(g_weights_1,n_weights,'same');        
    % update
    g_old = g_weights_1;
    g_weights_1 = g_weights_1 .* g_likelihood;
    g_weights_1 = g_weights_1 ./ sum(g_weights_1);
    
    % noise
    g_weights_2 = conv(g_weights_2,n_weights,'same');        
    % update
    g_weights_2 = g_weights_2 .* g_likelihood;
    g_weights_2 = g_weights_2 ./ sum(g_weights_2);
    
    
    % noise
    gt_weights = conv(gt_weights,n_weights,'same');        
    % update
    gt_weights = gt_weights .* g_likelihood .* g_likelihood;
    gt_weights = gt_weights ./ sum(gt_weights);

    plot_steps = 3;
    if rem(runs,plot_steps) == 1 || plot_steps == 1
        gt_weights = gt_weights ./ sum(gt_weights);
        g_pf_weights_prior = g_pf_weights_prior ./ sum(g_pf_weights_prior);
        g_pf_weights_noise = g_pf_weights_noise ./ sum(g_pf_weights_noise);
        
        sb_count = ceil(run_count/plot_steps);
        idx =(ceil(runs/plot_steps));
        
        figure(11);
        subplot(3,sb_count,idx); cla; hold on;
        axis([limits(1) limits(2) 0 0.05]);
        subtitle("Filter State (Step: "+runs+")");
        plot(g_samples, gt_weights, 'Color', [0 0 1], 'LineWidth', 1.2);
        subplot(3,sb_count,idx+sb_count); cla; hold on;
        axis([limits(1) limits(2) 0 0.02]);
        subtitle("Prior Information (Step: "+runs+")");
        plot(g_samples, (g_pf_weights_prior), 'Color', [0 0.5 0], 'LineWidth', 1.2);
        subplot(3,sb_count,idx+2*sb_count); cla; hold on;
        axis([limits(1) limits(2) 0 0.02]);
        subtitle("Common Noise (Step: "+runs+")");
        plot(g_samples, (g_pf_weights_noise), 'Color', [0 0.5 0], 'LineWidth', 1.2);
        
    end
    
end

g_pf_weights = (g_pf_weights_noise);
g_pf_weights = g_pf_weights ./ sum(g_pf_weights);

g_fused_weights = (g_weights_1.*g_weights_2);
g_fused_weights = exp( log(g_fused_weights)-log(g_pf_weights) );    
g_fused_weights = g_fused_weights ./ sum(g_fused_weights);

g_naive = sqrt(g_weights_1.*g_weights_2);
%g_naive = g_naive./g_pf_weights;    
g_naive = g_naive ./ sum(g_naive);

% plot
if doPlot
       
    figure(1); clf; hold on;
       
    plot(g_samples, g_weights_org, 'Color', [0.5 0 0.6], 'LineWidth', 1.2);
    plot(g_samples, g_weights_1, 'Color', [0 0 1], 'LineWidth', 1.2);
    plot(g_samples, g_weights_2, 'Color', [0 0 1], 'LineWidth', 1.2, 'LineStyle', '-.');
    plot(g_samples, n_weights, 'Color', [0 1 0], 'LineWidth', 1.2);
    plot(g_samples, gt_weights, 'Color', [1 0 1], 'LineWidth', 1.2);
    plot(g_samples, g_fused_weights, 'Color', [0 0 0], 'LineWidth', 1.2);
    plot(g_samples, g_naive, 'Color', [0 0 0], 'LineWidth', 1.2, 'LineStyle', '-.');
    plot(g_samples, g_pf_weights, 'Color', [0 1 1], 'LineWidth', 1.2);
    
    legend("Original", "W1", "W2", "Noise", "GT", "Fused", "Naive", "Common");

%     figure(2); clf; hold on;
%     stem(g_samples, g_weights_1, 'Color', [0 0 1], 'LineWidth', 1.2);
%     stem(g_samples, g_weights_2, 'Color', [1 0 0], 'LineWidth', 1.2, 'LineStyle', '-.');
%     
%     figure(3); clf; hold on;
%     plot(g_samples, g_likelihood_1.*g_weights_1, 'Color', [0 0 1], 'LineWidth', 1.2);
%     plot(g_samples, g_likelihood_2.*g_weights_1, 'Color', [1 0 0], 'LineWidth', 1.2, 'LineStyle', '-.');
end

end



