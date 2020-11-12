%*********************************************
% do monte-carlo runs of tracker
%*********************************************

addpath(genpath('../'));

doPlot = true;

Nc = 101;
limits = [-5 5];

gauss_mean = [-1.5, 1.5];
gauss_sig = diag(0.25);

noise_sig = diag(0.25);

gauss_count = size(gauss_mean, 2);
Ds = size(gauss_mean, 1);

% samples
g_samples = linspace(limits(1,1), limits(1,2), Nc);
g_weights = zeros(1,Nc);
for i=1:gauss_count
   g_weights = g_weights + mvnpdf(g_samples', gauss_mean(:,i)', gauss_sig')';
end
g_weights = g_weights ./ sum(g_weights);

% noise
n_weights = mvnpdf(g_samples', zeros(Ds,1)', noise_sig');
n_weights = n_weights ./ sum(n_weights);


% conv
g_noised_weights = conv(g_weights,n_weights,'same');
g_noised_weights = g_noised_weights ./ sum(g_noised_weights);

% conv gaussian
g_conv_weights = zeros(1,Nc);
for i=1:gauss_count
   g_conv_weights = g_conv_weights + conv(mvnpdf(g_samples', gauss_mean(:,i)', gauss_sig')',n_weights,'same');
end
g_conv_weights = g_conv_weights ./ sum(g_conv_weights);


% plot
if doPlot
   
    figure(1); clf; hold on;
    plot(g_samples, g_weights, 'Color', [0 0 1], 'LineWidth', 1.2);
    plot(g_samples, n_weights, 'Color', [0 1 0], 'LineWidth', 1.2);
    plot(g_samples, g_noised_weights, 'Color', [1 0 1], 'LineWidth', 1.2);
    plot(g_samples, g_conv_weights, 'Color', [0 0 0], 'LineWidth', 1.2);
    legend("Original", "Noise", "Noised PDF");
    
end