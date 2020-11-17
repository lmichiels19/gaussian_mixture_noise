function plot_grid_data(fig,g_samples,g_weights_org)
%PLOT_GRID_DATA Plots estimation grid data

    Np = size(g_samples,2);

    x_g = reshape(g_samples(1,:),[sqrt(Np),sqrt(Np)]);
    y_g = reshape(g_samples(2,:),[sqrt(Np),sqrt(Np)]);
    z_g = reshape(g_weights_org,[sqrt(Np),sqrt(Np)]);
    
    figure(fig);
    surf(x_g,y_g,z_g,'LineStyle','none');
end

