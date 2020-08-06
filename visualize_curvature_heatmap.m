function visualize_curvature_heatmap(X, itr, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, ...
                                    x_limits, y_limits, z_limits, color_limits, curvature, display_particles)

    fig = figure('visible', 'off');
    surf(mesh_x, mesh_y, mesh_z, curvature, 'EdgeColor', 'interp', 'FaceColor', 'interp', ...
        'FaceAlpha', 0.75, 'FaceLighting', 'gouraud', 'LineStyle', 'none');
    shading interp
    hold on;
    mesh(vis_x, vis_y, vis_z, 'EdgeColor', 'none', 'FaceColor', 'none', 'FaceAlpha', 0, 'linestyle', '-');
    daspect([1 1 1])
    xlim(x_limits)
    ylim(y_limits)
    zlim(z_limits)
    if display_particles
        scatter3(X(:,1), X(:,2), X(:,3), 30, 'ro', 'filled');
    end
    colorbar
    caxis(color_limits)
    zoom(1.2)
    hold off;
    set(gca, 'visible', 'off')
    fname = strcat('sim_', sprintf('%03d', itr), '.png');
    saveas(fig, fname, 'png');
    close
    
end