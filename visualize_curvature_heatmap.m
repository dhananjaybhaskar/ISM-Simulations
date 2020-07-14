function visualize_curvature_heatmap(X, itr, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, ...
                                    x_limits, y_limits, z_limits, c_limits, curvature)

%     fig = figure('visible', 'off');
    figure(1);
    mesh(vis_x, vis_y, vis_z, 'EdgeColor', 'k', 'FaceColor', 'none', 'FaceAlpha', 0, 'linestyle', '-');
    hold on;
    surf(mesh_x, mesh_y, mesh_z, curvature, 'EdgeColor', 'none', 'FaceAlpha', 0.75, 'FaceLighting', 'gouraud', 'LineStyle', 'none');
    shading interp
    daspect([1 1 1])
    xlim(x_limits)
    ylim(y_limits)
    zlim(z_limits)
    scatter3(X(:,1), X(:,2), X(:,3), 30, 'filled');
    colorbar
    caxis(c_limits)
    zoom(1.2)
    hold off;
    pause;
%     fname = strcat('curvature', sprintf('%03d',itr), '.png');
%     saveas(fig, fname, 'png');
    
end