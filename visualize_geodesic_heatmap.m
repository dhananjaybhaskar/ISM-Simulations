function visualize_geodesic_heatmap(X, itr, mesh_x, mesh_y, mesh_z, pt_index, x_limits, y_limits, z_limits, dist_mat, color_axis)

    source_pt_coords = X(pt_index, :);
    [min_dist, mesh_idx] = closest_mesh_point(source_pt_coords, mesh_x,mesh_y,mesh_z);
    sz = size(mesh_x);
    heatmap = zeros(sz);
    for i = 1:sz(1)
        for j = 1:sz(2)
            dest_mesh_idx = (i - 1)*sz(2) + j;
            heatmap(i,j) = min_dist + dist_mat(mesh_idx, dest_mesh_idx);
        end
    end
    fig = figure('visible', 'off');
%     figure(1);
    surf(mesh_x, mesh_y, mesh_z, heatmap);
    caxis(color_axis);
    colorbar
    hold on;
%     mesh(vis_x, vis_y, vis_z, 'edgecolor', 'k');
    alpha 0.75;
    daspect([1 1 1])
    xlim(x_limits)
    ylim(y_limits)
    zlim(z_limits)
    scatter3(X(:,1), X(:,2), X(:,3));
    plot3(source_pt_coords(1), source_pt_coords(2), source_pt_coords(3), 'or');
    hold off;
    fname = strcat('sim_', sprintf('%03d',itr), '.png');
    saveas(fig, fname, 'png');
end