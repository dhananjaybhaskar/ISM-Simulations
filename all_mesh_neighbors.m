function [indices, dists] = all_mesh_neighbors(X, mesh_x, mesh_y, mesh_z, num_neighbors)
    new_len = numel(mesh_x);
    new_mat = [reshape(mesh_x, new_len, 1) reshape(mesh_y, new_len, 1) reshape(mesh_z, new_len, 1)];
    [indices, dists] = knnsearch(new_mat, X, 'K', num_neighbors);
end