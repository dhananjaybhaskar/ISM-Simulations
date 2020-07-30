function [indexes, dists] = all_mesh_neighbors(X, mesh_x, mesh_y, mesh_z)
    new_len = numel(mesh_x);
    new_mat = [reshape(mesh_x, new_len, 1) reshape(mesh_y, new_len, 1) reshape(mesh_z, new_len, 1)];
    [indexes, dists] = knnsearch(new_mat, X);
end