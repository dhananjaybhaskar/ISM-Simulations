function [adj_mat] = build_hexagonal_mesh(stabilized_floaters)
    len = length(stabilized_floaters);
    adj_mat = inf*ones(len);
    for i = 1:len
        dists = zeros(len, 1);
        for j = 1:len
            dists(j) = norm(stabilized_floaters(i, :) - stabilized_floaters(j, :));
        end
        [dists, idxes] = mink(dists, 7);
        for j = 1:7
            adj_mat(i, idxes(j)) = dists(j);
%             adj_mat(idxes(j), i) = dists(j);
        end
    end
%     G = graph(adj_mat);
%     plot(G);
%     figure(1); clf; gplot3(adj_mat, stabilized_floaters);
    figure(1); scatter3(stabilized_floaters(:, 1), stabilized_floaters(:, 2), stabilized_floaters(:, 3));
    assert(issymmetric(adj_mat));
end