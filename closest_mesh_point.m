function [min_dist, mesh_idx] = closest_mesh_point(point_coords, x, y, z)

    sz = size(x);
    min_dist = inf;
    
    for i = 1:sz(1)
        for j = 1:sz(2)
            d = pdist([point_coords ; x(i,j) y(i,j) z(i,j)]);
            if(d < min_dist)
                min_dist = d;
                min_i = i;
                min_j = j;
            end
        end
    end
    
    mesh_idx = (min_i - 1)*sz(2) + min_j;
    
end