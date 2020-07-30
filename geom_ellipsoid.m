%
% Agent-based model of constrained particle motion on an ellipsoid
% Authors: Tej Stead, Dhananjay Bhaskar
% Last Modified: Jul 14, 2020
%

% clear memory
close all; clear all;
% seed RNG
rng(4987)

% number of particles
N = 50;

% general params
deltaT = 0.1;
totT = 60;


% toggle interaction forces
FORCE_EUCLIDEAN_REPULSION_ON = false;
FORCE_ATTR_REPULSION_ON = true;
FORCE_RANDOM_POLARITY_ON = true;
FORCE_CUCKER_SMALE_POLARITY_ON = true;

% init positions
X = zeros(N, 3);

% Euclidean repulsion force params 
Alpha = 2;
Sigma = 0.5;
phi = 1;

% Attraction-repulsion force params
C_a = 1;
C_r = 1;
l_a = 6;
l_r = 0.5;
% random polarization params
walk_amplitude = 0.3;
walk_stdev = pi/4;
walk_direction = rand(N, 1) * 2 * pi;
num_repolarization_steps = 10;
init_repolarization_offset = floor(rand(N, 1) * num_repolarization_steps);
num_trailing_positions = 40;

% Cucker-Smale polarization/flocking params
CS_K = 2;
CS_Sigma = 1;
CS_Gamma = 1.5;
CS_threshold = 5;
use_nearest_neighbors = true;

% preallocate state variables
P = zeros(N, 3);
prev_EF_buffer = zeros(N, 3);
prev_AR_buffer = zeros(N, 3);
prev_RP_buffer = zeros(N, 3);
prev_CS_buffer = zeros(N, 3);
prev_EF = zeros(N, 3);
prev_AR = zeros(N, 3);
prev_RP = zeros(N, 3);
prev_CS = zeros(N, 3);
q = [20, 15, 20];
Q = [0.0, 0.0, 0.0];
F = zeros(N, 1);
dFdX = zeros(N, 3);
dFdq = zeros(N, 3);
dXdt = zeros(N, 3);

% preallocate particle trajectories
prev_paths = nan * ones(N, num_trailing_positions, 3);
path_colors = hsv(N);

% pick trajectory colors
for i = 1:N
    j = floor(rand() * N) + 1;
    temp = path_colors(j, :);
    path_colors(j, :) = path_colors(i, :);
    path_colors(i, :) = temp;
end

% pick random particles
pt_1_idx = floor(rand()*N) + 1;
pt_2_idx = floor(rand()*N) + 1;
pt_3_idx = floor(rand()*N) + 1;
pt_4_idx = floor(rand()*N) + 1;



% uniform distribution of (Theta, Phi) in [0, 2pi] for initial position
cnt = 0;
a = q(1);
b = q(2);
c = q(3);
while cnt < N
    U = rand();
    V = rand();
    Theta = 2*pi*U;
    Phi = pi*V;
    cnt = cnt + 1;
    X(cnt, :) = [a*cos(Theta)*sin(Phi), b*sin(Theta)*sin(Phi), c*cos(Phi)];
end

% preload pairwise geodesic distances between mesh points (for static surfaces)
if(isfile("ellipsoid_mesh.mat"))
    load("ellipsoid_mesh.mat");
else
    mesh_theta_num = 80;
    mesh_phi_num = 40;
    theta_grid = linspace(0, 2*pi, mesh_theta_num);
    phi_grid = linspace(0, pi, mesh_phi_num);
    [Phi_mesh_fine, Theta_mesh_fine] = meshgrid(phi_grid, theta_grid); 
    mesh_x = a.*cos(Theta_mesh_fine).*sin(Phi_mesh_fine);
    mesh_y = b.*sin(Theta_mesh_fine).*sin(Phi_mesh_fine);
    mesh_z = c.*cos(Phi_mesh_fine);
    mat = adj_mat_ellipsoid(mesh_x, mesh_y, mesh_z);
    [dist_mat, next] = FloydWarshall(mat);
    save ellipsoid_mesh.mat mesh_theta_num mesh_phi_num mesh_x mesh_y ...
        Phi_mesh_fine Theta_mesh_fine mesh_z mat dist_mat next;
end
dist_range = [0 max(dist_mat(:))];

% preload coarse mesh for visualization
theta_num = 36;
phi_num = 18;
theta_grid = linspace(0, 2 * pi, theta_num);
phi_grid = linspace(0, pi, phi_num);
[Phi_mesh, Theta_mesh] = meshgrid(phi_grid, theta_grid); 
vis_x = a .* cos(Theta_mesh) .* sin(Phi_mesh);
vis_y = b .* sin(Theta_mesh) .* sin(Phi_mesh);
vis_z = c .* cos(Phi_mesh);
    
% compute mean and gaussian curvature
G_curvature = gaussian_curvature_ellipsoid(Theta_mesh_fine, Phi_mesh_fine, q);
G_color_limits = [0 max(max(G_curvature))];
M_curvature = mean_curvature_ellipsoid(Theta_mesh_fine, Phi_mesh_fine, q);
M_color_limits = [min(min(M_curvature)) max(max(M_curvature))];

% visualize IC
% visualize_surface(X, 0, vis_x, vis_y, vis_z, [-30 30], [-30 30], [-30 30]);
% visualize_geodesic_path(X, 0, [pt_1_idx pt_3_idx], [pt_2_idx pt_4_idx], vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, mesh_phi_num, next, [-30 30], [-30 30], [-30 30]);
% visualize_geodesic_heatmap(X, 0, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, pt_1_idx, [-30 30], [-30 30], [-30 30], dist_range, dist_mat);
% visualize_curvature_heatmap(X, 0, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, [-30 30], [-30 30], [-30 30], G_color_limits, G_curvature, false);
visualize_trajectories(X, 0, prev_paths, path_colors, vis_x, vis_y, vis_z, [-30 30], [-30 30], [-30 30]);
t = 0;
itr = 0;

while t < totT
    [indexes, dists] = all_mesh_neighbors(X, mesh_x, mesh_y, mesh_z);
    % compute updated state vectors
    for i = 1 : N

        F(i) = (X(i,1)^2/a^2) + (X(i, 2)^2/b^2) + (X(i,3)^2/c^2) - 1;

        dFdX_i_x = 2*X(i,1)/(a^2);
        dFdX_i_y = 2*X(i,2)/(b^2);
        dFdX_i_z = 2*X(i,3)/(c^2);
        dFdX(i,:) = [dFdX_i_x, dFdX_i_y, dFdX_i_z];
        
        if (FORCE_EUCLIDEAN_REPULSION_ON)
            for j = setdiff(1:N, i)
                Fij = Alpha*exp(-1.0*(norm((X(i,:)-X(j,:)))^2)/(2*Sigma^2));
                deltaP = (X(i,:) - X(j,:))*Fij;
                prev_EF_buffer(i, :) = deltaP;
                P(i, :) = P(i, :) + deltaP;
            end
        end
        
        if (FORCE_ATTR_REPULSION_ON) 
            % since this has a long-range attraction force, we don't only
            % want to use nearest neighbors
            dPdt = 0;
            if(N > 1)
                for j = setdiff(1:N, i) % skips element i 
                    diff = X(i, :) - X(j, :);
                    dist = norm(diff);
                    grad_x = ((C_a * diff(1) * exp(-1 * (dist/l_a)))/(l_a * dist)) ...
                        - ((C_r * diff(1) * exp(-1 * (dist/l_r)))/(l_r * dist));
                    grad_y = ((C_a * diff(2) * exp(-1 * (dist/l_a)))/(l_a * dist)) ...
                        - ((C_r * diff(2) * exp(-1 * (dist/l_r)))/(l_r * dist));
                    grad_z = ((C_a * diff(3) * exp(-1 * (dist/l_a)))/(l_a * dist)) ...
                        - ((C_r * diff(3) * exp(-1 * (dist/l_r)))/(l_r * dist));
                    dPdt = dPdt - (1/(N - 1)) * [grad_x grad_y grad_z];
                end
            end
            deltaP = deltaT * dPdt;
            prev_AR_buffer(i, :) = deltaP;
            P(i, :) = P(i, :) + deltaP;
        end
        if (FORCE_RANDOM_POLARITY_ON)
            nullspace = [dFdX(i,:); zeros(2,3)];
            assert(numel(nullspace) == 9, "Nullspace computation error.");
            nullspace = null(nullspace);
            nullspace = nullspace';
            if(mod(itr, num_repolarization_steps) == init_repolarization_offset(i))
                temp = rand();
                walk_direction(i) = walk_direction(i) + norminv(temp, 0, walk_stdev);
            end
            deltaP =  cos(walk_direction(i)) * nullspace(1, :) * walk_amplitude;
            deltaP = deltaP + sin(walk_direction(i)) * nullspace(2, :) * walk_amplitude;
            prev_RP_buffer(i, :) = deltaP;
            P(i, :) = P(i, :) + deltaP;
        end
        
        if (FORCE_CUCKER_SMALE_POLARITY_ON)
            dPdt = [0 0 0];
            % compute the Cucker-Smale polarity
            if(use_nearest_neighbors)
                particle_indices = find_neighbors(X, indexes, dists, dist_mat, i, CS_threshold);
            else
                particle_indices = setdiff(1:N, i);
            end
            sz = numel(particle_indices);
            prev_p_i = prev_EF(i, :) + prev_RP(i, :) + prev_CS(i, :);
            for j = 1:sz
                idx = particle_indices(j);
                dist = norm(X(i, :) - X(idx, :));
                CS_H = CS_K/((CS_Sigma^2) + (dist^2))^CS_Gamma;
                prev_p_j = prev_EF(idx, :) + prev_AR(idx, :) + prev_RP(idx, :) + prev_CS(idx, :);
                dPdt = dPdt + (1/sz) .* CS_H .* (prev_p_j - prev_p_i);
            end
            deltaP = deltaT * dPdt;
            prev_CS(i, :) = deltaP;
            P(i, :) = P(i, :) + deltaP;
        end

        dFdq_i_a = -2*X(i,1)^2/(a^3);
        dFdq_i_b = -2*X(i,2)^2/(b^3);
        dFdq_i_c = -2*X(i,3)^2/(c^3);
        dFdq(i,:) = [dFdq_i_a, dFdq_i_b, dFdq_i_c];

        correction = (dot(dFdX(i,:), P(i,:)) + dot(dFdq(i,:), Q) + phi*F(i))/(norm(dFdX(i,:))^2);
        dXdt(i,:) = P(i,:) - correction*dFdX(i,:);
        if(itr < num_trailing_positions)
            prev_paths(i, itr + 1, :) = X(i, :);
        else
            prev_paths(i, 1:(num_trailing_positions - 1), :) = prev_paths(i, 2:num_trailing_positions, :);
            prev_paths(i, num_trailing_positions, :) = X(i, :);
        end
    end
    
    % update position
    for i = 1 : N
        X(i,:) = X(i,:) + deltaT*dXdt(i,:);
    end
    prev_EF = prev_EF_buffer;
    prev_AR = prev_AR_buffer;
    prev_RP = prev_RP_buffer;
    prev_CS = prev_CS_buffer;
    t = t + deltaT;
    itr = itr + 1;
    
    % visualize_surface(X, itr, vis_x, vis_y, vis_z, [-30 30], [-30 30], [-30 30]);
    % visualize_geodesic_path(X, itr, [pt_1_idx pt_3_idx], [pt_2_idx pt_4_idx], vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, mesh_phi_num, next, [-30 30], [-30 30], [-30 30]);
%     visualize_geodesic_heatmap(X, itr, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, pt_1_idx, [-30 30], [-30 30], [-30 30], dist_range, dist_mat);
    % visualize_curvature_heatmap(X, itr, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, [-30 30], [-30 30], [-30 30], G_color_limits, G_curvature, true);
    visualize_trajectories(X, itr, prev_paths, path_colors, vis_x, vis_y, vis_z, [-30 30], [-30 30], [-30 30]);
end

function [adj_mat] = adj_mat_ellipsoid(x, y, z)

    sz = size(x);
    height = sz(1);
    width = sz(2);
    adj_mat = inf*ones(height*width, height*width);
    
    for i = 1:height
        for j = 1:width
            dx = [-1 -1 -1 0 0 0 1 1 1];
            dy = [-1 0 1 -1 0 1 -1 0 1];
            for k = 1:numel(dx)
                new_i = mod(i+dy(k) - 1, height) + 1; 
                new_j = j + dx(k);
                if(new_j < 1 || new_j > width)
                    new_j = j;
                end
                distance = pdist([x(i,j) y(i,j) z(i,j) ; x(new_i, new_j) y(new_i, new_j) z(new_i, new_j)]);
                adj_mat((i-1)*width + j,(new_i - 1)*width + new_j) = distance;
            end
        end
    end
    
end

function [curvature] = gaussian_curvature_ellipsoid(Theta_mesh_fine, Phi_mesh_fine, q)

    % Reference: https://mathworld.wolfram.com/Ellipsoid.html 
    % Using polar coordinates for consistency with mean curvature, implementation is
    % tidier with Cartesian coordinates
    a = q(1);
    b = q(2);
    c = q(3);
    
    num = (a^2) * (b^2) * (c^2);
    
    denom = (((a^2) * (b^2) * (cos(Phi_mesh_fine).^2)) ... 
        + ((c^2) * (((b^2) * (cos(Theta_mesh_fine).^2)) ... 
        + ((a^2) * (sin(Theta_mesh_fine).^2))) .* (sin(Phi_mesh_fine).^2))).^2;
    
    curvature = num./denom;
    
end

function [curvature] = mean_curvature_ellipsoid(Theta_mesh_fine, Phi_mesh_fine, q)

    a = q(1);
    b = q(2);
    c = q(3);
    
    num = (a * b * c) * ((3 * ((a^2) + (b^2))) + (2 * (c^2)) ...
        + (((a^2) + (b^2) - 2 * (c^2)) * cos(2 * Phi_mesh_fine)) ...
        - (2 * ((a^2) - (b^2)) * cos(2 * Theta_mesh_fine) .* (sin(Phi_mesh_fine).^2)));
    
    denom = 8*(((a^2) * (b^2) * (cos(Phi_mesh_fine).^2)) ... 
        + ((c^2) * (((b^2) * (cos(Theta_mesh_fine).^2)) ... 
        + ((a^2) * (sin(Theta_mesh_fine).^2))) .* (sin(Phi_mesh_fine).^2))).^(3/2);
    
    curvature = num./denom;
    
end