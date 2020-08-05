%
% Agent-based model of constrained particle motion on a sinusoidal surface
% Authors: Tej Stead, Dhananjay Bhaskar
% Last Modified: Jul 24, 2020
%

% clear memory
close all; clear all;
% seed RNG
rng(4987)

% number of particles
N = 50;

% general params
deltaT = 0.1;
totT = 30;


% toggle interaction forces
FORCE_EUCLIDEAN_REPULSION_ON = false;
FORCE_ATTR_REPULSION_ON = false;
FORCE_RANDOM_POLARITY_ON = false;
FORCE_CUCKER_SMALE_POLARITY_ON = false;
FORCE_CURVATURE_ALIGNMENT_ON = true;

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
walk_amplitude = 0.5;
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

if(FORCE_CURVATURE_ALIGNMENT_ON)
    num_neighbors = 8; % used for curvature alignment
else
    num_neighbors = 1; % used as argument to knnsearch for closest mesh pt
end

% Supported modes:
% 'gauss-min' - align in direction of minimum Gaussian curvature
% 'gauss-max' - align in direction of maximum Gaussian curvature
% 'gauss-zero' - align in direction of lowest absolute Gaussian curvature
% 'mean-min' - align in direction of minimum mean curvature
% 'mean-max' - align in direction of maximum mean curvature
% 'mean-zero' - align in direction of lowest absolute mean curvature
alignment_mode = 'gauss-min';
alignment_magnitude = 0.02;

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
q = [2, 4];
Q = [0.0, 0.0];
F = zeros(N, 1);
dFdX = zeros(N, 3);
dFdq = zeros(N, 2);
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
while cnt < N
    U = rand();
    V = rand();
    Theta = 2*pi*U*a + 2*pi*a;
    Phi = 2*pi*V*a + 2*pi*a;
    cnt = cnt + 1;
    X(cnt, :) = [Theta, Phi, b*sin(Theta/a)*sin(Phi/a)];
end

% preload pairwise geodesic distances between mesh points (for static surfaces)
if(isfile("sinusoidal_mesh.mat"))
    load("sinusoidal_mesh.mat");
else
    mesh_theta_num = 90;
    mesh_phi_num = 45;
    theta_grid = linspace(0, 6*pi*a, mesh_theta_num);
    phi_grid = linspace(0, 6*pi*a, mesh_phi_num);
    [Phi_mesh_fine, Theta_mesh_fine] = meshgrid(phi_grid, theta_grid); 
    mesh_x = Theta_mesh_fine;
    mesh_y = Phi_mesh_fine;
    mesh_z = b .* sin(Theta_mesh_fine ./ a) .* sin(Phi_mesh_fine ./ a);
    mat = adj_mat_alternating_mountain(mesh_x, mesh_y, mesh_z);
    [dist_mat, next] = FloydWarshall(mat);
    save sinusoidal_mesh.mat mesh_theta_num mesh_phi_num mesh_x mesh_y ...
        Phi_mesh_fine Theta_mesh_fine mesh_z mat dist_mat next;
end
dist_range = [0 max(dist_mat(:))];

% preload coarse mesh for visualization
theta_num = 40;
phi_num = 20;
theta_grid = linspace(0, 6*pi*a, theta_num);
phi_grid = linspace(0, 6*pi*a, phi_num);
[Phi_mesh, Theta_mesh] = meshgrid(phi_grid, theta_grid); 
vis_x = Theta_mesh;
vis_y = Phi_mesh;
vis_z = b .* sin(Theta_mesh ./ a) .* sin(Phi_mesh ./ a);

% compute mean and gaussian curvature
G_curvature = gaussian_curvature_sinusoidal(Theta_mesh_fine, Phi_mesh_fine, q);
G_color_limits = [min(min(G_curvature)) max(max(G_curvature))];
M_curvature = mean_curvature_sinusoidal(Theta_mesh_fine, Phi_mesh_fine, q);
M_color_limits = [min(min(M_curvature)) max(max(M_curvature))];

% visualize IC
% visualize_surface(X, 0, vis_x, vis_y, vis_z, [-5 30], [-5 30], [-10 10]);
% visualize_geodesic_path(X, 0, pt_1_idx, pt_2_idx, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, mesh_phi_num, next, [-5 30], [-5 30], [-10 10]);
% visualize_geodesic_heatmap(X, 0, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, pt_1_idx, [-5 30], [-5 30], [-10 10], dist_range, dist_mat);
visualize_curvature_heatmap(X, 0, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, [-5 30], [-5 30], [-10 10], G_color_limits, G_curvature, true);
% visualize_trajectories(X, 0, prev_paths, path_colors, vis_x, vis_y, vis_z, [-5 30], [-5 30], [-10 10]);
t = 0;
itr = 0;

while t < totT
    [indices, dists] = all_mesh_neighbors(X, mesh_x, mesh_y, mesh_z, num_neighbors);
    % compute updated state vectors
    for i = 1 : N
        F(i) = (X(i,3)) - b * sin(X(i,2)/a) * sin(X(i, 1)/a);

        dFdX_i_x = -1 * (b/a) * sin(X(i,2)/a) * cos(X(i,1)/a);
        dFdX_i_y = -1 * (b/a) * sin(X(i,1)/a) * cos(X(i,2)/a);
        dFdX_i_z = 1;
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
                particle_indices = find_neighbors(X, indices(:, 1), dists(:, 1), dist_mat, i, CS_threshold);
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
            prev_CS_buffer(i, :) = deltaP;
            P(i, :) = P(i, :) + deltaP;
        end
        if(FORCE_CURVATURE_ALIGNMENT_ON)
            neighbor_indices = indices(i, :);
            switch alignment_mode
                case 'gauss-min'
                    [~, neighbor_idx] = min(G_curvature(neighbor_indices));
                case 'gauss-max'
                    [~, neighbor_idx] = max(G_curvature(neighbor_indices));
                case 'gauss-zero'
                    [~, neighbor_idx] = min(abs(G_curvature(neighbor_indices)));
                case 'mean-min'
                    [~, neighbor_idx] = min(M_curvature(neighbor_indices));
                case 'mean-max'
                    [~, neighbor_idx] = max(M_curvature(neighbor_indices));
                case 'mean-zero'
                    [~, neighbor_idx] = min(abs(M_curvature(neighbor_indices)));
                otherwise
                    error("Invalid curvature alignment mode.");
            end
            neighbor_point_idx = neighbor_indices(neighbor_idx);
            direction = [mesh_x(neighbor_point_idx) mesh_y(neighbor_point_idx) mesh_z(neighbor_point_idx)] - X(i, :);
            direction = direction/norm(direction);
            deltaP = direction * alignment_magnitude;
            P(i, :) = P(i, :) + deltaP;
        end
        
        dFdq_i_a = (b/a^2) * ((X(i,2)*sin(X(i,1)/a)*cos(X(i,2)/a)) + X(i,1)*cos(X(i,1)/a)*sin(X(i,2)/a));
        dFdq_i_b = -1 * sin(X(i,1)/a) * sin(X(i,2)/a);
        dFdq(i,:) = [dFdq_i_a, dFdq_i_b];

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
    
%     visualize_surface(X, itr, vis_x, vis_y, vis_z, [-5 30], [-5 30], [-10 10]);
    % visualize_geodesic_path(X, itr, pt_1_idx, pt_2_idx, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, mesh_phi_num, next, [-5 30], [-5 30], [-10 10]);
    % visualize_geodesic_heatmap(X, itr, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, pt_1_idx, [-5 30], [-5 30], [-10 10], dist_range, dist_mat);
    visualize_curvature_heatmap(X, itr, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, [-5 30], [-5 30], [-10 10], G_color_limits, G_curvature, true);
%     visualize_trajectories(X, itr, prev_paths, path_colors, vis_x, vis_y, vis_z, [-5 30], [-5 30], [-10 10]);
end

function [adj_mat] = adj_mat_alternating_mountain(x, y, z)

    sz = size(x);
    height = sz(1);
    width = sz(2);
    adj_mat = inf*ones(height*width, height*width);
    
    for i = 1:height
        for j = 1:width
            dx = [-1 -1 -1 0 0 0 1 1 1];
            dy = [-1 0 1 -1 0 1 -1 0 1];
            for k = 1:numel(dx)
                new_i = i + dy(k);
                new_j = j + dx(k);
                if(new_i < 1 || new_i > height)
                    new_i = i;
                end
                if(new_j < 1 || new_j > width)
                    new_j = j;
                end
                distance = pdist([x(i,j) y(i,j) z(i,j) ; x(new_i, new_j) y(new_i, new_j) z(new_i, new_j)]);
                adj_mat((i-1)*width + j,(new_i - 1)*width + new_j) = distance;
            end
        end
    end
    
end

function [curvature] = gaussian_curvature_sinusoidal(Theta_mesh_fine, Phi_mesh_fine, q)

    % Reference: https://en.wikipedia.org/wiki/Parametric_surface#Curvature
    a = q(1);
    b = q(2);
    u = Theta_mesh_fine;
    v = Phi_mesh_fine;
    
    K = sqrt(a^2 + b^2 .* ((cos(u/a).^2).*(sin(v/a).^2) + (sin(u/a).^2) .* (cos(v/a).^2)));
    L = (-b/a) .* sin(u/a) .* sin(v/a) ./ K;
    M = (b/a) .* cos(u/a) .* cos(v/a) ./ K;
    % Note: N == L
    E = (b/a)^2 .*(cos(u/a).^2) .* (sin(v/a).^2) + 1;
    F = (b/a)^2 .* cos(u/a) .* sin(u/a) .* cos(v/a) .* sin(v/a);
    G = (b/a)^2 .*(sin(u/a).^2) .* (cos(v/a).^2) + 1;
    
    A = (E .* G - L.^2);
    B = (2 .* M .* F - L .* G - E .* L);
    C = (L.^2 - M.^2);
    
    assert(min(min(B.^2 - 4 .* A .* C)) > 0);
    
    curv_1 = ((-1 .* B) + sqrt(B.^2 - 4 .* A .* C))./(2 .* A);
    curv_2 = ((-1 .* B) - sqrt(B.^2 - 4 .* A .* C))./(2 .* A);
    
    curvature = curv_1 .* curv_2;
    
end

function [curvature] = mean_curvature_sinusoidal(Theta_mesh_fine, Phi_mesh_fine, q)

    a = q(1);
    b = q(2);
    u = Theta_mesh_fine;
    v = Phi_mesh_fine;
    
    K = sqrt(a^2 + b^2 .* ((cos(u/a).^2).*(sin(v/a).^2) + (sin(u/a).^2) .* (cos(v/a).^2)));
    L = (-b/a) .* sin(u/a) .* sin(v/a) ./ K;
    M = (b/a) .* cos(u/a) .* cos(v/a) ./ K;
    % Note: N == L
    E = (b/a)^2 .*(cos(u/a).^2) .* (sin(v/a).^2) + 1;
    F = (b/a)^2 .* cos(u/a) .* sin(u/a) .* cos(v/a) .* sin(v/a);
    G = (b/a)^2 .*(sin(u/a).^2) .* (cos(v/a).^2) + 1;
    
    A = (E .* G - L.^2);
    B = (2 .* M .* F - L .* G - E .* L);
    C = (L.^2 - M.^2);
    
    assert(min(min(B.^2 - 4 .* A .* C)) > 0);
    
    curv_1 = ((-1 .* B) + sqrt(B.^2 - 4 .*A .* C))./(2 .* A);
    curv_2 = ((-1 .* B) - sqrt(B.^2 - 4 .* A .* C))./(2 .* A);
    
    curvature = (curv_1 + curv_2) ./ 2;
    
end