%
% Agent-based model of particle movement on implicit surface (torus)
% Authors: Dhananjay Bhaskar, Tej Stead
% Last Modified: Aug 5, 2020
% Reference: Using Particles to Sample and Control Implicit Surfaces
% Andrew P. Witkin and Paul S. Heckbert, Proc. SIGGRAPH '94
%

close all; clear all;

% seed RNG
rng(42)

% number of particles
N = 120;

% simulation params
deltaT = 0.1;
totT = 80;
phi = 1;

% toggle interaction forces
FORCE_EUCLIDEAN_REPULSION_ON = false;
FORCE_ATTR_REPULSION_ON = true;
FORCE_RANDOM_POLARITY_ON = true;
FORCE_CUCKER_SMALE_POLARITY_ON = false;
FORCE_CURVATURE_ALIGNMENT_ON = false;

% init positions
X = zeros(N, 3);

% euclidean repulsion params
Alpha = 2;
Sigma = 0.5;

% attraction-repulsion params
C_a = 200;
C_r = 200;
l_a = 1;
l_r = 0.3;

% random polarization params
walk_amplitude = 0.1;
walk_stdev = pi/4;
walk_direction = rand(N, 1) * 2 * pi;
num_repolarization_steps = 10;
num_trailing_positions = 40;

init_repolarization_offset = floor(rand(N, 1) * num_repolarization_steps);

% Cucker-Smale params
CS_K = 30;
CS_Sigma = 1;
CS_Gamma = 0.8;
CS_NN_threshold = 5;
USE_NEAREST_NEIGHBORS = false;

% Curvature alignment params
num_neighbors = 1;                  % argument in knnsearch to find closest mesh pt
if(FORCE_CURVATURE_ALIGNMENT_ON)
    num_neighbors = 8;
end

% Supported modes:
% 'gauss-min' - align in direction of minimum Gaussian curvature
% 'gauss-max' - align in direction of maximum Gaussian curvature
% 'gauss-zero' - align in direction of lowest absolute Gaussian curvature
% 'mean-min' - align in direction of minimum mean curvature
% 'mean-max' - align in direction of maximum mean curvature
% 'mean-zero' - align in direction of lowest absolute mean curvature
alignment_mode = 'gauss-zero';
alignment_magnitude = 0.3;

% preallocate state variables
P = zeros(N, 3);
PV = zeros(N, 3);
q = [2, 5];
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

% use rejection sampling for initial position
% https://math.stackexchange.com/questions/2017079/uniform-random-points-on-a-torus
cnt = 0;
r = q(1);
R = q(2);
while cnt < N
    U = rand();
    V = rand();
    Theta = 2*pi*U;
    Phi = 2*pi*V;
    thresh = (R + r*cos(Theta))/(R + r);
    cnt = cnt + 1;
    X(cnt, :) = [(R + r*cos(Theta))*cos(Phi), (R + r*cos(Theta))*sin(Phi), r*sin(Theta)]; 
end

% preload pairwise geodesic distance between mesh points (for static surfaces)
if(isfile("torus_mesh.mat"))
    load("torus_mesh.mat");
else
    mesh_theta_num = 80;
    mesh_phi_num = 40;
    theta_grid = linspace(0, 2*pi, mesh_theta_num);
    phi_grid = linspace(0, 2*pi, mesh_phi_num);
    [Phi_mesh_fine, Theta_mesh_fine] = meshgrid(phi_grid, theta_grid); 
    mesh_x = (R+r.*cos(Theta_mesh_fine)).*cos(Phi_mesh_fine);
    mesh_y = (R+r.*cos(Theta_mesh_fine)).*sin(Phi_mesh_fine);
    mesh_z = r.*sin(Theta_mesh_fine);
    mat = adj_mat_torus(mesh_x,mesh_y,mesh_z);
    [dist_mat, next] = FloydWarshall(mat);
    save torus_mesh.mat Theta_mesh_fine Phi_mesh_fine mesh_theta_num mesh_phi_num mesh_x mesh_y mesh_z mat dist_mat next;
end
dist_range = [0 max(dist_mat(:))];

% create coarse mesh for visualization
theta_num = 36;
phi_num = 18;
theta_grid = linspace(0, 2*pi, theta_num);
phi_grid = linspace(0, 2*pi, phi_num);
[Phi_mesh, Theta_mesh] = meshgrid(phi_grid, theta_grid); 
vis_x = (R+r.*cos(Theta_mesh)).*cos(Phi_mesh);
vis_y = (R+r.*cos(Theta_mesh)).*sin(Phi_mesh);
vis_z = r.*sin(Theta_mesh);

% compute mean and gaussian curvature
G_curvature = gaussian_curvature_torus(Theta_mesh_fine, Phi_mesh_fine, q);
G_color_limits = [min(min(G_curvature)) max(max(G_curvature))];
M_curvature = mean_curvature_torus(Theta_mesh_fine, Phi_mesh_fine, q);
M_color_limits = [min(min(M_curvature)) max(max(M_curvature))];

% visualize IC
% visualize_surface(X, 0, vis_x, vis_y, vis_z, [-10 10], [-10 10], [-3 3]);
% visualize_geodesic_path(X, 0, pt_1_idx, pt_2_idx, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, mesh_phi_num, next, [-10 10], [-10 10], [-3 3]);
% visualize_geodesic_heatmap(X, 0, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, pt_1_idx, [-10 10], [-10 10], [-3 3], dist_range, dist_mat);
% visualize_curvature_heatmap(X, 0, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, [-10 10], [-10 10], [-3 3], G_color_limits, G_curvature, true);
visualize_trajectories(X, 0, prev_paths, path_colors, vis_x, vis_y, vis_z, [-10 10], [-10 10], [-3 3]);

t = 0;
itr = 0;

while t < totT
    
    % initialize nearest neighbors array
    [indices, dists] = all_mesh_neighbors(X, mesh_x, mesh_y, mesh_z, num_neighbors);

    % compute updated state vectors
    for i = 1 : N
        
        F(i) = (X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)^2 - 4*q(2)^2*(X(i,1)^2 + X(i,2)^2);

        dFdX_i_x = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*X(i,1) - 8*q(2)^2*X(i,1);
        dFdX_i_y = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*X(i,2) - 8*q(2)^2*X(i,2);
        dFdX_i_z = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*X(i,3);
        dFdX(i,:) = [dFdX_i_x, dFdX_i_y, dFdX_i_z];
        
        if (FORCE_EUCLIDEAN_REPULSION_ON)
            for j = setdiff(1:N, i)
                Fij = Alpha*exp(-1.0*(norm((X(i,:)-X(j,:)))^2)/(2*Sigma^2));
                deltaP = (X(i,:) - X(j,:))*Fij;
                P(i, :) = P(i, :) + deltaP;
            end
        end
        
        if (FORCE_ATTR_REPULSION_ON) 
            dPdt = 0;
            if(N > 1)
                for j = setdiff(1:N, i)                     % skip element i 
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
            P(i, :) = P(i, :) + deltaP;
        end
        
        if (FORCE_CUCKER_SMALE_POLARITY_ON)
            if itr > 0
                dPdt = [0 0 0];
                if(USE_NEAREST_NEIGHBORS)
                    particle_indices = find_neighbors(X, indices(:, 1), dists(:, 1), dist_mat, i, CS_NN_threshold);
                else
                    particle_indices = setdiff(1:N, i);
                end
                sz = numel(particle_indices);
                for j = 1:sz
                    idx = particle_indices(j);
                    dist = norm(X(i, :) - X(idx, :));
                    CS_H = CS_K/((CS_Sigma^2) + (dist^2))^CS_Gamma;
                    dPdt = dPdt + (1/sz) .* CS_H .* (PV(idx,:) - PV(i,:));
                end
                deltaP = deltaT * dPdt;
                P(i, :) = P(i, :) + deltaP;
            end
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

        dFdq_i_a = -4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*q(1);
        dFdq_i_R = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*q(2) - 8*q(2)*(X(i,1)^2 + X(i,2)^2); 
        dFdq(i,:) = [dFdq_i_a, dFdq_i_R];

        correction = (dot(dFdX(i,:), P(i,:)) + dot(dFdq(i,:), Q) + phi*F(i))/(norm(dFdX(i,:))^2);
        dXdt(i,:) = P(i,:) - correction*dFdX(i,:);
        
        % store trajectories
        if(itr < num_trailing_positions)
            prev_paths(i, itr + 1, :) = X(i, :);
        else
            prev_paths(i, 1:(num_trailing_positions - 1), :) = prev_paths(i, 2:num_trailing_positions, :);
            prev_paths(i, num_trailing_positions, :) = X(i, :);
        end
        
    end
    
    % update position and forces
    PV = dXdt;
    for i = 1 : N
        X(i,:) = X(i,:) + deltaT*dXdt(i,:);
    end
    
    P = zeros(N, 3);
    t = t + deltaT;
    itr = itr + 1;
    
    % visualize_surface(X, itr, vis_x, vis_y, vis_z, [-10 10], [-10 10], [-3 3]);
    % visualize_geodesic_path(X, itr, pt_1_idx, pt_2_idx, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, mesh_phi_num, next, [-10 10], [-10 10], [-3 3]);
    % visualize_geodesic_heatmap(X, itr, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, pt_1_idx, [-10 10], [-10 10], [-3 3], dist_range, dist_mat);
    % visualize_curvature_heatmap(X, itr, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, [-10 10], [-10 10], [-3 3], G_color_limits, G_curvature, true);
    visualize_trajectories(X, itr, prev_paths, path_colors, vis_x, vis_y, vis_z, [-10 10], [-10 10], [-3 3]);

end


function [adj_mat] = adj_mat_torus(x, y, z)

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
                new_j = mod(j+dx(k) - 1, width) + 1;
                distance = pdist([x(i,j) y(i,j) z(i,j) ; x(new_i, new_j) y(new_i, new_j) z(new_i, new_j)]);
                adj_mat((i-1)*width + j,(new_i - 1)*width + new_j) = distance;
            end
        end
    end
    
end

function [curvature] = gaussian_curvature_torus(Theta_mesh_fine, ~, q)

    % Adapted from: https://mathworld.wolfram.com/Torus.html
    r = q(1);
    R = q(2);
    num = cos(Theta_mesh_fine);
    denom = r * (R + (r * cos(Theta_mesh_fine)));
    curvature = num./denom;
    
end

function [curvature] = mean_curvature_torus(Theta_mesh_fine, ~ , q)
    
    r = q(1);
    R = q(2);
    num = R + (2 * r * cos(Theta_mesh_fine));
    denom = 2 * r * (R + (r * cos(Theta_mesh_fine)));
    curvature = num./denom;
    
end