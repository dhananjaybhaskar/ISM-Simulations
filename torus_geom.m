%
% Agent-based model of particle movement on implicit surface (torus)
% Authors: Dhananjay Bhaskar, Tej Stead
% Last Modified: Jul 14, 2020
% Reference: Using Particles to Sample and Control Implicit Surfaces
% Andrew P. Witkin and Paul S. Heckbert, Proc. SIGGRAPH '94
%

% seed RNG
rng(1337);
% number of particles
N = 70; %testing for single cell repolarization

% params
Alpha = 2;
Sigma = 0.2;
phi = 1;
deltaT = 0.1;
totT = 10;
num_repolarization_steps = 10;
num_trailing_positions = 100;
% init positions
X = zeros(N, 3);
initial_repolarization_times = floor(rand(N, 1) * 10);

% preallocate state variables
P = zeros(N, 3);
q = [2, 5];
Q = [0.0, 0.0];
F = zeros(N, 1);
dFdX = zeros(N, 3);
dFdq = zeros(N, 2);
dXdt = zeros(N, 3);

% pick a random particle
pt_1_idx = floor(rand()*N) + 1;
pt_2_idx = floor(rand()*N) + 1;

walk_amplitude = 0.3;
walk_stdev = pi/4;
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

% visualize IC
% visualize_surface(X, 0, vis_x, vis_y, vis_z, [-10 10], [-10 10], [-3 3]);
G_curvature = gaussian_curvature_torus(Theta_mesh_fine, Phi_mesh_fine, q);
G_color_limits = [min(min(G_curvature)) max(max(G_curvature))];
M_curvature = mean_curvature_torus(Theta_mesh_fine, Phi_mesh_fine, q);
M_color_limits = [min(min(M_curvature)) max(max(M_curvature))];
% visualize_curvature_heatmap(X,0,vis_x,vis_y,vis_z,mesh_x,mesh_y,mesh_z, [-10 10], [-10 10], [-10 10], G_color_limits, G_curvature);
% visualize_curvature_heatmap(X,1,vis_x,vis_y,vis_z,mesh_x,mesh_y,mesh_z, [-10 10], [-10 10], [-10 10], M_color_limits, M_curvature);
% visualize_geodesic_path(X, 0, pt_1_idx, pt_2_idx, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, mesh_phi_num, next, [-10 10], [-10 10], [-10 10]);
% visualize_geodesic_heatmap(X, 0, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, pt_1_idx, [-10 10], [-10 10], [-3 3], dist_range, dist_mat);
prev_paths = nan * ones(N,20,3);
path_colors = hsv(N);

%randomize color layout
for i = 1:N
    j = floor(rand() * N) + 1;
    temp = path_colors(j, :);
    path_colors(j, :) = path_colors(i, :);
    path_colors(i, :) = temp;
end
visualize_particle_path(X,0, prev_paths, path_colors, vis_x,vis_y,vis_z, [-10 10], [-10 10], [-10 10]);
t = 0;
itr = 0;

walk_direction = rand(N, 1) * 2 * pi;
while t < totT

    % compute updated state vectors
    for i = 1 : N

        P(i,:) = [0, 0, 0];

        F(i) = (X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)^2 - 4*q(2)^2*(X(i,1)^2 + X(i,2)^2);

        dFdX_i_x = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*X(i,1) - 8*q(2)^2*X(i,1);
        dFdX_i_y = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*X(i,2) - 8*q(2)^2*X(i,2);
        dFdX_i_z = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*X(i,3);
        dFdX(i,:) = [dFdX_i_x, dFdX_i_y, dFdX_i_z];
        
        % generate orthonormal basis for tangent plane 
        nullspace = [dFdX(i,:)' zeros(3,2)];
        assert(numel(nullspace) == 9);
        nullspace = null(nullspace);
        nullspace = nullspace';
        
        if(1) %turn back on to enable simple repulsion
            for j = 1 : N %calculate repulsion force
                Fij = Alpha*exp(-1.0*(norm((X(i,:)-X(j,:)))^2)/(2*Sigma^2));
                P(i,:) = P(i,:) + (X(i,:) - X(j,:))*Fij;
            end
        end
        
        
        if(1) % turn off to disable random motion
            if(mod(itr, num_repolarization_steps) == initial_repolarization_times(i))
                temp = rand();
                walk_direction(i) = walk_direction(i) + norminv(temp,0,walk_stdev);
            end
            P(i,:) = P(i,:) + cos(walk_direction(i)) * nullspace(1, :) * walk_amplitude;
            P(i,:) = P(i,:) + sin(walk_direction(i)) * nullspace(2, :) * walk_amplitude;
        end

        dFdq_i_a = -4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*q(1);
        dFdq_i_R = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*q(2) - 8*q(2)*(X(i,1)^2 + X(i,2)^2); 
        dFdq(i,:) = [dFdq_i_a, dFdq_i_R];

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
    
    t = t + deltaT;
    itr = itr + 1;
    
    visualize_particle_path(X,itr, prev_paths, path_colors, vis_x,vis_y,vis_z, [-10 10], [-10 10], [-10 10]);
%     visualize_geodesic_path(X, itr, pt_1_idx, pt_2_idx, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, mesh_phi_num, next, [-10 10], [-10 10], [-10 10]);
    % visualize_geodesic_heatmap(X, itr, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, pt_1_idx, [-10 10], [-10 10], [-3 3], dist_range, dist_mat);
    
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
    % Second argument is not needed but is taken for consistency with other
    % surfaces
    % https://mathworld.wolfram.com/Torus.html, adapted
    r = q(1);
    R = q(2);
    num = cos(Theta_mesh_fine);
    denom = r * (R + (r * cos(Theta_mesh_fine)));
    curvature = num ./ denom;
end

function [curvature] = mean_curvature_torus(Theta_mesh_fine, ~ , q)
    % Second argument is not needed but is taken for consistency with other
    % surfaces
    % https://mathworld.wolfram.com/Torus.html, adapted
    r = q(1);
    R = q(2);
    num = R + (2 * r * cos(Theta_mesh_fine));
    denom = 2 * r * (R + (r * cos(Theta_mesh_fine)));
    curvature = num ./ denom;
end