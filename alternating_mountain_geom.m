%
% Agent-based model of constrained particle motion on an ellipsoid
% Authors: Tej Stead, Dhananjay Bhaskar
% Last Modified: Jul 10, 2020
%

% seed RNG
rng(1337);

% number of particles
N = 50;

% params
Alpha = 2;
Sigma = 1.5;
phi = 1;
deltaT = 0.1;
totT = 5;

% init positions
X = zeros(N, 3);

% preallocate state variables
P = zeros(N, 3);
q = [2,4];
Q = [0.0, 0.0];
F = zeros(N, 1);
dFdX = zeros(N, 3);
dFdq = zeros(N, 2);
dXdt = zeros(N, 3);

% pick random particle
pt_1_idx = floor(rand()*N) + 1;
pt_2_idx = floor(rand()*N) + 1;

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
if(isfile("alternating_mountain_mesh_2.mat"))
    load("alternating_mountain_mesh_2.mat");
else
    mesh_theta_num = 80;
    mesh_phi_num = 40;
    theta_grid = linspace(0, 6*pi*a, mesh_theta_num);
    phi_grid = linspace(0, 6*pi*a, mesh_phi_num);
    [Phi_mesh_fine, Theta_mesh_fine] = meshgrid(phi_grid, theta_grid); 
    mesh_x = Theta_mesh_fine;
    mesh_y = Phi_mesh_fine;
    mesh_z = b .* sin(Theta_mesh_fine ./ a) .* sin(Phi_mesh_fine ./ a);
    mat = adj_mat_alternating_mountain(mesh_x, mesh_y, mesh_z);
    [dist_mat, next] = FloydWarshall(mat);
    save alternating_mountain_mesh_2.mat mesh_theta_num mesh_phi_num mesh_x mesh_y ...
        Phi_mesh_fine Theta_mesh_fine mesh_z mat dist_mat next;
end
dist_range = [0 max(dist_mat(:))];

% preload coarse mesh for visualization
theta_num = 36;
phi_num = 18;
theta_grid = linspace(0, 6*pi*a, theta_num);
phi_grid = linspace(0, 6*pi*a, phi_num);
[Phi_mesh, Theta_mesh] = meshgrid(phi_grid, theta_grid); 
vis_x = Theta_mesh;
vis_y = Phi_mesh;
vis_z = b .* sin(Theta_mesh ./ a) .* sin(Phi_mesh ./ a);
    
% visualize IC
G_curvature = gaussian_curvature_alternating_mountain(Theta_mesh_fine, Phi_mesh_fine, q);
G_color_limits = [0 max(max(G_curvature))];
M_curvature = mean_curvature_alternating_mountain(Theta_mesh_fine, Phi_mesh_fine, q);
M_color_limits = [min(min(M_curvature)) max(max(M_curvature))];
% visualize_curvature_heatmap(X,0,vis_x,vis_y,vis_z,mesh_x,mesh_y,mesh_z, [-30 30], [-30 30], [-30 30], G_color_limits, G_curvature);
% visualize_curvature_heatmap(X,1,vis_x,vis_y,vis_z,mesh_x,mesh_y,mesh_z, [-30 30], [-30 30], [-30 30], M_color_limits, M_curvature);
% pause;
visualize_geodesic_path(X, 0, pt_1_idx, pt_2_idx, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, mesh_phi_num, next, [-10 20*pi], [-10 20*pi], [-10 10]);
% visualize_geodesic_heatmap(X, 0, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, pt_1_idx, [-30 30], [-30 30], [-30 30], dist_range, dist_mat);

t = 0;
itr = 0;

while t < totT

    % compute updated state vectors
    for i = 1 : N

        P(i,:) = [0, 0, 0]; 
        for j = 1 : N
            Fij = Alpha*exp(-1.0*(norm((X(i,:)-X(j,:)))^2)/(2*Sigma^2));
            P(i,:) = P(i,:) + (X(i,:) - X(j,:))*Fij;
        end

        F(i) = (X(i,3)) - b * sin(X(i,2)/a) * sin(X(i, 1)/a);

        dFdX_i_x = -1 * (b/a) * sin(X(i,2)/a) * cos(X(i,1)/a);
        dFdX_i_y = -1 * (b/a) * sin(X(i,1)/a) * cos(X(i,2)/a);
        dFdX_i_z = 1;
        dFdX(i,:) = [dFdX_i_x, dFdX_i_y, dFdX_i_z];

        dFdq_i_a = (b/a^2) * ((X(i,2)*sin(X(i,1)/a)*cos(X(i,2)/a)) + X(i,1)*cos(X(i,1)/a)*sin(X(i,2)/a));
        dFdq_i_b = -1 * sin(X(i,1)/a) * sin(X(i,2)/a);
        dFdq(i,:) = [dFdq_i_a, dFdq_i_b];

        correction = (dot(dFdX(i,:), P(i,:)) + dot(dFdq(i,:), Q) + phi*F(i))/(norm(dFdX(i,:))^2);
        dXdt(i,:) = P(i,:) - correction*dFdX(i,:);

    end
    
    % update position
    for i = 1 : N
        X(i,:) = X(i,:) + deltaT*dXdt(i,:);
    end
    
    t = t + deltaT;
    itr = itr + 1;
    
    visualize_geodesic_path(X, itr, pt_1_idx, pt_2_idx, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, mesh_phi_num, next, [-10 20*pi], [-10 20*pi], [-10 10]);
    % visualize_geodesic_heatmap(X, itr, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, pt_1_idx, [-30 30], [-30 30], [-30 30], dist_range, dist_mat);

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

function [curvature] = gaussian_curvature_alternating_mountain(Theta_mesh_fine, Phi_mesh_fine, q)
    curvature = 0;
end

function [curvature] = mean_curvature_alternating_mountain(Theta_mesh_fine, Phi_mesh_fine, q)
    curvature = 0;
end