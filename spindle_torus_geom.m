%
% Agent-based model of particle movement on implicit surface (torus)
% Authors: Dhananjay Bhaskar, Tej Stead
% Last Modified: Jul 14, 2020
% Reference: Using Particles to Sample and Control Implicit Surfaces
% Andrew P. Witkin and Paul S. Heckbert, Proc. SIGGRAPH '94
%

% seed RNG
rng(1337)

% number of particles
N = 200;

% params
Alpha = 2;
Sigma = 0.5;
phi = 1;
deltaT = 0.1;
totT = 5;

% init positions
X = zeros(N, 3);

% preallocate state variables
P = zeros(N, 3);
q = [5, 2];
Q = [0.0, 0.0];
F = zeros(N, 1);
dFdX = zeros(N, 3);
dFdq = zeros(N, 2);
dXdt = zeros(N, 3);

% pick a random particle
pt_1_idx = floor(rand()*N) + 1;
pt_2_idx = floor(rand()*N) + 1;

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
if(isfile("spindle_torus_mesh.mat"))
    load("spindle_torus_mesh.mat");
else
    mesh_theta_num = 80;
    mesh_phi_num = 40;
    theta_grid = linspace(0, 2*pi, mesh_theta_num);
    phi_grid = linspace(0, 2*pi, mesh_phi_num);
    [Phi_mesh_fine, Theta_mesh_fine] = meshgrid(phi_grid, theta_grid); 
    mesh_x = (R+r.*cos(Theta_mesh_fine)).*cos(Phi_mesh_fine);
    mesh_y = (R+r.*cos(Theta_mesh_fine)).*sin(Phi_mesh_fine);
    mesh_z = r.*sin(Theta_mesh_fine);
    mat = adj_mat_spindle_torus(mesh_x,mesh_y,mesh_z);
    [dist_mat, next] = FloydWarshall(mat);
    save spindle_torus_mesh.mat Theta_mesh_fine Phi_mesh_fine mesh_theta_num mesh_phi_num mesh_x mesh_y mesh_z mat dist_mat next;
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
G_curvature = gaussian_curvature_spindle_torus(Theta_mesh_fine, Phi_mesh_fine, q);
G_color_limits = [min(min(G_curvature)) max(max(G_curvature))];
M_curvature = mean_curvature_torus(Theta_mesh_fine, Phi_mesh_fine, q);
M_color_limits = [min(min(M_curvature)) max(max(M_curvature))];
% visualize_curvature_heatmap(X,0,vis_x,vis_y,vis_z,mesh_x,mesh_y,mesh_z, [-10 10], [-10 10], [-10 10], G_color_limits, G_curvature);
% visualize_curvature_heatmap(X,1,vis_x,vis_y,vis_z,mesh_x,mesh_y,mesh_z, [-10 10], [-10 10], [-10 10], M_color_limits, M_curvature);
visualize_geodesic_path(X, 0, pt_1_idx, pt_2_idx, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, mesh_phi_num, next, [-10 10], [-10 10], [-10 10]);
% visualize_geodesic_heatmap(X, 0, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, pt_1_idx, [-10 10], [-10 10], [-3 3], dist_range, dist_mat);

t = 0;
itr = 0;

while t < totT

    % compute updated state vectors
    for i = 1 : N

        P(i,:) = [0, 0, 0]; 
        for j = 1 : N
            Fij = Alpha*exp(-1.0*norm((X(i,:)-X(j,:)))/(2*Sigma^2));
            P(i,:) = P(i,:) + (X(i,:) - X(j,:))*Fij;
        end

        F(i) = (X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)^2 - 4*q(2)^2*(X(i,1)^2 + X(i,2)^2);

        dFdX_i_x = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*X(i,1) - 8*q(2)^2*X(i,1);
        dFdX_i_y = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*X(i,2) - 8*q(2)^2*X(i,2);
        dFdX_i_z = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*X(i,3);
        dFdX(i,:) = [dFdX_i_x, dFdX_i_y, dFdX_i_z];

        dFdq_i_a = -4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*q(1);
        dFdq_i_R = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*q(2) - 8*q(2)*(X(i,1)^2 + X(i,2)^2); 
        dFdq(i,:) = [dFdq_i_a, dFdq_i_R];

        correction = (dot(dFdX(i,:), P(i,:)) + dot(dFdq(i,:), Q) + phi*F(i))/(norm(dFdX(i,:))^2);
        dXdt(i,:) = P(i,:) - correction*dFdX(i,:);

    end
    
    % update position
    for i = 1 : N
        X(i,:) = X(i,:) + deltaT*dXdt(i,:);
    end
    
    t = t + deltaT;
    itr = itr + 1;
    
    visualize_geodesic_path(X, itr, pt_1_idx, pt_2_idx, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, mesh_phi_num, next, [-10 10], [-10 10], [-10 10]);
    % visualize_geodesic_heatmap(X, itr, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, pt_1_idx, [-10 10], [-10 10], [-3 3], dist_range, dist_mat);
    
end


function [adj_mat] = adj_mat_spindle_torus(x, y, z)

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

function [curvature] = gaussian_curvature_spindle_torus(Theta_mesh_fine, ~, q)
    % Second argument is not needed but is taken for consistency with other
    % surfaces
    % https://mathworld.wolfram.com/Torus.html, adapted
    r = q(1);
    R = q(2);
    num = cos(Theta_mesh_fine);
    denom = r * (R + (r * cos(Theta_mesh_fine)));
    curvature = num ./ denom;
end

function [curvature] = mean_curvature_spindle_torus(Theta_mesh_fine, ~ , q)
    % Second argument is not needed but is taken for consistency with other
    % surfaces
    % https://mathworld.wolfram.com/Torus.html, adapted
    r = q(1);
    R = q(2);
    num = R + (2 * r * cos(Theta_mesh_fine));
    denom = 2 * r * (R + (r * cos(Theta_mesh_fine)));
    curvature = num ./ denom;
end