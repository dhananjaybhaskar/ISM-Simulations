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
q = [20,15,20];
Q = [0.0, 0.0, 0.0];
F = zeros(N, 1);
dFdX = zeros(N, 3);
dFdq = zeros(N, 3);
dXdt = zeros(N, 3);

% pick random particle
pt_1_idx = floor(rand()*N) + 1;
pt_2_idx = floor(rand()*N) + 1;

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
    [Phi_mesh, Theta_mesh] = meshgrid(phi_grid, theta_grid); 
    mesh_x = a.*cos(Theta_mesh).*sin(Phi_mesh);
    mesh_y = b.*sin(Theta_mesh).*sin(Phi_mesh);
    mesh_z = c.*cos(Phi_mesh);
    mat = adj_mat_ellipsoid(mesh_x, mesh_y, mesh_z);
    [dist_mat, next] = FloydWarshall(mat);
    save ellipsoid_mesh.mat mesh_theta_num mesh_phi_num mesh_x mesh_y mesh_z mat dist_mat next;
end
dist_range = [0 max(dist_mat(:))];

% preload coarse mesh for visualization
theta_num = 36;
phi_num = 18;
theta_grid = linspace(0, 2*pi, theta_num);
phi_grid = linspace(0, pi, phi_num);
[Phi_mesh, Theta_mesh] = meshgrid(phi_grid, theta_grid); 
vis_x = a .*cos(Theta_mesh).*sin(Phi_mesh);
vis_y = b.*sin(Theta_mesh).*sin(Phi_mesh);
vis_z = c.*cos(Phi_mesh);
    
% visualize IC
visualize_geodesic_path(X, 0, pt_1_idx, pt_2_idx, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, phi_num, next, [-10 10], [-10 10], [-3 3]);
% visualize_geodesic_heatmap(X, 0, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, pt_1_idx, [-30 30], [-30 30], [-30 30], dist_range, dist_mat);

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

        F(i) = (X(i,1)^2/a^2) + (X(i, 2)^2/b^2) + (X(i,3)^2/c^2) - 1;

        dFdX_i_x = 2*X(i,1)/(a^2);
        dFdX_i_y = 2*X(i,2)/(b^2);
        dFdX_i_z = 2*X(i,3)/(c^2);
        dFdX(i,:) = [dFdX_i_x, dFdX_i_y, dFdX_i_z];

        dFdq_i_a = -2*X(i,1)^2/(a^3);
        dFdq_i_b = -2*X(i,2)^2/(b^3);
        dFdq_i_c = -2*X(i,3)^2/(c^3);
        dFdq(i,:) = [dFdq_i_a, dFdq_i_b, dFdq_i_c];

        correction = (dot(dFdX(i,:), P(i,:)) + dot(dFdq(i,:), Q) + phi*F(i))/(norm(dFdX(i,:))^2);
        dXdt(i,:) = P(i,:) - correction*dFdX(i,:);

    end
    
    % update position
    for i = 1 : N
        X(i,:) = X(i,:) + deltaT*dXdt(i,:);
    end
    
    t = t + deltaT;
    itr = itr + 1;
    
    visualize_geodesic_path(X, itr, pt_1_idx, pt_2_idx, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, phi_num, next, [-10 10], [-10 10], [-3 3]);
    % visualize_geodesic_heatmap(X, itr, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, pt_1_idx, [-30 30], [-30 30], [-30 30], dist_range, dist_mat);

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