%
% Agent-based model of constrained particle motion on an ellipsoid
% Author: Tej Stead
% Last Modified: Jul 10, 2020
%

%
% TODO: simulate simple repulsion model on implicit surface
% TODO: compute geodesic distance between particles
% TODO: plot surface curvature
%

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

pt_1_idx = floor(rand()*N) + 1;
pt_2_idx = floor(rand()*N) + 1;

%uniform distribution of (Theta, Phi) in [0,2pi] for initial position
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
    X(cnt, :) = [a*cos(Theta)*sin(Phi), b*sin(Theta)*sin(Phi), c*cos(Phi)]; %ellipsoid parametrization
end

%preload visualizer (for static surfaces)
    theta_num = 36;
    phi_num = 18;
    theta_grid = linspace(0, 2*pi, theta_num);
    phi_grid = linspace(0, pi, phi_num);
    [Phi_mesh, Theta_mesh] = meshgrid(phi_grid, theta_grid); 
    x = a .*cos(Theta_mesh).*sin(Phi_mesh);
    y = b.*sin(Theta_mesh).*sin(Phi_mesh);
    z = c.*cos(Phi_mesh);
    prev_x = x;
    mat = adj_mat(x,y,z);
    [dist_mat, next] = FloydWarshall(mat);
    
% visualize
visualize(X, 0, pt_1_idx, pt_2_idx,x,y,z,phi_num, next, [-30 30], [-30 30], [-30 30]);
t = 0;
itr = 0;

while t < totT

    % compute updated state vectors
    for i = 1 : N

        P(i,:) = [0, 0, 0]; 
        for j = 1 : N
            Fij = Alpha*exp(-1.0*norm((X(i,:)-X(j,:)))/(2*Sigma^2));
            P(i,:) = P(i,:) + (X(i,:) - X(j,:))*Fij; %directional energy
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
    visualize(X, itr, pt_1_idx, pt_2_idx,x,y,z,phi_num, next, [-30 30], [-30 30], [-30 30]);
end

function [adj_mat] = adj_mat(x,y,z)
    sz = size(x);
    height = sz(1);
    width = sz(2);
    adj_mat = inf*ones(height*width, height*width);
    for i = 1:height
        for j = 1:width
            dx = [-1 -1 -1 0 0 0 1 1 1];
            dy = [-1 0 1 -1 0 1 -1 0 1];
%               dx = [-1 0 0 0 1];
%               dy = [0 -1 0 1 0];
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