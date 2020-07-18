%
% Agent-based model of particle movement on a genus 2 surface (double
% torus)
% Author: Tej Stead, Dhananjay Bhaskar
% Last Modified: Jul 14, 2020
% Reference: Using Particles to Sample and Control Implicit Surfaces
% Andrew P. Witkin and Paul S. Heckbert, Proc. SIGGRAPH '94
% Additional References:
% https://math.stackexchange.com/questions/152256/implicit-equation-for-double-torus-genus-2-orientable-surface
% http://stanwagon.com/wagon/mathimages/htmllinks/mathimages_1.html
%

% number of particles
N = 100;

% size of mesh
num_current_floaters = 100;
num_max_floaters = 20000;
% params
Alpha = 0.5;
Sigma = 0.15;
Sigma_floaters = Sigma*ones(num_max_floaters, 1);
phi = 0.1;
rho = 1;
beta = 0.1;
deltaT = 0.1;
totT = 10;

% init positions
X = zeros(N, 3);
X_floaters = zeros(num_max_floaters, 3);
% preallocate state variables
P = zeros(N, 3);
P_floaters = zeros(num_max_floaters, 3);
q = 0.2;
Q = 0.0;
F = zeros(N, 1);
F_floaters = zeros(num_max_floaters, 1);
D_floaters = zeros(num_max_floaters, 1);
dFdX = zeros(N, 3);
dFdq = zeros(N, 1);
dXdt = zeros(N, 3);
dFdX_floaters = zeros(num_max_floaters, 3);
dFdq_floaters = zeros(num_max_floaters, 1);
dXdt_floaters = inf*ones(num_max_floaters, 3); % not zero because of threshold

% pick a random particle
pt_1_idx = floor(rand()*N) + 1;
pt_2_idx = floor(rand()*N) + 1;

% use rejection sampling for initial position of particles
cnt = 0;
while cnt < N
    U = -1+2*rand();
    V = -1+2*rand();
    W = rand();
    det = -U^4 + 2*U^6 - U^8 + 2*U^2*V^2 - 2*U^4*V^2 - V^4 + q^2;
    if det > 0
        if W <= 0.5
            cnt = cnt + 1;
            X(cnt, :) = [U, V, sqrt(det)]; 
        else
            cnt = cnt + 1;
            X(cnt, :) = [U, V, -1*sqrt(det)];
        end
    end
end

% preload visualizer (for static surfaces)
if(isfile("double_torus_mesh.mat"))
    load("double_torus_mesh.mat")
else

% use rejection sampling for initial position of floaters
    [xs, ys] = meshgrid(-1.5:.05:1.5);
    zs = sqrt(-1.*xs.^4 + 2.*xs.^6 - xs.^8 + 2.*xs.^2.*ys.^2 - 2.*xs.^4.*ys.^2 - ys.^4 + q^2);
    zs(~(zs==real(zs))) = nan;
    cnt = 0;
    while cnt < num_current_floaters
        U = -1+2*rand();
        V = -1+2*rand();
        W = rand();
        det = -U^4 + 2*U^6 - U^8 + 2*U^2*V^2 - 2*U^4*V^2 - V^4 + q^2;
        if det > 0
            if W <= 0.5
                cnt = cnt + 1;
                X_floaters(cnt, :) = [U, V, sqrt(det)]; 
            else
                cnt = cnt + 1;
                X_floaters(cnt, :) = [U, V, -1*sqrt(det)];
            end
        end
    end
    % allow floaters to stabilize
%     thresh = Sigma_floaters * 0.00000005;
    is_stable = 0;
    itr_count = 0;
    E_ideal = 6 * exp(-2) * Alpha;
    while(is_stable == 0)
    % compute updated state vectors
        for i = 1 : num_current_floaters
            
            Di = 0;
            dDdSigma = 0;
            for j = 1 : num_current_floaters
                dist = norm((X_floaters(i,:)-X_floaters(j,:)));
                Fij = Alpha*exp(-1.0*dist^2/(2*Sigma_floaters(i)^2));
                Di = Di + Fij;
                Fji = Alpha*exp(-1.0*dist^2/(2*Sigma_floaters(i)^2));
                P_floaters(i,:) = P_floaters(i,:) + (X_floaters(i,:)-X_floaters(j,:))*Fij;
                P_floaters(i,:) = P_floaters(i,:) + (X_floaters(j,:)-X_floaters(i,:))*Fji ...
                    * (Sigma_floaters(i)/Sigma_floaters(j))^2;
                dDdSigma = dDdSigma + (dist^2 * Fij);
            end
            dDdt = -1 * rho * (Di - E_ideal);
            D_floaters(i) = Di;
            dDdSigma = dDdSigma/(Sigma_floaters(i)^3);
            F_floaters(i) = X_floaters(i,1)^8 - 2*X_floaters(i,1)^6 + X_floaters(i,1)^4 - 2*X_floaters(i,1)^2*X_floaters(i,2)^2 + 2*X_floaters(i,1)^4*X_floaters(i,2)^2 + X_floaters(i,2)^4 + X_floaters(i,3)^2 -q^2;

            dFdX_i_x = 8*X_floaters(i,1)^7 - 12*X_floaters(i,1)^5 + 4*X_floaters(i,1)^3 -4*X_floaters(i,1)*X_floaters(i,2)^2 + 8*X_floaters(i,1)^3*X_floaters(i,2)^2;
            dFdX_i_y = -4*X_floaters(i,1)^2*X_floaters(i,2) + 4*X_floaters(i,1)^4*X_floaters(i,2) + 4*X_floaters(i,2)^3;
            dFdX_i_z = 2*X_floaters(i,3);
            dFdX_floaters(i,:) = [dFdX_i_x, dFdX_i_y, dFdX_i_z];

            dFdq_floaters(i) = -2*q;

            correction = (dot(dFdX_floaters(i,:), P_floaters(i,:)) + dot(dFdq_floaters(i,:), Q) + phi*F_floaters(i))/(norm(dFdX_floaters(i,:))^2);
            dXdt_floaters(i,:) = P_floaters(i,:) - correction*dFdX_floaters(i,:);
            
            dSigmadt = dDdt/(dDdSigma + beta);
            Sigma_floaters(i) = Sigma_floaters(i) + deltaT * dSigmadt;
        end
        P_floaters = zeros(size(P_floaters));

        % update position
        for i = 1 : N
            X_floaters(i,:) = X_floaters(i,:) + deltaT*dXdt_floaters(i,:);
        end
        
%         % perform fission 
%         for i = 1:num_current_floaters
%             if((norm(dXdt_floaters(i)) < Sigma_floaters(i)/4) && (Sigma_floaters(i) > 2 * Sigma))
%                 Sigma_floaters(i) = Sigma_floaters(i)/sqrt(2);
%                 num_current_floaters = num_current_floaters + 1;
%                 X_floaters(num_current_floaters, :) = X_floaters(i, :);
%                 Sigma_floaters(num_current_floaters) = Sigma_floaters(i);
% %                 angles = [rand() rand() rand() rand()] * 2 * pi;
% %                 dir_1 = Sigma_floaters(i) * 0.05 *  [cos(angles(1))*sin(angles(2)) ...
% %                     sin(angles(1))*sin(angles(2)) cos(angles(2))];
% %                 dir_2 = Sigma_floaters(i) * 0.05 * [cos(angles(3))*sin(angles(4)) ...
% %                     sin(angles(3))*sin(angles(4)) cos(angles(4))];
% %                 P_floaters(i, :) = dir_1;
% %                 P_floaters(num_current_floaters, :) = dir_2;
%             end
%         end
%         
%         %perform death
%         kill_points = zeros(num_current_floaters, 1);
%         kill_count = 0;
%         for i = 1:num_current_floaters
%             if((norm(dXdt_floaters(i)) < Sigma_floaters(i)/4) && (Sigma_floaters(i) < 0.4 * Sigma))
%                 temp = rand();
%                 if(temp > Sigma_floaters(i)/(0.4 * Sigma))
%                     kill_count = kill_count + 1;
%                     kill_points(kill_count) = i;
%                 end
%             end
%         end
%         kill_points = kill_points(kill_points ~= 0);
%         num_current_floaters = num_current_floaters - numel(kill_points);
%         for i = 1: numel(kill_points)
%             kill_idx = kill_points(i);
%             X_floaters(kill_idx, :) = [];
%             Sigma_floaters(kill_idx) = [];
%             P_floaters(kill_idx, :) = [];
%         end
        for i = 1:num_current_floaters
            if(norm(X_floaters(i, :)) > 3)
                X_floaters(i, :)
                i
                Sigma_floaters(i)
            end
        end
        figure(1); clf; scatter3(X_floaters(:, 1), X_floaters(:, 2), X_floaters(:, 3)); pause;
    end
   build_hexagonal_mesh(X_floaters); pause; 
end

% visualize_genus_2_mesh(X, q, 0);

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

        F(i) = X(i,1)^8 - 2*X(i,1)^6 + X(i,1)^4 - 2*X(i,1)^2*X(i,2)^2 + 2*X(i,1)^4*X(i,2)^2 + X(i,2)^4 + X(i,3)^2 -q^2;

        dFdX_i_x = 8*X(i,1)^7 - 12*X(i,1)^5 + 4*X(i,1)^3 -4*X(i,1)*X(i,2)^2 + 8*X(i,1)^3*X(i,2)^2;
        dFdX_i_y = -4*X(i,1)^2*X(i,2) + 4*X(i,1)^4*X(i,2) + 4*X(i,2)^3;
        dFdX_i_z = 2*X(i,3);
        dFdX(i,:) = [dFdX_i_x, dFdX_i_y, dFdX_i_z];
        
        dFdq(i) = -2*q;

        correction = (dot(dFdX(i,:), P(i,:)) + dot(dFdq(i,:), Q) + phi*F(i))/(norm(dFdX(i,:))^2);
        dXdt(i,:) = P(i,:) - correction*dFdX(i,:);

    end
    
    % update position
    for i = 1 : N
        X(i,:) = X(i,:) + deltaT*dXdt(i,:);
    end
    
    t = t + deltaT;
    itr = itr + 1;
    visualize_genus_2_mesh(X, q, itr);
    
end

function [] = visualize_genus_2_mesh(X, q, itr)

%     fig = figure('visible', 'off');
    figure(1);
    [xs, ys] = meshgrid([-1.2:.01:1.2], [-0.68:0.001:-0.651 -0.65:0.01:0.65 0.651:0.001:0.68]);
    zs = sqrt(-1.*xs.^4 + 2.*xs.^6 - xs.^8 + 2.*xs.^2.*ys.^2 - 2.*xs.^4.*ys.^2 - ys.^4 + q^2);
    zs(~(zs==real(zs))) = nan;
    surf(xs, ys, zs, 'edgecolor', 'none', 'facecolor', 'interp');
    hold on;
    surf(xs, ys, -zs, 'edgecolor', 'none', 'facecolor', 'interp');
    pause;
    alpha 0.5;
    [xs, ys] = meshgrid(-1.2:.06:1.2);
    zs = sqrt(-1.*xs.^4 + 2.*xs.^6 - xs.^8 + 2.*xs.^2.*ys.^2 - 2.*xs.^4.*ys.^2 - ys.^4 + q^2);
    zs(~(zs==real(zs))) = nan;
    surf(xs, ys, zs, 'edgecolor', 'k', 'facecolor', 'none', 'edgealpha', 1);
    surf(xs, ys, -zs, 'edgecolor', 'k', 'facecolor', 'none', 'edgealpha', 1);
    %daspect([1 1 1])
    hold on;
    scatter3(X(:,1), X(:,2), X(:,3), 14, 'filled');
    xlim([-1.2 1.2])
    ylim([-1.2 1.2])
    zlim([-0.3 0.3])
    pause;
%     fname = strcat('sim_', sprintf('%03d',itr), '.png');
%     saveas(fig, fname, 'png');

end