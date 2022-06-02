close all; clear;
% Create base attachment points

Br = 300;               % Base attachment radius (mm)
Pr = 400;               % Platform attachment radius (mm)
delta_p = 6 * pi/180;   % Spacing between attachment points (rad) 
delta_b = 24 * pi/180;   % Spacing between attachment points (rad) 


[k] = (1:6); % Actuators
phi_pk = 2*pi/3 * floor(k./2) - (-1).^k * delta_p /2 + pi/3; % Angular placement of actuators
phi_bk = 2*pi/3 * floor((k+1)./2) + (-1).^k * delta_b /2;

p_k = Pr * [cos(phi_pk') sin(phi_pk') zeros(6,1)]; % Vector placement of actuators
b_k = Br * [cos(phi_bk') sin(phi_bk') zeros(6,1)];

%% Calculate actuator lengths based on platform position and rotation
t_end = 18;
t = linspace(0, t_end, 500)';      % Time
Tt = [0+200*sin(2*pi*1/3.*t).*(t<1*t_end/6) 0+200*sin(2*pi*1/3.*t).*(t>1*t_end/6 & t<2*t_end/6) 500+200*sin(2*pi*1/3.*t).*(t>2*t_end/6 & t<3*t_end/6)];    % Translation [x y z] as function of time
alpha = 20 * pi/180;%; * sin(2*pi*1*t); % Rotation angle as function of time
v_rot = [sin(2*pi*1/3.*t).*(t>3*t_end/6 & t<4*t_end/6) sin(2*pi*1/3.*t).*(t>4*t_end/6 & t<5*t_end/6) sin(2*pi*1/3.*t).*(t>5*t_end/6)];    % Rotation direction vector
%v_rot = v_rot ./ norm(v_rot); % Normalize
Rt = [cos(alpha/2)+0.*v_rot(:,1) sin(alpha/2).*v_rot]; % Rotation as quaternions

% Length of actuator at point in time
lk = zeros(size(t, 1), 6);
for n = 1:size(t, 1)
    % Calculate for each timestep
    T = Tt(n,:);
    quat = quaternion(Rt(n,:));
    lk(n,:) = vecnorm( (T + rotatepoint(quat, p_k) - b_k), 2, 2);
end

% Plot
hold on     
plot(t, lk(:, 1))
plot(t, lk(:, 2))
plot(t, lk(:, 3))
plot(t, lk(:, 4))
plot(t, lk(:, 5))
plot(t, lk(:, 6))
hold off

%% Plot fun
figure
k = 1;
f = getframe(gcf);
for n = 1:size(t, 1)
    quat = quaternion(Rt(n,:));
    p_t = rotatepoint(quat, p_k) + Tt(n, :);
    circ_alpha = linspace(0, 2*pi, 1000)';
    base_circ = [Br*cos(circ_alpha), Br*sin(circ_alpha), 0.*circ_alpha];
    plat_circ = [Pr*cos(circ_alpha), Pr*sin(circ_alpha), 0.*circ_alpha];
    plat_circ = Tt(n, :) + rotatepoint(quat, plat_circ);
    X(1, :) = b_k(:,1);
    X(2, :) = p_t(:,1);
    Y(1, :) = b_k(:,2);
    Y(2, :) = p_t(:,2);
    Z(1, :) = b_k(:,3);
    Z(2, :) = p_t(:,3);
    p = plot3(X, Y, Z, " o", ...
            [X(1, 1) X(2, 1)], [Y(1, 1) Y(2, 1)], [Z(1, 1) Z(2, 1)], "-r", ...
            [X(1, 2) X(2, 2)], [Y(1, 2) Y(2, 2)], [Z(1, 2) Z(2, 2)], "-r", ...
            [X(1, 3) X(2, 3)], [Y(1, 3) Y(2, 3)], [Z(1, 3) Z(2, 3)], "-r", ...
            [X(1, 4) X(2, 4)], [Y(1, 4) Y(2, 4)], [Z(1, 4) Z(2, 4)], "-r", ...
            [X(1, 5) X(2, 5)], [Y(1, 5) Y(2, 5)], [Z(1, 5) Z(2, 5)], "-r", ...
            [X(1, 6) X(2, 6)], [Y(1, 6) Y(2, 6)], [Z(1, 6) Z(2, 6)], "-r", ...
            base_circ(:,1), base_circ(:,2), base_circ(:,3), "-b", ...
            plat_circ(:,1), plat_circ(:,2), plat_circ(:,3), "-b");
     xlim([-600, 600])
     ylim([-600, 600])
     zlim([0, 800])
     for i = 7:12
        p(i).LineWidth = 3; 
     end
     f = getframe(gcf);
     if k == 1 
        [im, map] = rgb2ind(f.cdata, 256, 'nodither');
     else
        im(:,:,1,k) = rgb2ind(f.cdata, map, 'nodither');
     end
     k = k+1
end

deltaT = t(2)-t(1);
imwrite(im, map, 'platformani.gif', 'DelayTime', deltaT, 'LoopCount', inf)

%% 
