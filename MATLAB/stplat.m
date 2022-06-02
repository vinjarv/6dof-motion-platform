classdef stplat < handle
    properties
        % Notation and equations based on Harib, Khalifa & Srinivasan, Krishnaswamy. (2003)
        % Kinematic and dynamic analysis of Stewart platform-based machine tool structures. 
        % Robotica. 21. 541-554. 10.1017/S0263574703005046. 
        % https://www.researchgate.net/publication/220103568_Kinematic_and_dynamic_analysis_of_Stewart_platform-based_machine_tool_structures

        Br; Pr; % Base and platform radii
        Bk; Pk; % Vectors to actuator attachment points
        Act_fixed_l; Act_max_l; % 
        Act_m1; Act_m2; % Mass components of actuators
        Act_l1; Act_l2  % COM locations of actuators 
        Act_Inn1; Act_Inn2; % Actuator mass moments of inertia, radial, for fixed and moving parts
        
        Plat_m; % Mass of platform
        Plat_I; % Moment of inertia tensor of platform
        Plat_COM; % Center of mass vector
        
        home_pos; % Default position of actuator
        
    end
    methods
        function obj = stplat()
            % Constructor
            % TODO: values in function call?
            obj.Br = 1.700 / 2;           % Base attachment radius (m)
            obj.Pr = 0.989;               % Platform attachment radius (m)
            delta_p = 8.45 * pi/180;   % Spacing between attachment points (rad) 
            delta_b = 16 * pi/180;  % Spacing between attachment points (rad) 
            
            obj.Act_fixed_l = 1.125;      % Fixed part of actuator (m)
            obj.Act_max_l = obj.Act_fixed_l+0.380;       % Max length of actuator (attachment point to attachment point) (m)
            obj.home_pos = [0 1.0763 0 0 -pi/2 0]'; % Default position of platform in world frame
            
            % Mass of actuators
            motor_weight = 3.9;
            rho_steel = 8000; %kg/m3
            ballscrew_dia = 0.020;
            ballscrew_weight = (ballscrew_dia/2)^2 *obj.Act_fixed_l * rho_steel;
            obj.Act_m1 = motor_weight + ballscrew_weight;
            obj.Act_m2 = ballscrew_weight; % very rough estimate
            
            % Inertia of actuator parts
            motor_diameter = 0.095;
            motor_length = 0.140;
            motor_inertia = motor_weight * (1/4*(motor_diameter/2)^2 + 1/12*motor_length^2);
            ballscrew_length = obj.Act_fixed_l;
            ballscrew_inertia = ballscrew_weight * 1/12 * ballscrew_length^2;

            % Center of mass
            obj.Act_l1 = (motor_length/2*motor_weight + (motor_length+ballscrew_length/2)*ballscrew_weight) / obj.Act_m1;
            % Assuming upper part is one fixed part, no translation of
            % rotation axis needed
            obj.Act_l2 = ballscrew_length;
            
            % Actuator inertia around COM
            obj.Act_Inn1 = motor_inertia + motor_weight*(obj.Act_l1 - motor_length/2)^2 ...
                           + ballscrew_inertia + ballscrew_weight*(obj.Act_l1 - (motor_length + ballscrew_length/2))^2;
            obj.Act_Inn2 = obj.Act_m2 * (1/4*(0.100/2)^2 + 1/12*(obj.Act_l2/2)^2); % Assume Ø100
            
            % Platform and payload mass/inertia
            obj.Plat_m = 100; % Total mass
            obj.Plat_I = [-4.98 2.45 0.14; -0.14 0.07 9.24; -2.45 5.94 0.07];
            obj.Plat_COM = [-0.200 0 0.2]';
            
            [k] = (1:6); % Actuators
            phi_pk = 2*pi/3 * floor(k./2) - (-1).^k * delta_p /2 + pi/3; % Angular placement of actuators
            phi_bk = 2*pi/3 * floor((k+1)./2) + (-1).^k * delta_b /2;

            obj.Pk = obj.Pr * [cos(phi_pk') sin(phi_pk') zeros(6,1)]'; % Vector placement of actuators
            obj.Bk = obj.Br * [cos(phi_bk') sin(phi_bk') zeros(6,1)]'; 
            obj.Bk = rotx(-90)*obj.Bk;   % In X-Z plane to avoid singularities
        end
        
        function [l, n, a_w, L, wRp] = ikine(obj, q)
            % Solve inverse kinematics for a pose
            % Returns:
            % l - scalar lengths of actuators
            % n - normalized actuator vectors
            % a_w - attachment point coordinates in world frame
            % L - complete length vectors of actuators
            % wRp - rotation matrix from platform to world coordinates
            % Input:
            % q - pose vector - [x, y, z, phi, theta, psi]'
            x = q(1:3);         % Get translational part of pose
            angles = q(4:6);    % Get rotational part of pose
            wRp = obj.wRp(angles); % Rotation matrix
            
            l = zeros(6,1);n=zeros(3,6);a_w = zeros(3,6);L=zeros(3,6);
            for i = 1:6
                a_w(:,i) = x + wRp*obj.Pk(:,i);     % Translate and rotate attachment point around tool frame
                L(:,i) = a_w(:,i) - obj.Bk(:,i);    % Vector sum
                l(i) = vecnorm(L(:,i));             % Pythagoras
                n(:,i)=L(:,i)/vecnorm(L(:,i));      % Normalized length vectors
            end
        end
        
        function [wRp] = wRp(~, angles)
            % Get rotation matrix from platform to world
            phi = angles(1); theta = angles(2); psi = angles(3);
%             wRpold = [ cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi) -sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi) sin(theta)*sin(phi)
%                     cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi) -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi) -sin(theta)*cos(phi)
%                     sin(psi)*sin(theta)                            cos(psi)*sin(theta)                             cos(theta) ]            
            wRp = rotz(phi*180/pi)*rotx(theta*180/pi)*rotz(psi*180/pi);
        end
        
        function [Jinv, J1inv, J2inv] = ijac(obj, q)
            % Calculate inverse jacobian 
            % Also returns the components Jinv1 and Jinv2 as have been used
            % in the article this code is been based on
            
            [~, n, a_w, ~, wRp] = obj.ikine(q); % Get inverse kinematics parameters
            
            J1inv=zeros(6,6);               % Calculate first jacobian
            for i = 1:6
                J1inv(i,:) = [n(:,i)' (wRp*cross(a_w(:,i), n(:,i)))'];
            end
            
            angles=q(4:6);                  % Calculate second jacobian
            phi = angles(1); theta = angles(2);
            tmp = [0 cos(phi) sin(phi)*sin(theta)
                   0 sin(phi) -cos(phi)*sin(theta)
                   1 0        cos(theta)];
            J2inv = [eye(3)     zeros(3)
                     zeros(3)   tmp];
            
            Jinv = J2inv*J1inv;             % Combine
        end
        
        function [q, qd] = acc_to_v_p(~, qdd, dt)
            % Numerical integration of an acceleration vector to velocity 
            % and position vectors
            % -- Input --
            % qdd - List of carthesian acceleration vector 
            % dt  - Time step between vectors
            % -- Output --
            % q   - List of position vectors
            % qd  - List of velocity vectors
            qd = dt*cumtrapz(qdd, 2);
            q  = dt*cumtrapz(qd,  2);
        end
        
        function [fa, j_vel] = idyn(obj, qdd, dt)
            G = [0 -9.81 0]';
            [q, qd] = obj.acc_to_v_p(qdd, dt);
            % Move around home position
            q = q + obj.home_pos;
            % Initialize force vector
            fa = zeros(6, length(qdd(1,:)));
            j_vel = zeros(6, length(qdd(1,:)));
            % For each timestep
            for ts = 1:length(qdd(1,:))
                % Inverse kinematics data
                [l, n, ~, ~, wRp] = obj.ikine(q(:,ts));
                % Inverse jacobians
                [Jinv, Jinv1, ~] = obj.ijac(q(:,ts));
                
                j_vel(:,ts) = Jinv*qd(:,ts);
                
                Ip = wRp*obj.Plat_I*wRp'; % Transform platform moment of inertia to world
                r = wRp*obj.Plat_COM; % Transform platform center of mass to world
                
                % Platform rotational velocity and acceleration
                % Get rotational position and rates
                phi=q(4,ts); phid=qd(4,ts); phidd=qdd(4,ts);
                theta=q(5,ts); thetad=qd(5,ts); thetadd=qdd(5,ts);
                               psid=qd(6,ts); psidd=qdd(6,ts);
                
                % Prepare trig values
                sphi = sin(phi); cphi=cos(phi);
                stheta = sin(theta); ctheta= cos(theta);
                
                % Rotational velocity in platform frame
                omega = [0 cphi sphi*ctheta
                         0 sphi -cphi*stheta 
                         1 0    ctheta]...
                      * [phid thetad psid]';
                
                % Rotational acceleration in platform frame
                alpha = [0 cphi sphi*stheta
                         0 sphi -cphi*stheta
                         1 0    ctheta]...
                         * [phidd thetadd psidd]'...
                      + [0 -phid*sphi phid*cphi*stheta+thetad*sphi*ctheta
                         0 phid*cphi  phid*sphi*stheta-thetad*cphi*ctheta
                         0 0          -thetad*stheta]...
                         * [phid thetad psid]';
                
                % Linear acceleration of platform COM (7, eq. 65)
                xgdd = qdd(1:3,ts) + cross(alpha,r) + cross(omega, cross(omega,r));
                
                % Calculate normal moment N for each joint
                Fn = zeros(3,6);
                % Store a1, will be used to calculate actuator force
                a1 = zeros(3,6);
                for i = 1:6
                    % Get some variables for specific joint
                    nj = n(:,i); lj = l(i); aj = wRp*obj.Pk(:,i);
                    xd = qd(1:3,ts); xdd=qdd(1:3,ts);
                    l1 = obj.Act_l1; l2 = obj.Act_l2;
                    Inn1 = obj.Act_Inn1; Inn2 = obj.Act_Inn2;
                    m1 = obj.Act_m1; m2 = obj.Act_m2;
                    
                    % Calculations
                    % Linear velocity of attachment point (4.2, eq. 10)
                    ajd = xd + cross(omega, aj);
                    % Linear acceleration (4.3, eq. 17)
                    ajdd = xdd + cross(alpha, wRp*aj) + cross(omega, cross(omega, wRp*aj));
                    
                    % Joint velocity (6, eq. 36)
                    ljd = dot(ajd, nj);
                    
                    % Joint rotational velocity (6, eq. 49)
                    omega_j = cross(nj, ajd)*lj; 
                    % Joint rotational acceleration (6, eq. 50)
                    alpha_j = (cross(nj, ajd) - 2*ljd*omega_j)/lj;
                    
                    % Joint acceleration (6, eq. 42)
                    ljdd = dot(ajdd, nj) - dot(lj*cross(omega_j, cross(omega_j,nj)), nj);
                    
                    % Linear acceleration of actuator parts (6, eq. 51 and 52)
                    aj1 = (lj-l1)*cross(omega_j, cross(omega_j,nj))...
                          + (lj-l1)*cross(alpha_j, nj)...
                          + 2*cross(omega_j, ljd*nj) + ljdd*nj;
                    aj2 = l2*cross(omega_j, cross(omega_j, nj)) + l2*cross(alpha_j, nj);
                    
                    % Moment balance (7, eq. 57)
                    % Uses simplification - assuming no rotation around
                    % actuation axis
                    Nj = -m1*(l(i)-l1)*cross(nj, G) -m2*l2*cross(nj, G)...
                             +(Inn1+Inn2)*cross(nj, cross(alpha_j, nj))...
                             -(-Inn1-Inn2)*(dot(omega_j, nj))*cross(nj, omega_j)...
                             +m1*(lj-l1)*cross(nj, aj1) + m2*l2*cross(nj, aj2);                             
                    
                    % Calculate normal force for actuator (7, eq. 61)
                    Fn(:,i) = cross(Nj,nj)/lj;
                    
                    % Store a1
                    a1(:,i) = aj1;
                end
                
                % C vector (7, eq. 69)
                Clin = obj.Plat_m*(G-xgdd) - sum(Fn,2);  
                Crot = obj.Plat_m*(cross(r,G) - cross(r,xgdd))...
                       - Ip*alpha+Ip*cross(omega,omega) - sum(wRp*cross(obj.Pk,Fn), 2);
                C = [Clin; Crot];
                
                % Axial force (7, eq. 68)
                Fax = inv(Jinv1)'*C;
                
                % Calculate actuator forces (7, eq. 72)
                F = zeros(6,1);
                for i = 1:6
                    m1 = obj.Act_m1; nj = n(:,i); aj1 = a1(:,i);
                    F(i) = dot(m1*(aj1-G), nj) - Fax(i);
                end
                
                % Store force for current timestep
                fa(:,ts) = F;
            end
        end

        function ikplot(obj)
            % Plot the platform in 3D with sliders to control position and
            % rotation
            
            % Parameters
            lin_amp = 0.300;  % Max range of motion
            ang_amp = 0.5;  % Max range of motion angular 
            
            % Make the figure window
            f = figure;
            ax = axes('Parent', f, 'position', [0.15 0.3 0.75 0.75]);
            grid on; axis equal;
            
            % Sliders
            x=0;y=0;z=0;a=0;b=0;c=0;           
            s_w = 220; s_h = 20; s_gap = 10; s_0 = [50,50];
            bx = uicontrol('Parent',f,'Style','slider','Position',[s_0(1)+0*(s_w+s_gap), s_0(2)+2*(s_h+s_gap), s_w, s_h],...
              'value',x, 'min',-lin_amp, 'max',lin_amp);
            by = uicontrol('Parent',f,'Style','slider','Position',[s_0(1)+0*(s_w+s_gap), s_0(2)+1*(s_h+s_gap), s_w, s_h],...
              'value',y, 'min',-lin_amp, 'max',lin_amp);
            bz = uicontrol('Parent',f,'Style','slider','Position',[s_0(1)+0*(s_w+s_gap), s_0(2)+0*(s_h+s_gap), s_w, s_h],...
              'value',z, 'min',-lin_amp, 'max',lin_amp);
            ba = uicontrol('Parent',f,'Style','slider','Position',[s_0(1)+1*(s_w+s_gap), s_0(2)+2*(s_h+s_gap), s_w, s_h],...
              'value',a, 'min',-ang_amp, 'max',ang_amp);
            bb = uicontrol('Parent',f,'Style','slider','Position',[s_0(1)+1*(s_w+s_gap), s_0(2)+1*(s_h+s_gap), s_w, s_h],...
              'value',b, 'min',-ang_amp, 'max',ang_amp);
            bc = uicontrol('Parent',f,'Style','slider','Position',[s_0(1)+1*(s_w+s_gap), s_0(2)+0*(s_h+s_gap), s_w, s_h],...
              'value',c, 'min',-ang_amp, 'max',ang_amp);
            addlistener(bx,'Value','PreSet',@(~,~)updatePlot()); addlistener(by,'Value','PreSet',@(~,~)updatePlot());
            addlistener(bz,'Value','PreSet',@(~,~)updatePlot()); addlistener(ba,'Value','PreSet',@(~,~)updatePlot());
            addlistener(bb,'Value','PreSet',@(~,~)updatePlot()); addlistener(bc,'Value','PreSet',@(~,~)updatePlot());
            
            function [rx, ry, rz] = getAngles(q)
              [l, n_vec, ~, ~, wRp]  = obj.ikine(q);
              n_r = wRp*n_vec(:,1); % Choose first, doesn't matter)
              disp("Angles:")
              disp(n_r)
              % Rotate to attachment point space - x axis through home pos.
              % normalized actuator axis, y axis towards platform center
              [~, n_home, ~, ~, wRp_home] = obj.ikine(obj.home_pos);
              n_home_x = wRp_home*n_home(:,1);
              % Project to XY-plane, rotate -90 in Z, normalize
              n_home_y = rotz(-90)*[n_home_x(1); n_home_x(2); 0];
              n_home_y = n_home_y/norm(n_home_y);
              n_home_z = cross(n_home_x, n_home_y); % Y axis is just the cross product
              pRa = [n_home_x'; n_home_y'; n_home_z']; % Rotation matrix to actuator attachment point coordinates
              n_r = pRa*n_r;
              disp(n_r)
              rx = atan2(n_r(3), n_r(2));
              ry = atan2(-n_r(3), n_r(1));
              rz = atan2(n_r(2), n_r(1));
              min_length = min(l);
              max_length = max(l);
              if min_length < obj.Act_fixed_l || max_length > obj.Act_max_l
                rx = 0; ry = 0; rz = 0;
              end
            end
            
            function updatePlot()
                % Get pose from sliders
                x = get(bx, "Value"); y = get(by, "Value"); z = get(bz, "Value");
                a = get(ba, "Value"); b = get(bb, "Value"); c = get(bc, "Value");
                q = [x y z a b c]' + obj.home_pos;

                % Get IK positions and rotation matrix
                [l, n, a_w, ~, wRp] = obj.ikine(q);
                % Print out act. 1 angles relative to platform
                angles_act1 = zeros(1, 3);
                [angles_act1(1), angles_act1(2), angles_act1(3)] = getAngles(q);
                angles_act1 = angles_act1 * 180/pi;
                disp(angles_act1)
                % Only visualize if position is legal!
                if sum(l>obj.Act_max_l) || sum(l<obj.Act_fixed_l)
                    return;
                end
                hold off; [caz, cel] = view(ax); % Store camera view
                
                % Base points
                plot3(ax, obj.Bk(1,:), obj.Bk(3,:), obj.Bk(2,:), "O")
                hold on; view(ax, caz, cel); axis equal; % Enable drawing over
                
                % Top points
                plot3(ax, a_w(1,:), a_w(3,:), a_w(2,:), "O");
                
                % Joints
                for i = 1:6
                    % Full length
                    X = [obj.Bk(1,i) a_w(1,i)];
                    Y = [obj.Bk(3,i) a_w(3,i)];
                    Z = [obj.Bk(2,i) a_w(2,i)];
                    plot3(ax, X, Y, Z, "-r", "LineWidth", 3);
                    
                    % Fixed length
                    % Plot from platform along leg vector for fixed
                    % distance
                    X = [obj.Bk(1,i) obj.Bk(1,i)+obj.Act_fixed_l*n(1,i)];
                    Y = [obj.Bk(3,i) obj.Bk(3,i)+obj.Act_fixed_l*n(3,i)];
                    Z = [obj.Bk(2,i) obj.Bk(2,i)+obj.Act_fixed_l*n(2,i)];
                    if i == 1
                        plot3(ax, X, Y, Z, "-r", "LineWidth", 6);
                    elseif i == 2
                        plot3(ax, X, Y, Z, "-g", "LineWidth", 6);
                    else
                        plot3(ax, X, Y, Z, "-k", "LineWidth", 6);
                    end
                end
                
                % Platform/base circles
                theta = linspace(0, 2*pi, 90)'; 
                p_circ = obj.Pr * [cos(theta) sin(theta) zeros(length(theta),1)]';
                p_circ = q(1:3) + wRp*p_circ;
                b_circ = obj.Br * [cos(theta) sin(theta) zeros(length(theta),1)]';
                b_circ = rotx(90)*b_circ;
               
                plot3(ax, p_circ(1,:), p_circ(3,:), p_circ(2,:), "-r", "LineWidth", 2);
                plot3(ax, b_circ(1,:), b_circ(3,:), b_circ(2,:), "-b", "LineWidth", 2)

                % Set plot limits
                xlim([-0.600-1.2*lin_amp +0.600+1.2*lin_amp]);
                ylim([-0.600-1.2*lin_amp +0.600+1.2*lin_amp]);
                zlim([0 obj.home_pos(2)+1.5*lin_amp]);
                
                delete(findall(gcf,'Tag','stream'));
                annotation(f, 'textbox', [0.01, 0.05, 0, 0], ...
                            'FitBoxToText','on', ...
                            'string', sprintf('X%05.3f  Y%05.3f  Z%05.3f  A%05.2f  B%05.2f  C%05.2f', q-obj.home_pos), ...
                            'Tag','stream');
            end
            % Initialize plot
            updatePlot();
            view(3);
        end
                
        function plot_motion(obj, q, dt)
            % Visualize the platform moving between positions
            % -- Input --
            % q  - series of position vectors 
            %      [x y z phi theta psi]' relative to the set platform home
            %      position. Linear position in meters, angular in rad
            % dt - time step in seconds
            % Parameters
            lin_amp = 0.200;  % Max range of motion 
            ang_amp = 0.5;  % Max range of motion angular
            % Make the figure window
            f = figure;
            ax = axes('Parent', f, 'position', [0.1 0.1 0.9 0.9]);
            grid on; axis equal;
            
            pause(1)
            for i = 1:length(q(1,:))
                qi = q(:,i) + obj.home_pos;
                
                % Get IK positions and rotation matrix
                [~, n, a_w, ~, wRp] = obj.ikine(qi);
                
                hold off;
                % Base points
                plot3(ax, obj.Bk(1,:), obj.Bk(3,:), obj.Bk(2,:), "O")
                hold on; % Enable drawing over
                
                % Top points
                plot3(ax, a_w(1,:), a_w(3,:), a_w(2,:), "O");
                
                % Joints
                for j = 1:6
                    % Full length
                    X = [obj.Bk(1,j) a_w(1,j)];
                    Y = [obj.Bk(3,j) a_w(3,j)];
                    Z = [obj.Bk(2,j) a_w(2,j)];
                    plot3(ax, X, Y, Z, "-r", "LineWidth", 3);
                    
                    % Fixed length
                    % Plot from platform along leg vector for fixed
                    % distance
                    X = [obj.Bk(1,j) obj.Bk(1,j)+obj.Act_fixed_l*n(1,j)];
                    Y = [obj.Bk(3,j) obj.Bk(3,j)+obj.Act_fixed_l*n(3,j)];
                    Z = [obj.Bk(2,j) obj.Bk(2,j)+obj.Act_fixed_l*n(2,j)];
                    plot3(ax, X, Y, Z, "-k", "LineWidth", 6);
                end
                
                % Platform/base circles
                theta = linspace(0, 2*pi, 90)'; 
                p_circ = obj.Pr * [cos(theta) sin(theta) zeros(length(theta),1)]';
                p_circ = qi(1:3) + wRp*p_circ;
                b_circ = obj.Br * [cos(theta) sin(theta) zeros(length(theta),1)]';
                b_circ = rotx(90)*b_circ;
               
                plot3(ax, p_circ(1,:), p_circ(3,:), p_circ(2,:), "-r", "LineWidth", 2);
                plot3(ax, b_circ(1,:), b_circ(3,:), b_circ(2,:), "-b", "LineWidth", 2)

                % Set plot limits
                xlim([-0.600-1.2*lin_amp +0.600+1.2*lin_amp]);
                ylim([-0.600-1.2*lin_amp +0.600+1.2*lin_amp]);
                zlim([0 obj.home_pos(2)+1.5*lin_amp]);
                % Plot from default view
                view(3);
                
                % Wait for next frame
                pause(dt)
            end
        end  
        
        function ROMplot(obj, q_0, amp, N)
            % 3D plot of how far axes are from end of motion
            % -- Input --
            % q_0 - Starting pose in absolute world coordinates
            % amp - amplitude of linear motions that will be plotted.
            %       e.g. 0.5 will plot a volume of 0.5x0.5x0.5m around
            %       q_0
            % N   - Number of samples to take per axis
            s = zeros(N, N, N);
            for n = 1:N
                for m = 1:N
                    for h = 1:N
                        x = amp*2*(-N/2+n)/N;
                        y = amp*2*(-N/2+h)/N;
                        z = amp*2*(-N/2+m)/N;
                        q = [x y z 0 0 0]' + q_0;
                        s(n,m,h) = goodness(q);
                    end
                end
            end
            figure()
            sliceViewer(s)
            
            function [goodness] = goodness(q_i)
                % Naive approach, a better solution could be using the inverse
                % jacobian
                l = obj.ikine(q_i);
                l_2 = (obj.Act_max_l - obj.Act_fixed_l) / 2; % half length
                % Get all distances from ends
                dists = [l-obj.Act_fixed_l; obj.Act_max_l-l];
                min_length = min(dists);
                goodness = min_length / l_2 ; % Normalize the worst value
                %goodness = goodness * (goodness>0); % Constrain to positive only
                goodness = (0.2 + goodness/2)*(goodness>0);
            end
        end
    
        function AngleROM(obj, q_0, amp, N)
            % 3D plot of how far axes are from end of motion
            % -- Input --
            % q_0 - Starting pose in absolute world coordinates
            % amp - amplitude of linear motions that will be plotted.
            %       e.g. 0.5 will plot a volume of 0.5x0.5x0.5m around
            %       q_0
            % N   - Number of samples to take per axis
            angles_x = zeros(1, N^3);
            angles_y = zeros(1, N^3);
            angles_z = zeros(1, N^3);
            for n = 1:N
                for m = 1:N
                    for h = 1:N
                        x = amp*2*(-N/2+n)/N;
                        y = amp*2*(-N/2+h)/N;
                        z = amp*2*(-N/2+m)/N;
                        q = [x y z 0 0 0]' + q_0;
                        idx = n + m*N + h*N^2;
                        [angles_x(idx), angles_y(idx), angles_z(idx)] = getAngles(q);
                    end
                end
            end
            function [rx, ry, rz] = getAngles(q)
              [l, n_vec, ~, ~, wRp]  = obj.ikine(q);
              n_r = wRp*n_vec(:,1); % Choose first, doesn't matter)
              
              % Rotate to attachment point space - x axis through home pos.
              % normalized actuator axis, y axis towards platform center
              [~, n_home, ~, ~, wRp_home] = obj.ikine(q_0);
              n_home_x = wRp_home*n_home(:,1);
              % Project to XY-plane, rotate -90 in Z, normalize
              n_home_y = rotz(-90)*[n_home_x(1); n_home_x(2); 0];
              n_home_y = n_home_y/norm(n_home_y);
              n_home_z = cross(n_home_x, n_home_y); % Z axis is just the cross product
              pRa = [n_home_x'; n_home_y'; n_home_z']; % Rotation matrix to actuator attachment point coordinates
              n_r = pRa*n_r;
              
              rx = 180/pi * atan2(n_r(3), n_r(2));
              ry = 180/pi * atan2(-n_r(3), n_r(1));
              rz = 180/pi * atan2(n_r(2), n_r(1));
              min_length = min(l);
              max_length = max(l);
              if min_length < obj.Act_fixed_l || max_length > obj.Act_max_l
                rx = 0; ry = 0; rz = 0;
              end
            end
            angles_x = sort(nonzeros(angles_x)');
            %angles_x = angles_x - mean(angles_x);
            angles_y = sort(nonzeros(angles_y)');
            %angles_y = angles_y - mean(angles_y);
            angles_z = sort(nonzeros(angles_z)');
            %angles_z = angles_z - mean(angles_z);
            figure();
            subplot(3, 1, 1);
            histogram(angles_x);
            xlabel("X angles");
            subplot(3, 1, 2);
            histogram(angles_y);
            xlabel("Y angles");
            subplot(3, 1, 3);
            histogram(angles_z);
            xlabel("Z angles");
        end
    end
end