% Create symbolic form of rotation matrix for implementation in PLC
% Verify results by comparing against MATLAB functions
%% ZXZ Rotation matrix
clear; clc;
syms a b c
rz_a = [cos(a) -sin(a) 0;
        sin(a) cos(a)  0;
        0      0       1];
rx_b = [1   0      0;
        0   cos(b) -sin(b);
        0   sin(b) cos(b)];
rz_c = [cos(c) -sin(c) 0;
        sin(c) cos(c)  0;
        0      0       1];
disp("ZXZ Rotation matrix:")
pretty(rz_a*rx_b*rz_c)

%% Compare
clear; clc;
% Inputs:
phi = 2*pi/3;
theta = 2*pi/4;
psi = 5*pi/5;
testvec = [1 2 5]';

% Rotation matrices
wRp = rotz(phi*180/pi)*rotx(theta*180/pi)*rotz(psi*180/pi);
    
disp("MATLAB Rotation matrix:")
disp(wRp*testvec)

disp("Manual rotation:")
a=phi;b=theta;c=psi;
sa=sin(a);ca=cos(a);
sb=sin(b);cb=cos(b);
sc=sin(c);cc=cos(c);
disp([testvec(1)*(ca*cc-cb*sa*sc) + testvec(2)*(-ca*sc-cb*cc*sa) + testvec(3)*(sa*sb);
      testvec(1)*(cc*sa+ca*cb*sc) + testvec(2)*(ca*cb*cc-sa*sc) + testvec(3)*(-ca*sb);
      testvec(1)*(sb*sc) + testvec(2)*(cc*sb) + testvec(3)*(cb)])