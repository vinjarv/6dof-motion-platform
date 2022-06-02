% Define platform
s = stplat;

% Time step and time vector
ts = 0.005;
t = 0:ts:2.4;
len = length(t);

% Initialize acceleration vector
qdd = zeros(6,len);

% Motions to run
a = 2*9.81;
ra_deg = 250; %deg/s^2
ra = ra_deg * pi/180;
move_time = 0.2; % Time per accel move

% Create accelerations
qdd = motion1w(qdd, 0*move_time, 1*move_time, t, [a 0 0 0 0 0]');
qdd = motion1w(qdd, 1*move_time, 2*move_time, t, [-a 0 0 0 0 0]');

qdd = motion1w(qdd, 2*move_time, 3*move_time, t, [0 a 0 0 0 0]');
qdd = motion1w(qdd, 3*move_time, 4*move_time, t, [0 -a 0 0 0 0]');

qdd = motion1w(qdd, 4*move_time, 5*move_time, t, [0 0 a 0 0 0]');
qdd = motion1w(qdd, 5*move_time, 6*move_time, t, [0 0 -a 0 0 0]');

qdd = motion1w(qdd, 6*move_time, 7*move_time, t, [0 0 0 ra 0 0]');
qdd = motion1w(qdd, 7*move_time, 8*move_time, t, [0 0 0 -ra 0 0]');

qdd = motion1w(qdd, 9*move_time, 9*move_time, t, [0 0 0 0 ra 0]');
qdd = motion1w(qdd, 9*move_time, 10*move_time, t, [0 0 0 0 -ra 0]');

qdd = motion1w(qdd, 10*move_time, 11*move_time, t, [0 0 0 0 0 ra]');
qdd = motion1w(qdd, 11*move_time, 12*move_time, t, [0 0 0 0 0 -ra]');

% Simulate forces and velocities
[f, v] = s.idyn(qdd, ts);


figure(1)
plot(t, f)
grid()
legend("Actuator 1", "Actuator 2", "Actuator 3", "Actuator 4", "Actuator 5", "Actuator 6");
title("Actuator force")
xlabel("Time (s)")
ylabel("Force (N)")
ylim([-2000 2000]);

figure(2)
%ballscrew_lead = 10e-3;
%rpm = v/ballscrew_lead * 60
plot(t, v)
grid()
legend("Actuator 1", "Actuator 2", "Actuator 3", "Actuator 4", "Actuator 5", "Actuator 6");
title("Actuator velocity")
xlabel("Time (s)")
ylabel("Velocity (m/s)")

% Adds constant acceleration between t0 and t1
function [qdd] = motion1w(qdd, t0, t1, t, qdd0)
    dt = t1-t0;
    m1 = (t>t0).*(t<(t0+(t1-t0)/2));
    m2 = (t>(t0+(t1-t0)/2)).*(t<t1);
    for i = 1:length(qdd(1,:))
        qdd(:, i) = qdd(:, i) + qdd0*m1(i);
        qdd(:, i) = qdd(:, i) - qdd0*m2(i);
    end
end
