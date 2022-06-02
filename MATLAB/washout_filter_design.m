% Calculation of filter coefficients for 6DOF platform washout

% 3 HP filters + double integration for linear accelerations
% 3 LP filters for linear acceleration with scaling to generate fixed angle
% offsets
% 3 HP filters with integration for rotation rates

clear; clc
format shortG
%% Translational channels
% Highpass filter and double integration

fs = 250;% Hz
K = 1.8;
zeta = 0.95;
wn = 11;

% 1st try - z 0.75 wn 3.6 reg_k 0.18 not very responsive
% 2nd try - z 1.5 wn 15 reg_k 0.4 very bumpy (around 8cm at 2g step, seems good)
% 3rd try - z 1.1 wn 8 reg_k 0.4 not very responsive
% 4th try - z 0.95 wn 11 reg_k 0.5 pretty good

% Transfer function for filter
%g_hp = tf(num_hp, den_hp);
s = tf('s');
num = [K 0 0];
den = [1 2*zeta*wn wn^2];
% p3 = poly(10*min(real(roots(den))));
% den = conv(den, p3)
% den = conv(den, p3)
g_hp = tf(num, den)

figure(1)
step(g_hp)
%% Add feedback to return to zero
k_feedback = 0.5;

int1 = tf([1],[1 0]);
int2 = feedback(int1, k_feedback);

int = int1 * int2;

g_tot = g_hp * int
figure(2)
step(2*9.81*g_tot)

%% Discretize 
g_disc = c2d(g_tot, 1/fs, 'tustin')
b = g_disc.Numerator{1}
a = g_disc.Denominator{1}

fprintf("B:\n")
fprintf("%.20f,", b)
fprintf("\nA:\n")
fprintf("%.20f,", a)
figure(3)
freqz(b, a)

figure(4)
step(g_disc)

%% Load recorded data and run filter
load gamedata
sim("sim_with_dataset")