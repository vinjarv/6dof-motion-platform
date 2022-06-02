%% Manual HP filter
zeta = 0.7 ;
wn = 1.5; % rad/s
K = 1/20;

fs = 100; % Hz

G_hp = tf([K 0 0], [1 2*zeta*wn, wn^2])
G_int = tf(1, [1 0 0])
G_trans = G_hp*G_int

G_z_trans = c2d(G_trans, 1/fs)
b = G_z_trans.Numerator{1}
a = G_z_trans.Denominator{1}

%freqz(b, a)
step(18*G_z_trans, 5)