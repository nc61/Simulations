% Problem 4

c = 3e8;
eps0 = 8.85e-12;
lambda = 1.500*1e-6; %m
k0 = 2*pi/lambda;
D = 3*1e-12*1e-3*1e9; %s/(m m)

P_peak = 0.5; % W
n_fiber = 1.44;
A_eff = 60*(1e-6)^2; %m^2
eta_0 = 1/(c*eps0);
n2_I = 3e-20;
n2 = n2_I*n_fiber/(2*eta_0); %m^2/V^2

beta_0 = -(lambda)^2/(2*pi*c)*D;
phi_peak_squared = 2*P_peak/(n_fiber*c*eps0*A_eff)
Omega_max = sqrt(k0*n2*phi_peak_squared/(-1*beta_0))
g = 1/2*k0*n2*phi_peak_squared


% Code for estimating power needed
L_eff = 2.17e3;
P_peak_min = 6;
k0*n2/2*2*P_peak_min/(n_fiber*c*eps0*A_eff)*L_eff