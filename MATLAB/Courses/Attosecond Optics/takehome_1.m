c = 3e8;
C = 9.33e-14;
q = 1.6e-19;
m_e = 9.109e-31;
eps_0 = 8.85e-12;

%% Functions

Up_eV = @(I_L_W_per_cm2, lambda_um)C*I_L_W_per_cm2*lambda_um^2;
x_function = @(t, x_0, w_0, w_0t_p)x_0/2*(cos(w_0*t) - cos(w_0t_p) + sin(w_0t_p)*(w_0*t - w_0t_p));
Ep = @(K_theta, irradiance_W_per_cm2, lambda_um, ionization_eV)K_theta*Up_eV(irradiance_W_per_cm2, lambda_um) + ionization_eV;
C_theta = @(theta)-1./(4*pi*K(theta).*(cos(theta) + pi/6*cos(1/3*theta - pi/6).*sin(pi/2*sin(theta/3 - pi/6))));
C_as_per_eV = @(dK_dtheta, Up_eV, T)1./dK_dtheta/(2*pi)*T*1e18/Up_eV;
as_per_eV_to_as2 = @(C_as_per_eV)C_as_per_eV/1.516*1e3;


%% Question 1

E_cutoff_eV = 560;
ionization_eV = 24.5874;
lambda_um = 1.7;
w_0 = 2*pi*c/(lambda_um*1e-6);

irradiance_W_per_cm2 = (E_cutoff_eV - ionization_eV)/(3.17*C*lambda_um^2)
E_field_V_per_m = sqrt(2*irradiance_W_per_cm2*100^2/(eps_0*c))
%% Question 2

T = (lambda_um*1e-6)/c;
t = linspace(-T, 3*T, 301);
E_field = E_field_V_per_m*cos(2*pi*t/T).*(heaviside(t) - heaviside(t-2*T));
figure(1)
plot(t*1e15, E_field);
xlabel('t [fs]'), ylabel('E [V/m]'), title('IR electric field')
xlim(1e15*[min(t) max(t)])

%% Question 3

x_0 = 2*q*E_field_V_per_m/(m_e*w_0^2);
t2 = linspace(0, 2*T, 301);
x = x_function(t2, x_0, w_0, 0.05*2*pi);
figure(2)
plot(t2*1e15,x*1e9);
xlabel('t [fs]'), ylabel('x [nm]'), title('Position of electron realeased at t_p = 0.05T');

%% Question 4

% Calculate the kinetic energy as a function of return time (in cycles).
% This way we can see the chirp in short and long trajectories
K = @(w_0t)2*(sin(w_0t) - cos(pi/2*sin(1/3*w_0t - pi/6))).^2;
theta = linspace(0.3*2*pi,2*pi,1000);
K_theta = K(theta);
figure(3)
plot(theta/(2*pi), Ep(K_theta, irradiance_W_per_cm2, lambda_um, ionization_eV))
xlabel('Return time (cycles)'), ylabel('E_p (eV)'), title('Photon energy vs return time');

C_as_per_eV_short = 1e18*0.069*T/Up_eV(irradiance_W_per_cm2, lambda_um);
C_as2_short = as_per_eV_to_as2(C_as_per_eV_short)
C_as_per_eV_long = 1e18*(-0.059)*T/Up_eV(irradiance_W_per_cm2, lambda_um);
C_as2_long = as_per_eV_to_as2(C_as_per_eV_long)


%% Question 5

% Plot the derivitave of K vs theta
figure(4)
dK_dtheta = diff(K_theta)/(theta(2) - theta(1));
plot(theta(2:end)/(2*pi), dK_dtheta)
xlabel('Return time [cycles]'), ylabel('dK/d\theta'), title('Derivative of photon energy wrt return time (phase angle)')

% Find the short trajectory inflection point
[max_dK_dtheta, t_max_index] = max(dK_dtheta);
short_trajectory_chirp = 1/max_dK_dtheta/(2*pi)

% Find the long trajectory inflection point
[min_dK_dtheta, t_min_index] = min(dK_dtheta);
long_trajectory_chirp = 1/min_dK_dtheta/(2*pi)

% Separate theta into short and long trajectories so that K(t) is single
% valued
[~, peak_dK_dtheta_index] = min(abs(dK_dtheta));
theta_short = linspace(0.3*2*pi, theta(peak_dK_dtheta_index), 200);
theta_long = linspace(theta(peak_dK_dtheta_index), 2*pi, 200);

% compute the short trajectory chirp
K_short = K(theta_short);
dK_dtheta_short = diff(K_short)/(theta_short(2) - theta_short(1));
C_short_as_per_eV = C_as_per_eV(dK_dtheta_short, Up_eV(irradiance_W_per_cm2, lambda_um), T);
C_short_as2 = as_per_eV_to_as2(C_short_as_per_eV);

% Compute the long trajectory chirp
K_long = K(theta_long);
dK_dtheta_long = diff(K_long)/(theta_long(2) - theta_long(1));
C_long_as_per_eV = C_as_per_eV(dK_dtheta_long, Up_eV(irradiance_W_per_cm2, lambda_um), T);
C_long_as2 = as_per_eV_to_as2(C_long_as_per_eV);

% Plot short and long trajectory chirps together
figure(5)
plot(Ep(K_long(2:end), irradiance_W_per_cm2, lambda_um, ionization_eV), C_long_as2)
hold on
plot(Ep(K_short(2:end), irradiance_W_per_cm2, lambda_um, ionization_eV), C_short_as2, 'r')
hold off
ylim(1e4*[-1 1]), xlabel('E_p [eV]'), ylabel('C [as^2]'), legend('long', 'short'), title('Long and short trajectory chirp vs photon energy')


