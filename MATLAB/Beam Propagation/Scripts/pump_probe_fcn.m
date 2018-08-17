function Trans = pump_probe_fcn(alpha_2_nd_cm_per_GW, pump_pulse_width_HW_einv_max_fs, delay_offset_fs, probe_pulse_width_HW_einv_max_fs, pump_energy_J, sample_thickness_m, pump_spot_size_x_HW_einv2_max_m, pump_spot_size_y_HW_einv2_max_m, probe_spot_size_HW_einv2_max_m, t_d)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Sample Parameters

alpha_2 = alpha_2_nd_cm_per_GW*10 * 1e-12;                         % ND-2PA Coefficient [m/W]
alpha_3 = 0.0 * 1e-24;                        % ND-3PA Coefficient [m^3/W^2]

n_2 = 0 * 1e-19;                             % ND third order NLR index [m^2/W]
n_4 = 0.00 * 1e-36;                             % ND fifth order NLR index [I don't know]

Lin_Trans = 1;                                  % Linear Transmission of probe 

n_pump = 3.4258;                                   % Index at Pump            refractiveindex.info
ng_pump = 3.4441;                          % Dispersion at Pump [m^-1]
n_probe = 3.522;                                     % Index at Probe
ng_probe = 3.7620;


Conc = 0.0091;                                  	% Concentration [M = mol/L]
% 
% % Beam Parameters
% 
d = 0.25;                                       % Distance from Sample to Detector [m]
% 
lambda_pump = 4000e-9;                             % Pump wavelength [m]
lambda_probe = 1200e-9;                                % Probe wavelength [m]
tau_pump = pump_pulse_width_HW_einv_max_fs * 1e-15;          % Pump temporal pulse width (HW1/e of I) [s]
tau_probe = probe_pulse_width_HW_einv_max_fs * 1e-15;            % Probe temporal pulse width (HW1/e of I) [s]

X0 = 0;                                       % Pump-Probe x-axis displacement [w_p]

scale = 1;

% Simulation Parameters

X_max = 10;                                    	% Maximum of Normalized Spatial Vector in X-direction [w]
X_num = 50;                                  	% Number of points in T and T_d

T_max = 20;                                     % Maximum of Normalized Temporal Vector [tau_p] (should probably write this in terms of rho
T_num = 500;                                     % Number of points in T and T_d

%%
% CONSTANTS
%

c = 299792458;                      % Speed of light [m/s]
epsilon_0 = 8.854187817620e-12;     % Permittivity of free space [F/m]
h = 6.62606957e-34;                 % Plank's constant [J s]
N_A = 6.02214129e23;                % Avagadro's constant [mol^-1]

%%
% CALCULATED PARAMETERS
%
E_ph_pump = h*c/lambda_pump;                      % Pump Photon energy [J]
k_0_pump = 2*pi/lambda_pump;                      % Pump Wavenumber [m^-1]
k_0_probe = 2*pi/lambda_probe;                          % Probe Wavenumber [m^-1]

N = Conc * N_A * 10^3;                      % Concentration [m^-3]
delta = E_ph_pump*alpha_2/N;                   % 2PA cross-section [m^4 s]
delta_r = E_ph_pump*n_2*k_0_pump/N;               % NLR cross-section [m^4 s]

alpha = -log(Lin_Trans)/sample_thickness_m;                          % Linear Absorption Coefficient of Pump [m^-1]

sigma_probe = alpha * sample_thickness_m/2;                                % Linear Absorption parameter

R_pump = ( (n_pump-1) / (n_pump+1) )^2;                      % Refectivity of Pump at first interface
% R_p = 0;
I_0_p = 2*pump_energy_J/(pi^(3/2)*pump_spot_size_x_HW_einv2_max_m*pump_spot_size_y_HW_einv2_max_m*tau_pump) * (1-R_pump);   % Peak Irradiance of Pump [W/m^2]

t = tau_probe/tau_pump;                                      % Normalized Probe Pulse Width
W_0_pump_x = pump_spot_size_x_HW_einv2_max_m/probe_spot_size_HW_einv2_max_m;                                  % Normalized Pump Radius
W_0_pump_y = pump_spot_size_y_HW_einv2_max_m/probe_spot_size_HW_einv2_max_m;
D = d/probe_spot_size_HW_einv2_max_m;
LAMBDA = lambda_probe/probe_spot_size_HW_einv2_max_m;
K_0 = 2*pi/LAMBDA;
z_0 = pi*probe_spot_size_HW_einv2_max_m^2/lambda_probe;

rho = sample_thickness_m/(tau_pump*c)*(ng_probe-ng_pump);                        % GVM Parameter

eta = 4*pi * sample_thickness_m/lambda_probe * I_0_p*n_2;
Gamma = sample_thickness_m * I_0_p*alpha_2;
chi = eta + 1i*Gamma;                               % Third-Order Nonlinear Parameter

dX = 2*X_max/X_num;                                 % Differential Unit of Normalized Spatial Vector in X-direction [w]
dT =  2*(T_max+abs(rho))/T_num;                     % Differential Unit of Normalized Temporal Vector [tau_p]

%%
% Initializeation
%

X = -X_max:dX:X_max;                                % Normalized Spatial Vector in X-direction [w]
Y = -0:dX:X_max;                                    % Normalized Spatial Vector in Y-direction [w]

if rho == 0
    T_L = -T_max;
    T_U = T_max;
elseif rho > 0
    T_L = -T_max-abs(rho);
    T_U = T_max;
else
    T_L = -T_max;
    T_U = T_max+abs(rho);
end
T = T_L:dT:T_U;                                     % Normalized Temporal Vector [tau_p]
T2 = T;
T_d = t_d*1e-15;
T_d = T_d/tau_pump;                                 % Normalized Temporal Delay Vector [tau_p]
[Y,X,T] = ndgrid(Y,X,T);                    % Normalized Spatio-Temporal Matrix [w, tau_p] (y,x,t,t_d)
[T2, T_d2] = ndgrid(T2, T_d);

a_p_0 =  exp(-((X+X0*W_0_pump_x).^2/W_0_pump_x^2 + Y.^2/W_0_pump_y^2));
%%
% CALCULATIONS
%

if rho ~= 0
    a_out = exp(-X.^2-Y.^2) .* exp(- sigma_probe ...
                           +1i*chi*sqrt(pi)/(2*rho).*a_p_0.^2.*(erf(T)-erf(T-rho)) );
else                    
    a_out = exp(-X.^2-Y.^2) .* exp(-1/2.*(T-T_d).^2 / t^2 - sigma_probe ...
                          +1i*chi.*a_p_0.^2.*exp(-T.^2) ...
                          +1i*psi.*a_p_0.^4.*exp(-2*T.^2));
end

a_det = fftshift(fftshift(fft(fft(a_out.*exp(1i*K_0/(2*D)*((X).^2+(Y).^2)),[],2),[],1),2),1);

E_tot = squeeze((trapz(trapz(abs(a_det).^2,1),2)));
T_d_term = abs(exp(-1/2*(T2 - T_d2 - rho - delay_offset_fs*1e-15/tau_pump).^2/t^2)).^2;

E_tot = bsxfun(@times, E_tot, T_d_term);
E_tot = squeeze(trapz(E_tot,1));

E_left = squeeze(trapz(trapz(abs(a_det(:,1:ceil(end/2),:)).^2,1),2));
E_left = bsxfun(@times, E_left, T_d_term);
E_left = squeeze(trapz(E_left,1));

E_right = squeeze(trapz(trapz(abs(a_det(:,ceil(end/2):end,:)).^2,1),2));
E_right = bsxfun(@times, E_right, T_d_term);
E_right = squeeze(trapz(E_right,1));

Trans = transpose(E_tot / max(E_tot));
Delta_E = E_left - E_right;
S = Delta_E ./ E_tot;

t_d = (T_d)*tau_pump*1e15;

end

