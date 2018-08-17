
alpha_2_nd_m_per_W = 10 * 1e-12;                         % ND-2PA Coefficient [m/W]

linear_transmission = 1;                                  % Linear Transmission of probe 
sample_thickness = 0.28e-3;                                    % Thickness [m]

pump_refractive_index = 3.4229;                                   % Index at Pump            refractiveindex.info
pump_group_index = 3.444; 
     % Dispersion at Pump [m^-1]
probe_refractive_index = 3.522;                                     % Index at Probe
probe_group_index = 3.7604;
   


% 
% % Beam Parameters
% 
% 
pump_wavelength = 4000e-9;                             % Pump wavelength [m]
probe_wavelength = 1250e-9;                                % Probe wavelength [m]
pump_pulsewidth = 105 * 1e-15;          % Pump temporal pulse width (HW1/e of I) [s]
probe_pulsewidth = 137 * 1e-15;            % Probe temporal pulse width (HW1/e of I) [s]

pump_energy = 4e-6;                                   % Pump input energy [J]

pump_spot_size = sqrt(2)*350e-6;                                	% Pump spot size (HW1/e^2 of I) [m]
probe_spot_size = sqrt(2)*90e-6;                                  	% Probe spot size (HW1/e^2 or I) [m]


% Simulation Parameters

X_max = 10;                                    	% Maximum of Normalized Spatial Vector in X-direction [w]
X_num = 50;                                  	% Number of points in T and T_d

T_max = 45;                                     % Maximum of Normalized Temporal Vector [tau_p] (should probably write this in terms of rho
T_num = 400;                                     % Number of points in T and T_d

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

alpha = -log(linear_transmission)/sample_thickness;                          % Linear Absorption Coefficient of Pump [m^-1]

sigma_probe = alpha * sample_thickness/2;                                % Linear Absorption parameter

R_pump = ( (pump_refractive_index-1) / (pump_refractive_index+1) )^2;                      % Refectivity of Pump at first interface
% R_p = 0;
pump_peak_irradiance = 2*pump_energy/(pi^(3/2)*pump_spot_size^2*pump_pulsewidth) * (1-R_pump)   % Peak Irradiance of Pump [W/m^2]

t = probe_pulsewidth/pump_pulsewidth;                                      % Normalized Probe Pulse Width
W_0_pump = pump_spot_size/probe_spot_size;                                  % Normalized Pump Radius
D = d/probe_spot_size;
LAMBDA = probe_wavelength/probe_spot_size;
K_0 = 2*pi/LAMBDA;
z_0 = pi*probe_spot_size^2/probe_wavelength;

rho = sample_thickness/(pump_pulsewidth*c)*(probe_group_index-pump_group_index);                        % GVM Parameter

Gamma = sample_thickness * pump_peak_irradiance*alpha_2_nd_m_per_W;
chi = 1i*Gamma;                               % Third-Order Nonlinear Parameter

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
T_d = 1e-15*linspace(-2500, 1500, 400);
T_d = T_d/pump_pulsewidth;                                 % Normalized Temporal Delay Vector [tau_p]
[Y,X,T] = ndgrid(Y,X,T);                    % Normalized Spatio-Temporal Matrix [w, tau_p] (y,x,t,t_d)
[T2, T_d2] = ndgrid(T2, T_d);

a_p_0 =  exp(-((X).^2 + Y.^2) / W_0_pump^2 );
%%
% CALCULATIONS
%

if rho ~= 0
    a_out = exp(-X.^2-Y.^2) .* exp(-sigma_probe ...
                           +1i*chi*sqrt(pi)/(2*rho).*a_p_0.^2.*(erf(T)-erf(T-rho)) );
else
    a_out = exp(-X.^2-Y.^2) .* exp( sigma_probe ...
                          +1i*chi.*a_p_0.^2.*exp(-T.^2) ...
                          +1i*psi.*a_p_0.^4.*exp(-2*T.^2));
end

a_det = fftshift(fftshift(fft(fft(a_out.*exp(1i*K_0/(2*D)*((X).^2+(Y).^2)),[],2),[],1),2),1);

E_tot = squeeze((trapz(trapz(abs(a_det).^2,1),2)));
T_d_offset = 0;
T_d_term = abs(exp(-1/2*(T2 - T_d2 - rho - T_d_offset/pump_pulsewidth).^2/t^2)).^2;

E_tot = bsxfun(@times, E_tot, T_d_term);
E_tot = squeeze(trapz(E_tot,1));

E_left = squeeze(trapz(trapz(abs(a_det(:,1:ceil(end/2),:)).^2,1),2));
E_left = bsxfun(@times, E_left, T_d_term);
E_left = squeeze(trapz(E_left,1));

E_right = squeeze(trapz(trapz(abs(a_det(:,ceil(end/2):end,:)).^2,1),2));
E_right = bsxfun(@times, E_right, T_d_term);
E_right = squeeze(trapz(E_right,1));

Trans = E_tot / max(E_tot);
Delta_E = E_left - E_right;
S = Delta_E ./ E_tot;

t_d = (T_d')*pump_pulsewidth*1e15;

plot(t_d, Trans)