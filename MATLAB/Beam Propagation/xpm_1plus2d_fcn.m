function [normalized_probe_energy_out, varargout] = xpm_1plus2d_fcn(x_mm, y_mm, t_fs, phi_pump_sqrt_GW_per_sqrt_cm2, ... 
    phi_probe_sqrt_GW_per_sqrt_cm2, sample_thickness_mm, num_z_steps, alpha_2_nd_cm_per_GW, alpha_2_d_cm_per_GW, ... 
    pump_refractive_index, probe_refractive_index, pump_wavelength_um, probe_wavelength_um, ... 
    pump_group_index, probe_group_index, ...
    pump_dispersion_fs2_per_mm, probe_dispersion_fs2_per_mm, output_normalized_energy_only)

c = 2.998e8;

xmax_mm = max(x_mm);
ymax_mm = max(y_mm);
tmax_fs = max(t_fs);

num_time_points = length(t_fs);
num_x_points = length(x_mm);
num_y_points = length(y_mm);

% Fourier space for x and y
kx_invmm = 2*pi/xmax_mm*[(0:num_x_points/2-1) (-num_x_points/2:-1)];
ky_invmm = 2*pi/ymax_mm*[(0:num_y_points/2-1) (-num_y_points/2:-1)];
omega_rad_per_fs = 2*pi/tmax_fs*[(0:num_time_points/2-1) (-num_time_points/2:-1)];

[Kx_invmm, Ky_invmm, OMEGA_rad_per_fs] = ndgrid(kx_invmm,ky_invmm,omega_rad_per_fs);

% nonlinear
alpha_2_nd_mm_per_GW = 10*alpha_2_nd_cm_per_GW;
alpha_2_d_mm_per_GW = 10*alpha_2_d_cm_per_GW;

phi_probe_sqrt_GW_per_sqrt_mm2 = (1e-1)*phi_probe_sqrt_GW_per_sqrt_cm2;
phi_pump_sqrt_GW_per_sqrt_mm2 = (1e-1)*phi_pump_sqrt_GW_per_sqrt_cm2;

% Propagation distance
dz_mm = sample_thickness_mm/num_z_steps;
dx_mm = xmax_mm/num_x_points;
dy_mm = xmax_mm/num_x_points;
dt_s = 1e-15*tmax_fs/num_time_points;

probe_x_integrated = squeeze(trapz(x_mm, 1e9*abs(phi_probe_sqrt_GW_per_sqrt_mm2).^2, 1));
probe_power = squeeze(trapz(y_mm, probe_x_integrated, 1));
initial_probe_energy_J = squeeze(trapz(1e-15*t_fs, probe_power))

pump_x_integrated = squeeze(trapz(x_mm, 1e9*abs(phi_pump_sqrt_GW_per_sqrt_mm2).^2, 1));
pump_power = squeeze(trapz(y_mm, pump_x_integrated, 1));
initial_pump_energy_J = squeeze(trapz(1e-15*t_fs, pump_power))


% Calculation of propagation constants
k_pump_per_mm = 1e3*pump_refractive_index*2*pi/pump_wavelength_um;
k_probe_per_mm = 1e3*probe_refractive_index*2*pi/probe_wavelength_um;

% Calculate group velocity
pump_group_velocity_mm_per_fs = 1e-12*c/pump_group_index;
probe_group_velocity_mm_per_fs = 1e-12*c/probe_group_index;

group_velocity_dispersion_fs_per_mm = 1/pump_group_velocity_mm_per_fs - 1/probe_group_velocity_mm_per_fs;

% Linear propagation operator
D_pump = -(0.5*(-1i*1/(2*k_pump_per_mm)*(Kx_invmm.^2 + Ky_invmm.^2) + 1i*pump_dispersion_fs2_per_mm*OMEGA_rad_per_fs.^2 - 1i*(group_velocity_dispersion_fs_per_mm)*OMEGA_rad_per_fs)*0.5*dz_mm);
D_probe = -(0.5*(-1i*1/(2*k_probe_per_mm)*(Kx_invmm.^2 + Ky_invmm.^2) + 1i*probe_dispersion_fs2_per_mm*OMEGA_rad_per_fs.^2)*0.5*dz_mm);


if (~output_normalized_energy_only)
    phi_out_probe = zeros(num_z_steps, num_x_points, num_y_points, num_time_points);
    phi_out_pump = zeros(num_z_steps, num_x_points, num_y_points, num_time_points);
end

z_mm = dz_mm:dz_mm:sample_thickness_mm;

% Not symmetrized, will need to do that for better accuracy
for ind = 1:num_z_steps
    
    %Pump
    lin_step_pump_1 = ifftn(fftn(phi_pump_sqrt_GW_per_sqrt_mm2).*exp(D_pump));
    nonlin_step_pump = lin_step_pump_1.*exp((-2*alpha_2_d_mm_per_GW*abs(lin_step_pump_1).^2)*dz_mm);
    lin_step_pump_2 = ifftn(fftn(nonlin_step_pump).*exp(D_pump));
    phi_pump_sqrt_GW_per_sqrt_mm2 = lin_step_pump_2;
    
    % Probe
    lin_step_probe_1 = ifftn(fftn(phi_probe_sqrt_GW_per_sqrt_mm2).*exp(D_probe));
    nonlin_step_probe = lin_step_probe_1.*exp((-2*alpha_2_nd_mm_per_GW*abs(lin_step_pump_1).^2)*dz_mm);
    lin_step_probe_2 = ifftn(fftn(nonlin_step_probe).*exp(D_probe));
    phi_probe_sqrt_GW_per_sqrt_mm2 = lin_step_probe_2;
    
    if (~output_normalized_energy_only)
        phi_out_pump(ind, :, :, :) = phi_pump_sqrt_GW_per_sqrt_mm2;
        phi_out_probe(ind,:,:, :) = phi_probe_sqrt_GW_per_sqrt_mm2;
    end
end

probe_x_integrated = squeeze(trapz(x_mm, 1e9*abs(phi_probe_sqrt_GW_per_sqrt_mm2).^2, 1));
probe_power = squeeze(trapz(y_mm, probe_x_integrated, 1));
final_probe_energy_J = squeeze(trapz(1e-15*t_fs, probe_power));

normalized_probe_energy_out = (initial_probe_energy_J - final_probe_energy_J)/initial_probe_energy_J;

if ~output_normalized_energy_only
    varargout{1} = phi_out_probe;
    varargout{2} = phi_out_pump;
    varargout{3} = z_mm;
end

end

