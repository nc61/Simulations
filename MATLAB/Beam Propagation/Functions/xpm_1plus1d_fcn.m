function [normalized_probe_energy_out, varargout] = xpm_1plus1d_fcn(x_mm, t_fs, phi_pump, ... 
    phi_probe, sample_thickness_mm, num_z_steps, alpha_2_mm_per_W, alpha_2_d_mm_per_W, ... 
    pump_refractive_index, probe_refractive_index, pump_wavelength_um, probe_wavelength_um, ... 
    pump_group_index, probe_group_index, ...
    pump_dispersion_fs2_per_mm, probe_dispersion_fs2_per_mm, output_normalized_energy_only)

c = 2.998e8;

xmax_mm = max(x_mm);
tmax_fs = max(t_fs);

num_time_points = length(t_fs);
num_x_points = length(x_mm);

% Fourier space for x and y
kx = 2*pi/xmax_mm*[(0:num_x_points/2-1) (-num_x_points/2:-1)];
omega = 2*pi/tmax_fs*[(0:num_time_points/2-1) (-num_time_points/2:-1)];

[Kx, OMEGA] = ndgrid(kx,omega);

% nonlinear
Aeff_mm_squared = 20e-6;

gamma2 = alpha_2_mm_per_W/(2*Aeff_mm_squared);
gamma2_d = alpha_2_d_mm_per_W/(2*Aeff_mm_squared);

% Propagation distance
dz_mm = sample_thickness_mm/num_z_steps;
dx_m = 1e-3*xmax_mm/num_x_points;
dt_s = 1e-15*tmax_fs/num_time_points;
initial_probe_energy_W =  dx_m*dt_s*sum(sum(abs(phi_probe(:,:)).^2));

% Calculation of propagation constants
k_pump_per_mm = 1e3*pump_refractive_index*2*pi/pump_wavelength_um;
k_probe_per_mm = 1e3*probe_refractive_index*2*pi/probe_wavelength_um;

% Calculate group velocity
pump_group_velocity_mm_per_fs = 1e-12*c/pump_group_index;
probe_group_velocity_mm_per_fs = 1e-12*c/probe_group_index;

group_velocity_dispersion_fs_per_mm = 1/pump_group_velocity_mm_per_fs - 1/probe_group_velocity_mm_per_fs;

% Linear propagation operator
D_pump = -((-1i*1/(2*k_pump_per_mm)*Kx.^2 + 1i*pump_dispersion_fs2_per_mm*OMEGA.^2 - 1i*(group_velocity_dispersion_fs_per_mm)*OMEGA)*0.5*dz_mm);
D_probe = -((-1i*1/(2*k_probe_per_mm)*Kx.^2 + 1i*probe_dispersion_fs2_per_mm*OMEGA.^2)*0.5*dz_mm);


if (~output_normalized_energy_only)
    phi_out_probe = zeros(num_z_steps, num_x_points, num_time_points);
    phi_out_pump = zeros(num_z_steps, num_x_points, num_time_points);
end

z_mm = dz_mm:dz_mm:sample_thickness_mm;

% Not symmetrized, will need to do that for better accuracy
for ind = 1:num_z_steps
    
    %Pump
    lin_step_pump_1 = ifft2(fft2(phi_pump).*exp(D_pump));
    nonlin_step_pump = lin_step_pump_1.*exp((-2*gamma2_d*abs(lin_step_pump_1).^2)*dz_mm);
    lin_step_pump_2 = ifft2(fft2(nonlin_step_pump).*exp(D_pump));
    phi_pump = lin_step_pump_2;
    
    % Probe
    lin_step_probe_1 = ifft2(fft2(phi_probe).*exp(D_probe));
    nonlin_step_probe = lin_step_probe_1.*exp((-2*gamma2*abs(lin_step_pump_1).^2)*dz_mm);
    lin_step_probe_2 = ifft2(fft2(nonlin_step_probe).*exp(D_probe));
    phi_probe = lin_step_probe_2;
    
    if (~output_normalized_energy_only)
        phi_out_pump(ind, :, :) = phi_pump;
        phi_out_probe(ind,:,:) = phi_probe;
    end
end

final_probe_energy_W =  dx_m*dt_s*sum(sum(abs(phi_probe(:,:)).^2));
normalized_probe_energy_out = (initial_probe_energy_W - final_probe_energy_W)/initial_probe_energy_W;

if ~output_normalized_energy_only
    varargout{1} = phi_out_probe;
    varargout{2} = phi_out_pump;
    varargout{3} = z_mm;
end

end

