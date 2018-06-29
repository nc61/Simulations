function [normalized_probe_energy_out, phi_out_probe, phi_out_pump, x_mm, z_mm, t_fs] = xpm_1plus1d_fcn(xmax_mm, tmax_fs, sample_thickness_mm, num_z_steps, num_space_points, num_time_points, alpha_2_mm_per_W, ... 
    alpha_2_d_mm_per_W, pump_waist, probe_waist, pump_pulsewidth_fs, probe_pulsewidth_fs, pump_peak_power_W, probe_peak_power_W)


% points for x and y 
x_mm = (xmax_mm/num_space_points)*(0:num_space_points-1);
t_fs = (tmax_fs/num_time_points)*(0:num_time_points-1);
[X, T] = meshgrid(x_mm,t_fs);

% Fourier space for x and y
kx = 2*pi/xmax_mm*[(0:num_space_points/2-1) (-num_space_points/2:-1)];
omega = 2*pi/tmax_fs*[(0:num_time_points/2-1) (-num_time_points/2:-1)];


[Kx, OMEGA] = ndgrid(kx,omega);

% nonlinear
Aeff_mm_squared = 20e-6;

gamma2 = alpha_2_mm_per_W/(2*Aeff_mm_squared);
gamma2_d = alpha_2_d_mm_per_W/(2*Aeff_mm_squared);

% Propagation distance
dZ = sample_thickness_mm/num_z_steps;

% Initial envelope
phi_pump = sqrt(pump_peak_power_W).*exp(-((X - 0.5*xmax_mm).^2./pump_waist^2 + (T - 0.5*tmax_fs).^2./pump_pulsewidth_fs.^2));
phi_probe = sqrt(probe_peak_power_W).*exp(-((X - 0.5*xmax_mm).^2./probe_waist^2 + (T - 0.5*tmax_fs).^2./probe_pulsewidth_fs.^2));


dX = 1e-3*xmax_mm/num_space_points;
dT = 1e-15*tmax_fs/num_time_points;
initial_probe_energy =  dX*dT*sum(sum(abs(phi_probe(:,:)).^2));


% Linear propagation operator
D_pump = -(0.5*(-1i*0.01*Kx.^2 + 1i*5000*OMEGA.^2)*0.5*dZ);
D_probe = -(0.5*(-1i*0.01*Kx.^2 + 1i*5000*OMEGA.^2)*0.5*dZ);

phi_out_probe = zeros(num_z_steps, num_space_points, num_time_points);
phi_out_pump = zeros(num_z_steps, num_space_points, num_time_points);

z_mm = dZ:dZ:sample_thickness_mm;

% Not symmetrized, will need to do that for better accuracy
for ind = 1:num_z_steps
    
    %Pump
    lin_step_pump_1 = ifft2(fft2(phi_pump).*exp(D_pump));
    nonlin_step_pump = lin_step_pump_1.*exp((-2*gamma2_d*abs(lin_step_pump_1).^2)*dZ);
    lin_step_pump_2 = ifft2(fft2(nonlin_step_pump).*exp(D_pump));
    phi_pump = lin_step_pump_2;
    
    % Probe 
    lin_step_probe_1 = ifft2(fft2(phi_probe).*exp(D_probe));
    nonlin_step_probe = lin_step_probe_1.*exp((-2*gamma2*abs(lin_step_pump_1).^2)*dZ);    
    lin_step_probe_2 = ifft2(fft2(nonlin_step_probe).*exp(D_probe));    
    phi_probe = lin_step_probe_2;
    
    phi_out_pump(ind, :, :) = phi_pump;
    phi_out_probe(ind,:,:) = phi_probe;
end

final_probe_energy =  dX*dT*sum(sum(abs(phi_probe(:,:)).^2));
normalized_probe_energy_out = (initial_probe_energy - final_probe_energy)/initial_probe_energy;

end

