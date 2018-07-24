pump_waist_mm = 1.24;
probe_waist_mm = 0.0857;

probe_pulsewidth_fs = 182;
pump_pulsewidth_fs = 200;

pump_energy_J = 7e-6;
pump_peak_intensity_GW_per_cm2 = gaussianbeam_peak_intensity_fcn(pump_energy_J, pump_waist_mm, pump_pulsewidth_fs)
probe_peak_intensity_GW_per_cm2 = 0.01;

sample_thickness_mm = 1;
num_z_steps = 7;

alpha_2_nd_cm_per_GW = 30;
alpha_2_d_cm_per_GW = 0;

% Number of fourier modes
xmax_mm = 20*pump_waist_mm;
ymax_mm = 20*pump_waist_mm;
tmax_fs = 30*pump_pulsewidth_fs;

% Number of points to split up x into
num_x_points = 2^5;
num_y_points = 2^5;
num_time_points = 2^6;

% points for x and y 
x_mm = (xmax_mm/num_x_points)*(0:num_x_points-1);
y_mm = (ymax_mm/num_y_points)*(0:num_y_points-1);
t_fs = (tmax_fs/num_time_points)*(0:num_time_points-1);
[X_mm, Y_mm, T_fs] = meshgrid(x_mm, y_mm, t_fs);

pump_wavelength_um = 4;
pump_refractive_index = 3.4258;
pump_group_index = 3.4441;
pump_group_velocity_dispersion_fs2_per_mm = 389.25;

probe_wavelength_um = 1.25;
probe_refractive_index = 3.513;
probe_group_index = 3.7005;
probe_group_velocity_dispersion_fs2_per_mm = 2*389.25;


phi_pump = sqrt(pump_peak_intensity_GW_per_cm2).*exp(-((X_mm - 0.5*xmax_mm).^2./(2*pump_waist_mm^2) + (Y_mm - 0.5*ymax_mm).^2./(2*pump_waist_mm^2) + (T_fs - 0.4*tmax_fs).^2./(2*pump_pulsewidth_fs.^2)));
phi_probe = sqrt(probe_peak_intensity_GW_per_cm2).*exp(-(((X_mm - 0.5*xmax_mm).^2 + (Y_mm - 0.5*ymax_mm).^2)./(2*probe_waist_mm^2) + (T_fs - 0.5*tmax_fs).^2./(2*probe_pulsewidth_fs.^2)));

% [normalized_probe_energy_out, phi_out_probe, phi_out_pump, z_mm] = xpm_1plus2d_fcn(x_mm, y_mm, t_fs, phi_pump, phi_probe, sample_thickness_mm, num_z_steps, alpha_2_nd_cm_per_GW, ... 
%     alpha_2_d_cm_per_GW, pump_refractive_index, probe_refractive_index, pump_wavelength_um, probe_wavelength_um, pump_group_index, ... 
%     probe_group_index, pump_group_velocity_dispersion_fs2_per_mm, probe_group_velocity_dispersion_fs2_per_mm, 0);

signal = normalized_probe_energy_out


figure(1)
% surf(x,z,abs(phi_out(:,:,N_time/2)).^2)
surf(x_mm,z_mm,abs(squeeze(phi_out_pump(:,:,num_y_points/2+1,num_time_points/2+1))).^2)
hold on
surf(x_mm,z_mm,abs(squeeze(phi_out_probe(:,:,num_y_points/2+1,num_time_points/2+1))).^2)
hold off
shading interp, colormap('jet')

figure(2)
surf(t_fs,z_mm,abs(squeeze(phi_out_pump(:,num_x_points/2+1,num_y_points/2+1,:))).^2)
hold on
surf(t_fs,z_mm,abs(squeeze(phi_out_probe(:,num_x_points/2+1,num_y_points/2+1,:))).^2)
hold off
shading interp
colormap('jet')
colorbar
xlabel('X (a.u.)')
ylabel('Distance (a.u.)')
%xlim([min(x) max(x)])