pump_waist_mm = 0.3;
probe_waist_mm = 0.1;

probe_pulsewidth_fs = 222;
pump_pulsewidth_fs = 100;

pump_peak_intensity_GW_per_cm2 = 0.4494;
probe_peak_intensity_GW_per_cm2 = 0.01;

sample_thickness_mm = 1;
num_z_steps = 20;

alpha_2_nd_cm_per_GW = 20;
alpha_2_d_cm_per_GW = 5;

% Number of fourier modes
xmax_mm = 20*pump_waist_mm;
ymax_mm = 20*pump_waist_mm;
tmax_fs = 20*pump_pulsewidth_fs;

% Number of points to split up x into
num_x_points = 2^6;
num_y_points = 2^6;
num_time_points = 2^6;

% points for x and y 
x_mm = (xmax_mm/num_x_points)*(0:num_x_points-1);
y_mm = (ymax_mm/num_y_points)*(0:num_y_points-1);
t_fs = (tmax_fs/num_time_points)*(0:num_time_points-1);
[X_mm, Y_mm, T_fs] = meshgrid(x_mm, y_mm, t_fs);


phi_pump = sqrt(pump_peak_intensity_GW_per_cm2).*exp(-((X_mm - 0.5*xmax_mm).^2./(2*pump_waist_mm^2) + (Y_mm - 0.5*ymax_mm).^2./(2*pump_waist_mm^2) + (T_fs - 0.5*tmax_fs).^2./(2*pump_pulsewidth_fs.^2)));
phi_probe = sqrt(probe_peak_intensity_GW_per_cm2).*exp(-(((X_mm - 0.5*xmax_mm).^2 + (Y_mm - 0.5*ymax_mm).^2)./(2*probe_waist_mm^2) + (T_fs - 0.5*tmax_fs).^2./(2*probe_pulsewidth_fs.^2)));

tic
[normalized_probe_energy_out, phi_out_probe, phi_out_pump, z_mm] = xpm_1plus2d_fcn(x_mm, y_mm, t_fs, phi_pump, phi_probe, sample_thickness_mm, num_z_steps, alpha_2_nd_cm_per_GW, ... 
    alpha_2_d_cm_per_GW, 3.5110, 3.4190, 1.25, 7.75, 3.6985, 3.4239, 301, 201, 0);
toc

% tic
% normalized_probe_energy_out = xpm_1plus2d_fcn(x_mm, y_mm, t_fs, phi_pump, phi_probe, sample_thickness_mm, num_z_steps, alpha_2_nd_cm_per_GW, ... 
%     alpha_2_d_cm_per_GW, 3, 3, 1.5, 1.2, 4.2, 4, 5000, 5000, 1);
% toc

normalized_probe_energy_out

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
