pump_waist_mm = 0.01;
probe_waist_mm = 0.01;

probe_pulsewidth_fs = 150;
pump_pulsewidth_fs = 150;

pump_peak_intensity_GW_per_cm2 = 1;
probe_peak_intensity_GW_per_cm2 = 1;

sample_thickness_mm = 1;
num_z_steps = 30;

alpha_2_nd_cm_per_GW = 200e-8;
alpha_2_d_cm_per_GW = 50e-8;

% Number of fourier modes
xmax_mm = 10*pump_waist_mm;
ymax_mm = 10*pump_waist_mm;
tmax_fs = 10*pump_pulsewidth_fs;

% Number of points to split up x into
num_x_points = 2^5;
num_y_points = 2^5;
num_time_points = 2^5;

% points for x and y 
x_mm = (xmax_mm/num_x_points)*(0:num_x_points-1);
y_mm = (ymax_mm/num_y_points)*(0:num_y_points-1);
t_fs = (tmax_fs/num_time_points)*(0:num_time_points-1);
[X, Y, T] = meshgrid(x_mm, y_mm, t_fs);


phi_pump = sqrt(pump_peak_intensity_GW_per_cm2).*exp(-((X - 0.5*xmax_mm).^2./pump_waist_mm^2 + (Y - 0.5*ymax_mm).^2./pump_waist_mm^2 + (T - 0.5*tmax_fs).^2./pump_pulsewidth_fs.^2));
phi_probe = sqrt(probe_peak_power_GW_per_cm2).*exp(-(((X - 0.5*xmax_mm).^2 + (Y - 0.5*ymax_mm).^2)./probe_waist_mm^2 + (T - 0.5*tmax_fs).^2./probe_pulsewidth_fs.^2));

tic
[normalized_probe_energy_out, phi_out_probe, phi_out_pump, z_mm] = xpm_1plus2d_fcn(x_mm, y_mm, t_fs, phi_pump, phi_probe, sample_thickness_mm, num_z_steps, alpha_2_nd_cm_per_GW, ... 
    alpha_2_d_mm_per_W, 3, 3, 1.2, 1.2, 4, 4, 000, 000, 0);
toc

tic
normalized_probe_energy_out = xpm_1plus2d_fcn(x_mm, y_mm, t_fs, phi_pump, phi_probe, sample_thickness_mm, num_z_steps, alpha_2_nd_cm_per_GW, ... 
    alpha_2_d_mm_per_W, 3, 3, 1.5, 1.2, 4, 4.2, 5000, 5000, 1);
toc

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
