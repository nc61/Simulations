pump_waist_mm = 0.3;
probe_waist_mm = 0.1;

probe_pulsewidth_fs = 150;
pump_pulsewidth_fs = 150;

pump_peak_power_W = 1;
probe_peak_power_W = 1;

sample_thickness_mm = 15;
num_z_steps = 30;

alpha_2_mm_per_W = 0e-8;
alpha_2_d_mm_per_W = 0e-8;

% Number of fourier modes
xmax_mm = 20*pump_waist_mm;
tmax_fs = 20*pump_pulsewidth_fs;

% Number of points to split up x into
num_x_points = 2^8;
num_time_points = 2^8;

% points for x and y 
x_mm = (xmax_mm/num_x_points)*(0:num_x_points-1);
t_fs = (tmax_fs/num_time_points)*(0:num_time_points-1);
[X, T] = meshgrid(x_mm, t_fs);


phi_pump = sqrt(pump_peak_power_W).*exp(-((X - 0.5*xmax_mm).^2./pump_waist_mm^2 + (T - 0.5*tmax_fs).^2./pump_pulsewidth_fs.^2));
phi_probe = sqrt(probe_peak_power_W).*exp(-((X - 0.5*xmax_mm).^2./probe_waist_mm^2 + (T - 0.5*tmax_fs).^2./probe_pulsewidth_fs.^2));

tic
[normalized_probe_energy_out, phi_out_probe, phi_out_pump, z_mm] = xpm_1plus1d_fcn(x_mm, t_fs, phi_pump, phi_probe, sample_thickness_mm, num_z_steps, alpha_2_mm_per_W, ... 
    alpha_2_d_mm_per_W, 3, 3, 1.5, 1.2, 4, 4, 000, 000, 0);
toc

tic
normalized_probe_energy_out = xpm_1plus1d_fcn(x_mm, t_fs, phi_pump, phi_probe, sample_thickness_mm, num_z_steps, alpha_2_mm_per_W, ... 
    alpha_2_d_mm_per_W, 3, 3, 1.5, 1.2, 4, 4.05, 5000, 5000, 1);
toc

normalized_probe_energy_out

figure(1)
% surf(x,z,abs(phi_out(:,:,N_time/2)).^2)
surf(x_mm,z_mm,abs(squeeze(phi_out_pump(:,:,num_time_points/2+1))).^2)
hold on
surf(x_mm,z_mm,abs(squeeze(phi_out_probe(:,:,num_time_points/2+1))).^2)
hold off
shading interp, colormap('jet')

figure(2)
surf(t_fs,z_mm,abs(squeeze(phi_out_pump(:,num_x_points/2+1,:))).^2)
shading interp
colormap('jet')
colorbar
xlabel('X (a.u.)')
ylabel('Distance (a.u.)')
%xlim([min(x) max(x)])
