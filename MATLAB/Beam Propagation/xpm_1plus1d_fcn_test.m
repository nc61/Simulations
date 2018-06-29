
sample_thickness_mm = 5;
num_z_steps = 20;
num_space_points = 2^8;
num_time_points = 2^8;
alpha_2_mm_per_W = 20e-8;
alpha_2_d_mm_per_W = 50e-8;
pump_waist = 0.3;
probe_waist = 0.1;
pump_pulsewidth = 150;
probe_pulsewidth = 150;

xmax_mm = 10*max(pump_waist, probe_waist);
tmax_fs = 10*max(pump_pulsewidth, probe_pulsewidth);

pump_peak_power = 1;
probe_peak_power = 1;

[phi_out_probe, phi_out_pump, x, z, t] = xpm_1plus1d_fcn(xmax_mm, tmax_fs, sample_thickness_mm, num_z_steps, num_space_points, num_time_points, alpha_2_mm_per_W, alpha_2_d_mm_per_W, ... 
                                                pump_waist, probe_waist, pump_pulsewidth, probe_pulsewidth, pump_peak_power, probe_peak_power);

figure(1)
% surf(x,z,abs(phi_out(:,:,N_time/2)).^2)
surf(x,z,abs(squeeze(phi_out_pump(:,:,num_time_points/2+1))).^2)
hold on
surf(x,z,abs(squeeze(phi_out_probe(:,:,num_time_points/2+1))).^2)
hold off
shading interp, colormap('jet')

figure(2)
surf(t,z,abs(squeeze(phi_out_pump(:,num_space_points/2+1,:))).^2)
shading interp
colormap('jet')
colorbar
xlabel('X (a.u.)')
ylabel('Distance (a.u.)')
%xlim([min(x) max(x)])