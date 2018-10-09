c = 3e8;

probe_pulsewidth_fs = 200/(2*sqrt(log(2)));
pump_pulsewidth_fs = 200/(2*sqrt(log(2)));

sample_thickness_mm = 1;
num_z_steps = 30;

E_pump = 1e-12;
E_probe = 1e-12;

pump_peak_power_W = E_pump/(pump_pulsewidth_fs*1e-15*sqrt(pi));
probe_peak_power_W = E_probe/(probe_pulsewidth_fs*1e-15*sqrt(pi));

gamma_ee_minv_Winv =  0 + 1i*0;
gamma_ep_minv_W_inv = -1000 + 0i;


% Calculate group velocity
beta1_fs_per_mm_pump = 1.1209e4;
beta1_fs_per_mm_probe = 1.18778e4;

gvm_fs_per_mm = 668.26;

time_shift_fs = -gvm_fs_per_mm*sample_thickness_mm;

delay_fs = linspace(300, 300, 1);
signal = zeros(size(delay_fs));

for ind = 1:length(delay_fs)
    
% Number of fourier modes
tmax_fs = 10*pump_pulsewidth_fs + (2*abs(time_shift_fs) + abs(delay_fs(ind)));

% Number of points to split up x into
num_time_points = 2^9;

% points for x and y 
t_fs = (tmax_fs/num_time_points)*(0:num_time_points-1);


phi_pump = sqrt(pump_peak_power_W).*exp(-(t_fs - delay_fs(ind) - 0.5*tmax_fs).^2./(2*pump_pulsewidth_fs.^2));
phi_probe = sqrt(probe_peak_power_W).*exp(-(t_fs - 0.5*tmax_fs).^2./(2*probe_pulsewidth_fs.^2));


[normalized_probe_energy_out, phi_out_probe, phi_out_pump, z_mm] = xpm_1d_fcn(t_fs, phi_pump, phi_probe, sample_thickness_mm, num_z_steps, gamma_ep_minv_W_inv, ... 
    gamma_ee_minv_Winv, gvm_fs_per_mm, 8000, 8000, 0);


signal(ind) = normalized_probe_energy_out;
end

ratio = real(max(max(abs(phi_out_pump).^2))/max(max(abs(phi_out_probe).^2)));

figure(2)
surf(t_fs,z_mm,abs(phi_out_pump).^2)
hold on
surf(t_fs,z_mm,abs(phi_out_probe).^2*ratio)
hold off
shading interp
colormap('jet')
colorbar
xlabel('T (fs)')
ylabel('Distance (mm)')
%xlim([min(x) max(x)])

figure(3)
plot(delay_fs, signal)
xlabel('delay (fs)'), ylabel('\Delta T')
