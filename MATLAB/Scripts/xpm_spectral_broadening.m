c = 3e8;

pump_pulsewidth_fs_HW_einv_max = 70;
probe_pulsewidth_fs_HW_einv_max = 40/(2*sqrt(log(2)));

alpha_2_cm_per_GW = 25;
n2_m2_per_W = 1e-15;


num_t_points = 2^10;
t_max_fs = 100*probe_pulsewidth_fs_HW_einv_max;
sampling_period_fs = (2*t_max_fs)/num_t_points
t = (-t_max_fs+sampling_period_fs:sampling_period_fs:t_max_fs);
phi_probe = exp(-t.^2./probe_pulsewidth_fs_HW_einv_max.^2);
figure(1)
plot(t, phi_probe);

spectrum_phi_probe = fftshift(fft(phi_probe));
f_Hz = 1e15*1/sampling_period_fs*(-num_t_points/2+1:num_t_points/2)/num_t_points;
figure(2)

lambda_filter_HW_einv_nm = 20/(2*sqrt(log(2)));
lambda_center_nm = 1200;

frequency_width_filter_HW_einv_Hz = c*lambda_filter_HW_einv_nm*1e-9/(lambda_center_nm*1e-9).^2;
offset_Hz = c*0e-9/(lambda_center_nm*1e-9)^2;

filter_spectrum = exp(-(f_Hz - offset_Hz).^2./frequency_width_filter_HW_einv_Hz.^2);

filtered_probe_spectrum = abs(spectrum_phi_probe).*filter_spectrum;

lambda_nm = 1200;
lambda = 1e-9*lambda_nm^2*f_Hz/3e8;
figure(2)
plot(lambda, abs(spectrum_phi_probe))

figure(3)
plot(f_Hz, filter_spectrum)
hold on
plot(f_Hz, abs(spectrum_phi_probe)./max(abs(spectrum_phi_probe)))
plot(f_Hz,  abs(filtered_probe_spectrum)./max(abs(filtered_probe_spectrum)))
hold off

filtered_probe_time_domain = fftshift(ifft(filtered_probe_spectrum));
figure(4)
plot(t, abs(filtered_probe_time_domain))

