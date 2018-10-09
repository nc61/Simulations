
pump_pulsewidth_fs = 200;
probe_pulsewidth_fs = 200;
sample_thickness_mm = 5;
pump_peak_power_W = 1;
probe_peak_power_W = 1;
pump_energy = 1;

pump_waist_m = .3;
probe_waist_m = .7;

alpha_2_mm_per_W = 2500000e-8;

% Number of fourier modes
tmax_fs = 12*pump_pulsewidth_fs;

% Number of points to split up x into
num_time_points = 1000;

% points for x and y 
t_fs = linspace(-tmax_fs, tmax_fs, 200);
width = sqrt(pump_pulsewidth_fs^2 + probe_pulsewidth_fs^2);
delay_fs = linspace(-3*width, 3*width, 100); 

phi_pump = sqrt(pump_peak_power_W).*exp(-(t_fs).^2./(2*pump_pulsewidth_fs.^2));
phi_probe = sqrt(probe_peak_power_W).*exp(-(t_fs).^2./(2*probe_pulsewidth_fs.^2));
E_probe_initial = trapz(t_fs, abs(phi_probe).^2)*pi*probe_waist_m^2;

phi_probe_out = exp(-sample_thickness_mm*phi_pump.^2*alpha_2_mm_per_W*(1/(probe_waist_m^2 + pump_waist_m^2)));

cross_correlation = zeros(size(delay_fs));

for ind = 1:length(delay_fs)
    intensity_out = abs(phi_probe_out.*exp(-(t_fs - delay_fs(ind)).^2/(2*probe_pulsewidth_fs^2))).^2;
    cross_correlation(ind) = trapz(t_fs, intensity_out)*pi*probe_waist_m^2;
end

norm_trans_method1 = cross_correlation/E_probe_initial;
normalized_transmission = 1 - 2*alpha_2_mm_per_W*pump_energy*sample_thickness_mm*(pump_pulsewidth_fs^2/(pump_pulsewidth_fs^2 + probe_pulsewidth_fs^2))^(1/2)*exp(-delay_fs.^2/(pump_pulsewidth_fs^2 + probe_pulsewidth_fs^2)) ...
                           *1/(pump_waist_m^2 + probe_waist_m^2);
percent_error = 100*(min(normalized_transmission) - min(norm_trans_method1))/min(norm_trans_method1);


figure(1)
plot(delay_fs, norm_trans_method1)
hold on
plot(delay_fs, normalized_transmission, 'r')
hold off




