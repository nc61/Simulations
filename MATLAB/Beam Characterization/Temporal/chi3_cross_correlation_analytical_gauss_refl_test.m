pump_energy_J = 4.18e-6;
alpha_2_cm_per_GW = 1.43;
sample_thickness_m = 280e-6;
pump_spot_size_m_HWeM = 407e-6;
probe_spot_size_m_HWeM = 75.7e-6;
probe_pulsewidth_fs_HWeM = 135;
pump_pulsewidth_fs_HWeM = 42;
pump_group_index = 3.4441;
probe_group_index = 3.7620;
%probe_group_index = 3.4440;
pump_refractive_index = 3.4258;

scale = 8;
delay_fs = probe_pulsewidth_fs_HWeM*linspace(-scale, scale, 200);
delay_shift_fs = 0;
sc = 1;

% 
% [normalized_transmission] = chi3_cross_correlation_analytical_sech2_refl(delay_fs, delay_shift_fs, pump_refractive_index, pump_group_index, probe_group_index, pump_spot_size_m_HWeM, ...
%     probe_spot_size_m_HWeM, pump_pulsewidth_fs_HWeM, probe_pulsewidth_fs_HWeM, sample_thickness_m, pump_energy_J, alpha_2_cm_per_GW);

fit_normalized_transmission = @(fit_parameters, delay_fs) chi3_cross_correlation_analytical_gauss_refl(delay_fs, fit_parameters(3), pump_refractive_index, pump_group_index, probe_group_index, pump_spot_size_m_HWeM, ...
   probe_spot_size_m_HWeM, fit_parameters(2), probe_pulsewidth_fs_HWeM, sample_thickness_m, pump_energy_J, fit_parameters(1), fit_parameters(4), fit_parameters(5), 1);

fit_normalized_transmission_no_gvd = @(fit_parameters, delay_fs) chi3_cross_correlation_analytical_gauss_refl(delay_fs, fit_parameters(3), pump_refractive_index, pump_group_index, 0.999*pump_group_index, pump_spot_size_m_HWeM, ...
   probe_spot_size_m_HWeM, fit_parameters(2), probe_pulsewidth_fs_HWeM, sample_thickness_m, pump_energy_J, fit_parameters(1), fit_parameters(4), fit_parameters(5), 1);

[scan_delay_fs, scan_norm_trans] = read_scan('plot', 'off', 'normalizeToProbe', 'on', 'probeReading', 6.2, 'probeSensitivity', 20, 'scanSensitivity', 20, 'centerDataToPeak', 'on', 'micrometerUnits', 'in');
scan_delay_fs = -scan_delay_fs;
options = optimoptions('lsqcurvefit','FunctionTolerance', 1e-9);
fit_params = lsqcurvefit(fit_normalized_transmission, [1.55, 80, 0, 1, 0], scan_delay_fs, scan_norm_trans);

%fit_params = [1.55, 80e-15, 0];
[fit_results_total, fit_results_no_refl, reflection_results] = fit_normalized_transmission(fit_params, scan_delay_fs);
[fit_results_without_gvd, no_gvd_no_refl, no_gvd_refl] = fit_normalized_transmission_no_gvd(fit_params, scan_delay_fs);

scan_delay_fs = scan_delay_fs - fit_params(3);
plot(scan_delay_fs, fit_results_total, 'b', 'LineWidth', 2)
hold on
plot(scan_delay_fs, scan_norm_trans, 'ro')
plot(scan_delay_fs, fit_results_without_gvd, 'k--')
title_string = strcat('\alpha_2',  sprintf(' = %.2fcm/GW, E_{pump} = %.2f', fit_params(1), pump_energy_J*1e6), '\muJ, ', '\tau_{pump, fit} = ', sprintf('%.0ffs FWHM', fit_params(2)*2*sqrt(log(2))))
xlabel('delay (fs)'), ylabel('\DeltaT'), title(title_string);
hold off

legend('best fit', 'data', 'without gvd');

fwhm_gvd = get_fwhm(scan_delay_fs, fit_results_total)
fwhm_no_gvd = get_fwhm(scan_delay_fs, fit_results_without_gvd)

% plot(delay_fs, normalized_transmission);
% hold on
% hold off
% plot(delay_fs, reflection, 'r');
% plot(delay_fs, reflection + normalized_transmission, 'g');
% hold off