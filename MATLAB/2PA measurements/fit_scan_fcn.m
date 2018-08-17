function fit_params = fit_scan_fcn(probe_value, probe_sensitivity, scan_sensitivity, pump_energy_J, sample_thickness_m)

filename = file_selector();
probe = probe_value*probe_sensitivity;
scan_sens = scan_sensitivity;

scan_data = tdfread(filename);
    stage_position = scan_data.Position;
    absorp = scan_sens*scan_data.x1;
        
    [~, zero_delay_index] = max(absorp);
    delay = 2*(stage_position - stage_position(zero_delay_index))/10*2.54/100/2.998e8*1e15; % [fs]
    absorption = 1+(absorp - mean(absorp(1:20)))/probe;
    
    
probe_pulse_width_HW_einv_max_fs = 137;
pump_spot_size_x_HW_einv2_max_m = 530*1e-6;
pump_spot_size_y_HW_einv2_max_m = 626*1e-6;
probe_spot_size_HW_einv2_max_m = 107*1e-6;

pump_probe_fit = @(fit_parameters, t_d)pump_probe_fcn(fit_parameters(1), fit_parameters(2), fit_parameters(3), probe_pulse_width_HW_einv_max_fs, pump_energy_J, sample_thickness_m, pump_spot_size_x_HW_einv2_max_m, pump_spot_size_y_HW_einv2_max_m, probe_spot_size_HW_einv2_max_m, t_d); 

options = optimset('TolX', 1e-15);
fit_params = lsqcurvefit(pump_probe_fit, [1.5 ,45, 600],  delay, absorption, [1, 20, 0], [5, 200, 1600], options);
%fit_params = [1.4,45,600];
trans_fit = pump_probe_fit(fit_params, delay);

[~,peak_index] = min(trans_fit);
delay = delay - delay(peak_index);

figure(1)
plot(delay, absorption, 'ro')
hold on
plot(delay, trans_fit, 'b', 'Linewidth', 2)
hold off
xlabel('delay (fs)'), ylabel('T (a.u.)'), xlim([min(delay) max(delay)]), title(strcat(sprintf('Si Pump = 4000nm, Probe = 1200nm, E_{pump} = %.2f', pump_energy_J*1e6), '\muJ'))