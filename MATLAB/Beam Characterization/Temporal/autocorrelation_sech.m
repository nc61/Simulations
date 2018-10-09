c = 3e8;

filename = file_selector();

[stage_position, signal] = get_scan_data(filename);

fitfun = @(p, xdata)p(1) + p(2)*exp(-(xdata - p(3)).^2./p(4)^2);

[maximum, center_index] = max(signal);
offset = mean(signal(1:5));
amplitude = maximum - offset;
zero_delay = stage_position(center_index);
fwhm_positions = stage_position(signal - offset >= amplitude/2);
pulsewidth = abs(fwhm_positions(end) - fwhm_positions(1));

fit_parameters = lsqcurvefit(fitfun, [offset, amplitude, zero_delay, pulsewidth], stage_position, signal);

xfit = linspace(stage_position(1), stage_position(end), 1000);
yfit = fitfun(fit_parameters,xfit);

pulsewidth = 2*fit_parameters(4)/10*2.54/100/c*1e15/sqrt(2);

[~, fit_center_index] = max(yfit);
xfit_shifted = xfit - xfit(fit_center_index);
xfit_delay_shifted = 2*xfit_shifted/10*2.54/100/c*1e15;

data_pos_shifted = stage_position - xfit(fit_center_index);
delay = 2*data_pos_shifted/10*2.54/100/c*1e15;
plot(delay, signal, 'bo')
xlabel('delay (fs)'), ylabel('signal (a.u.)'), title(sprintf('Autocorrelation signal, t_p = %f', pulsewidth))

hold on
plot(xfit_delay_shifted, fitfun(fit_parameters,xfit), 'r');
hold off

