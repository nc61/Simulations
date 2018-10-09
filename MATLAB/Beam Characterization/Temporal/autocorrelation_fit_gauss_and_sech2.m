[delay_fs, signal] = read_scan('micrometerUnits', 'in', 'centerDataToPeak', 'on', 'convertToDelay', 'on');

fitfun_gaussian = @(p, xdata)p(1) + p(2)*exp(-(xdata - p(3)).^2./(2*p(4)^2));
fit_params_gaussian = lsqcurvefit(fitfun_gaussian, [pulsewidth, zero_delay, offset, amplitude], delay_fs, signal);
fit_params_sech2 = lsqcurvefit(@sech2_fit_function, [pulsewidth, zero_delay, offset, amplitude], delay_fs, signal);

gauss_fit = fitfun_gaussian(fit_params_gaussian, delay_fs);
sech2_fit = sech2_fit_function(fit_params_sech2, delay_fs);

figure(1)
plot(delay_fs, gauss_fit)
hold on
plot(delay_fs, sech2_fit)
hold off

gauss_fwhm = get_fwhm(delay_fs, gauss_fit)
sech2_fwhm = get_fwhm(delay_fs, sech2_fit)

figure(2)
gaussian_pulse = exp(-delay_fs.^2/fit_params_gaussian(4)^2);
sech2_pulse = sech(delay_fs/fit_params_sech2(1)).^2;
plot(delay_fs, gaussian_pulse)
hold on
plot(delay_fs, sech2_pulse)
hold off

get_fwhm(delay_fs, gaussian_pulse)
get_fwhm(delay_fs, sech2_pulse)