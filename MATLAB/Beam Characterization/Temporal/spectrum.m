c = 3e8;

data = dlmread('spec1440nm.txt');

wavelength_nm = data(:,1);
amplitude = data(:,2);

[maximum, center_index] = max(amplitude);
offset = min(amplitude);
center_wavelength = wavelength_nm(center_index);
fwhm_wavelengths = wavelength_nm(amplitude >= maximum/2);
width = abs(fwhm_wavelengths(end) - fwhm_wavelengths(1));


indices = find(wavelength_nm > center_wavelength - 4*width & wavelength_nm < center_wavelength + 4*width);
wavelength = wavelength_nm(indices);
amplitude = amplitude(indices);

fitfun = @(p, xdata)p(1) + p(2)*exp(-2.77*(xdata - p(3)).^2./p(4)^2);


xfit = linspace(wavelength(1), wavelength(end), 1000);
fit_parameters = lsqcurvefit(fitfun, [offset, maximum, center_wavelength, width], wavelength, amplitude);


plot(wavelength, amplitude, 'b')
xlabel('wavelength (nm)'), ylabel('amplitude (a.u.)'), title('Spectrum from Ocean Optics NIR spectrometer');
hold on
plot(xfit, fitfun(fit_parameters,xfit), 'r')
xlim([min(wavelength), max(wavelength)])
hold off

bandwidth_limit = 0.44*(1e-9*fit_parameters(3))^2/(3e8*(1e-9*fit_parameters(4)))*1e15
center_wavelength = fit_parameters(3)
bandwidth_fwhm = fit_parameters(4)
