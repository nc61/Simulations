wavelength_um = linspace(1, 2.2, 500);
[neff, beta1, beta2] = neff_dispersion('QW_dispersion_2um_ridge.mat', wavelength_um);

figure(1)

p1 = subplot(3,1,1);
plot(wavelength_um, neff);
grid(p1, 'on')
title('Effective index')
xlabel('\lambda [\mu m]'), ylabel('n_{eff}'), xlim([min(wavelength_um), max(wavelength_um)])
p2 = subplot(3,1,2);
plot(wavelength_um, 1./beta1*1e12);
grid(p2, 'on')
title('Group velocity')
xlabel('\lambda [\mu m]'), ylabel('v_g [m/s]'), xlim([min(wavelength_um), max(wavelength_um)])
p3 = subplot(3,1,3);
plot(wavelength_um, beta2);
grid(p3, 'on')

title('Group velocity dispersion')
xlabel('\lambda [\mu m]'), ylabel('\beta_2 [fs^2/mm]'), xlim([min(wavelength_um), max(wavelength_um)])
