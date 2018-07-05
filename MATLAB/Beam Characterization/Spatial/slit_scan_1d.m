energy_J = 0.5e-6;
pulsewidth_fs_HW_inve_max = 222;
waist_mm_HW_inve_max = 4;

peak_intensity_GW_per_cm2 = (1e-9)*energy_J/(pi^(3/2)*(1e-1*waist_mm_HW_inve_max)^2*(1e-15)*pulsewidth_fs_HW_inve_max);

num_x_points = 2^12;
xmax_mm = 20*waist_mm_HW_inve_max;
x_mm = (xmax_mm/num_x_points)*(0:num_x_points-1);

kx = 2*pi/num_x_points*[(0:num_x_points/2-1) (-num_x_points/2:-1)];

phi = sqrt(peak_intensity_GW_per_cm2)*sech(-((x_mm - 0.5*xmax_mm))/(0.1*waist_mm_HW_inve_max^2)).^2;
phi_fourier = fft(phi);

mask = rectangularPulse(45,55,x_mm);
mask_fourier = fft(mask);


figure(1)
plot(x_mm, phi)
hold on
plot(x_mm, mask)
hold off

figure(2)
plot(kx, phi_fourier.*exp(1i.*kx.*pi/2*xmax_mm))


figure(3)

output = ifft(mask_fourier.*abs(phi_fourier).^2);
plot(x_mm, output)

figure(4)
output_fft = fft(output)
plot(kx, output_fft)

signal_fourier = output_fft./mask_fourier;
signal_fourier = signal_fourier(~isnan(signal_fourier));
kx_nan_removed = kx(~isnan(signal_fourier))

ifft(signal_fourier)
plot(signal)






