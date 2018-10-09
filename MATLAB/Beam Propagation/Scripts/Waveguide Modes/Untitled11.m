phase = linspace(-5*pi,5*pi,20000);


I = 1/2*(1 - cos(pi/2 + pi*sin(phase)));

figure(1)
plot(1./phase, fftshift(fft(abs(I).^2)))