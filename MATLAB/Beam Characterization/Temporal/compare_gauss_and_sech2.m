delay = linspace(-500,500, 500);

gauss_tau = 132;
sech2_tau = 114;

gauss_pulse = exp(-delay.^2/gauss_tau^2);
sech2_pulse = sech(delay/sech2_tau).^2;

plot(delay, gauss_pulse)
hold on
plot(delay, sech2_pulse ,'r')
hold off