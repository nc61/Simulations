% Constants
eps0 = 8.85e-12;
c = 3e8;
eta_0 = 1/(c*eps0);

% Material/beam parameters
n = 1.44;
k0 = 2*pi/(1.5e-6);
n2_I = 3e-20;
n2 = n2_I*n/(2*eta_0);
w = 5e-6;
tau_0 = 1e-12;
P_peak = 1;

% Propagation distances
z = [0 5e3 10e3];

% Peak electric field
E0 = sqrt(2*P_peak/(n*c*eps0*pi*w^2))
% Sampling of Tau
N = 2^12;
tau = linspace(-50*tau_0, 50*tau_0, N);

% Power spectral density
phi_omega_abs = zeros(length(z), N);

% Define frequency range for FFT
f = 1/(tau(2) - tau(1))*(-N/2:N/2-1)/N;
linecolors = ['r','b','k'];

figure(1)
for ind = 1:length(z)
    phi_0 = E0*exp(-tau.^2./tau_0^2);
    kernel = phi_0.*exp(1i*k0*n2/2*z(ind)*abs(phi_0).^2);
    phi_omega_abs(ind,:) = abs(fftshift(fft(kernel))).^2;
    plot(f/1e12, phi_omega_abs(ind,:), 'color', linecolors(ind))
    hold on
end
xlabel('f (THz)'), ylabel('PSD'), title('Self phase modulation power spectrum'), legend('1km', '5km', '10km'), xlim([-3 3])
hold off
