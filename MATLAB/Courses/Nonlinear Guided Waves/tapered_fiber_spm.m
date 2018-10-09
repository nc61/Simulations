% Constants
eps0 = 8.85e-12;
c = 3e8;
eta_0 = 1/(c*eps0);

% Material and beam parameters
n = 1.44;
k0 = 2*pi/(1.5e-6);
n2_I = 3e-20;
n2 = n2_I*n/(2*eta_0);
w = 5e-6;

P_peak = 1;

% Propagation distances and z sampling
L = [0 5e3 10e3];
zsteps = 201;
z = zeros(length(L), zsteps);   

for ind = 1:length(L)
    z(ind, :) = linspace(0, L(ind), zsteps);
end
dz = z(:, 2) - z(:, 1);

% Tau sampling
N = 2^12;
tau_0 = 1e-12;
tau = linspace(-50*tau_0, 50*tau_0, N);

% initialize Power spectral density
phi_omega_abs = zeros(length(z), N);

% FFT frequency
f = 1/(tau(2) - tau(1))*(-N/2:N/2-1)/N;

% Function for calculating waist at disat
wz = @(z)w*exp(-z/1e4);

linecolors = ['r','b','k'];
for ind = 1:length(L)
    w_tapered = wz(z(ind,:));
    E0 = sqrt(2*P_peak./(n*c*eps0*pi*(w_tapered).^2));
    E0(end)
    abs_E0_integrated = sum(E0.^2, 2)*dz(ind)
    phi_squared_integrated = exp(-2*tau.^2./tau_0^2).*abs_E0_integrated;
    phase = k0*n2/2*phi_squared_integrated;
    phi_0 = E0(1)*exp(-tau.^2./tau_0^2);
    figure(1)
    plot(phi_0)
    kernel = phi_0.*exp(1i*phase);
    phi_omega_abs(ind,:) = abs(fftshift(fft(kernel))).^2;
    figure(2)
    plot(f/1e12, phi_omega_abs(ind,:), 'color', linecolors(ind))
    hold on
end
xlabel('f (THz)'), ylabel('PSD'), title('SPM power spectrum in tapered fiber'), legend('1km', '5km', '10km'), xlim([-7 7])
hold off

