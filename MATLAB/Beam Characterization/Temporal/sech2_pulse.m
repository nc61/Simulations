
tau_pump_s = 122e-15;
tau_probe_s = 80e-15;
tau_gvd_s = 500e-15;
alpha_2_m_per_W = 1.55*1e-11;
sample_thickness_m = 280e-6;
pump_energy_J = 5e-6;
pump_waist_m = 500e-6;
probe_waist_m = 100e-6;

z = linspace(0, sample_thickness_m, 100);

scale = 10;
t = scale*linspace(-tau_pump_s, tau_probe_s, 200);
[Z,T] = ndgrid(z,t);

delay = 1e-15*linspace(-1000, 600, 200);
signal = zeros(size(delay));

for ind = 1:length(delay)
G = -1e15*sech(T/tau_pump_s).^2.*sech((T - (Z/L)*tau_gvd_s - delay(ind))/tau_probe_s).^2;
signal(ind) = squeeze(trapz(t, trapz(z, G, 1), 2));
end

plot(delay, signal)


signal2 = zeros(size(delay));
E0 = 2*alpha_2_m_per_W*sample_thickness_m*pump_energy_J/(2*tau_pump_s)*1/(2*pi)*1/(pump_waist_m^2 + probe_waist_m^2)*1/tau_gvd_s;

for ind = 1:length(delay)
    int = sech(t/tau_pump_s).^2.*(tanh((t - delay(ind) - tau_gvd_s)/tau_probe_s) - tanh((t - delay(ind))/tau_probe_s));
    signal2(ind) = E0*trapz(delay, int);
end

hold on
plot(delay, signal2)
hold off