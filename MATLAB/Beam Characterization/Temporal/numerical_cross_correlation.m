
t1 = 10;
t2 = 10;
tau_max = 4*sqrt(t1^2 + t2^2);
tau = linspace(-tau_max,tau_max,100);
t = linspace(-tau_max,tau_max,100);

signal = zeros(size(tau));
for ind = 1:length(tau)
    y1 = sech((t - tau(ind))/t1).^2;
    y2 = sech(t/t2).^2;
    signal(ind) = trapz(t, y1.*y2);
end

plot(tau, signal)