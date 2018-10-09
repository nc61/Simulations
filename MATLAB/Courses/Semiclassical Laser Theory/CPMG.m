% Fit the peaks from the CPMG data to get T2

fid_indices = 4:2:num_blocks;
V_echo_peak = zeros(1, length(fid_indices));
t_peak = zeros(1, length(fid_indices));

jnd = 1;
for ind = 4:2:num_blocks
    t_step = t_cellarray{ind};
    t_peak(jnd) = t_step(round(length(t_step)/2));
    V_echo_peak(jnd) = abs(V_integrated(max(find(t_total <= t_peak(jnd)))));
    jnd = jnd + 1;
end

t_peak = [0 t_peak];
V_echo_peak = [1 V_echo_peak];
plot(t_peak*1e6, V_echo_peak, 'k+')
str = fit(t_peak(1:3)', V_echo_peak(1:3)', 'exp1');
t_fit = linspace(0, t_peak(end), 500);
hold on
plot(t_fit*1e6, str.a*exp(str.b*t_fit), 'r')
hold off
title(sprintf('CPMG sequence, T2 = %.1f us', -1/(1e-6*str.b)))
xlabel('t (\mus)'), ylabel('V')