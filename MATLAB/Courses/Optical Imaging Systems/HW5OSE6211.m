N = 128;
n = 0:N-1;
M = 8;
f = (-1).^n.*(exp(-(n-N/2).^2./M.^2));
F = fft(f);

freq = (n)/N;

figure(1)
plot(freq, abs(F));
xlabel('normalized frequency'), title('|F|')

figure(2)
plot(abs(f))
xlabel('n'), title('|f|')

sum(abs(f).^2)-sum(abs(F).^2)*(freq(2) - freq(1))