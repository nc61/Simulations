d = 10;
R1 = -50;
R2 = -100;

M1 = [1 d; 0 1];
M2 = [1 0; 2/R2 1];
M3 = [1 d; 0 1];
M4 = [1 0; 2/R1 1];

M = M4*M3*M2*M1;

n = 0:10;
y_out = zeros(size(n));

v_in = [0; 0.1*pi/180];

for ind = 1:length(n)
    v_out = M^n(ind)*v_in;
    y_out(ind) = v_out(1);
end

plot(n, y_out, 'b*');
hold on
n_sin = linspace(min(n), max(n), 200);

b = 1/2*trace(M);
phi = acos(b);
y_sin = d*2*sind(0.1)*sin(n_sin*phi);
plot(n_sin, y_sin, 'r-')
hold off
