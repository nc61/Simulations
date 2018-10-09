c = 3e8;
C = 9.33e-14;
q = 1.6e-19;
m_e = 9.109e-31;
eps_0 = 8.85e-12;

x_function = @(t, x_0, w_0, w_0t_p)x_0*(cos(w_0*t) - cos(w_0t_p) + sin(w_0t_p)*(w_0*t - w_0t_p));
T = 2.67e-15;
w_0 = 2*pi/T;
I = 5e14;
E_field_V_per_m = sqrt(2*I*100^2/(eps_0*c));

x_0 = q*E_field_V_per_m/(m_e*w_0^2);

theta_2 = pi/2;
t2 = linspace(theta_2/w_0, 3*T, 500);

theta_1 = pi/4;
t1 = linspace(theta_1/w_0, 3*T, 500);

x1 = x_function(t1, x_0, w_0, theta_1);
x2 = x_function(t2, x_0, w_0, theta_2);

zero_crossing_fs_1 = 1e15*fsolve(@(t)x_function(t, x_0, w_0, theta_1), 1.4)

figure(2)
plot(t1*1e15,x1*1e9);
hold on
plot(t2*1e15,x2*1e9);
hold off
xlabel('t [fs]'), ylabel('x [nm]'), title('Position of electrons for given release phase'), legend('\theta = \pi/4','\theta = \pi/2');
