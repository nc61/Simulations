L = 50;
N = 2^10;
T = (L/N)*(0:N-1);
w = 2*pi/L*[(0:N/2-1) (-N/2:-1)];

Z = 20;
steps = 1000;
dZ = Z/steps;
nl = 1;
alpha = 0;


u = sech(T-L/2);
v = 1i*u*dZ/2;

D = -0.5*(1i*w.^2 + alpha)*0.5*dZ;

u_out = zeros(steps, N);
v_out = zeros(steps, N);

for ind = 1:steps
    lin_step_1 = ifft(fft(u).*exp(D));
    nonlinear_step = lin_step_1.*exp(1i.*nl*(abs(lin_step_1).^2 + abs(v).^2 + v./u)*dZ);
    lin_step_2 = ifft(fft(nonlinear_step).*exp(D)); 
    u = lin_step_2;
    
    lin_step_1 = ifft(fft(v).*exp(D));
    nonlinear_step = lin_step_1.*exp(1i.*nl*(abs(lin_step_1).^2 + abs(u).^2 + u./v)*dZ);
    lin_step_2 = ifft(fft(nonlinear_step).*exp(D)); 
    v = lin_step_2;
  
    
    u_out(ind, :) = u;
    v_out(ind, :) = v;
end

z = linspace(0,Z,steps);
figure(1)
surf(T,z,abs(u_out).^2)
shading interp
view(2)
colormap('jet')
colorbar
xlabel('T (a.u.)')
ylabel('Distance (a.u.)')
xlim([min(T) max(T)])

figure(2)
surf(T,z,abs(v_out).^2)
shading interp
view(2)
colormap('jet')
colorbar
xlabel('T (a.u.)')
ylabel('Distance (a.u.)')
xlim([min(T) max(T)])

 
