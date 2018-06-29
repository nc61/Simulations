
% Number of fourier modes
L_space = 100;
L_time = 66;

% Number of points to split up x into
N_space = 2^8;
N_time = 2^8;

% points for x and y 
x = (L_space/N_space)*(0:N_space-1);
t = (L_time/N_time)*(0:N_time-1);
[X, T] = meshgrid(x,t);

% Fourier space for x and y
kx_pump = 2*pi/L_space*[(0:N_space/2-1) (-N_space/2:-1)];
omega_pump = 2*pi/L_time*[(0:N_time/2-1) (-N_time/2:-1)];


[Kx_pump, OMEGA_pump] = meshgrid(kx_pump,omega_pump);

% nonlinear
Aeff_mm_squared = 20e-6;

alpha_2_d_mm_per_W = 50e-8;

%alpha_2_d_mm_per_W = 0;

gamma2_d = alpha_2_d_mm_per_W/(2*Aeff_mm_squared);


% Propagation distance
Z = 5;
steps = 20;
dZ = Z/steps;

% Initial envelope
phi_pump = exp(-1i*X).*exp(-((X - 0.5*L_space).^2./5^2 + (T - 0.5*L_time).^2./6.617.^2));

% Linear propagation operator
D_pump = -(0.5*(-1i*5*Kx_pump.^2 + 1i*5*OMEGA_pump.^2)*0.5*dZ);

phi_out_pump = zeros(steps, N_time, N_space);


% Not symmetrized, will need to do that for better accuracy
for ind = 1:steps
    lin_step_pump_1 = ifft2(fft2(phi_pump).*exp(D_pump));

    nonlin_step_pump = lin_step_pump_1.*exp((-2*gamma2_d*abs(lin_step_pump_1).^2)*dZ);
    
    lin_step_pump_2 = ifft2(fft2(nonlin_step_pump).*exp(D_pump));
    
    phi_pump = lin_step_pump_2;

    phi_out_pump(ind, :, :) = phi_pump;
end

P = zeros(1,steps);
%Power
dX = L_space/N_space;
dT = L_time/N_time;

for ind = 1:steps
    P(ind) = dX*dT*sum(sum(abs(phi_out_pump(ind,:,:)).^2));
end
 
z = linspace(0,Z,steps);
figure(1)
% surf(x,z,abs(phi_out(:,:,N_time/2)).^2)
surf(t,z,abs(squeeze(phi_out_pump(:,:,N_space/2+1))).^2)
shading interp, colormap('jet')

figure(2)
surf(x,z,abs(squeeze(phi_out_pump(:,N_time/2+1,:))).^2)
shading interp
colormap('jet')
colorbar
xlabel('X (a.u.)')
ylabel('Distance (a.u.)')
%xlim([min(x) max(x)])



%Power conservation
figure(3)
plot(z,P,'Linewidth',2)
xlabel('Distance (a.u.)')
ylabel('P')
%ylim([P(1)-1 P(1)+1])


show_video = 0;
if show_video, for ind = 1:steps
        imshow(squeeze(abs(phi_out_pump(ind,:,:)).^2)./max(max(abs(phi_out_pump(ind,:,:)).^2)))
        pause(0.01)
        end
end

 
