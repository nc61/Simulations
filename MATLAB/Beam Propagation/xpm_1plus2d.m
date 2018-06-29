
% Number of fourier modes
L_space = 100;
L_time = 100;

% Number of points to split up x into
N_space = 2^5;
N_time = 2^5;

% points for x and y
x = (L_space/N_space)*(0:N_space-1);
y = (L_space/N_space)*(0:N_space-1);
t = (L_time/N_time)*(0:N_time-1);
[X, Y, T] = ndgrid(x,y,t);

% Fourier space for x and y
kx_pump = 2*pi/L_space*[(0:N_space/2-1) (-N_space/2:-1)];
ky_pump = 2*pi/L_space*[(0:N_space/2-1) (-N_space/2:-1)];
omega_pump = 2*pi/L_time*[(0:N_time/2-1) (-N_time/2:-1)];


[Kx, Ky, OMEGA] = ndgrid(kx_pump, ky_pump, omega_pump);

% nonlinear
Aeff_mm_squared = 20e-6;

alpha_2_d_mm_per_W = 250e-8;
alpha_2_mm_per_W = 200e-8;

%alpha_2_d_mm_per_W = 0;

gamma2 = alpha_2_mm_per_W/(2*Aeff_mm_squared);
gamma2_d = alpha_2_d_mm_per_W/(2*Aeff_mm_squared);


% Propagation distance
Z = 0.5;
steps = 10;
dZ = Z/steps;

delays = linspace(0.3, 0.7, 100);

signal = zeros(size(delays));

for jnd = 1:length(delays)
    
    % Initial envelope
    phi_pump = exp(-((X - 0.5*L_space).^2./5^2 + (Y - 0.5*L_space).^2./5^2 + (T - delays(jnd)*L_time).^2./6.617.^2));
    phi_probe = exp(-((X - 0.5*L_space).^2./5^2 + (Y - 0.5*L_space).^2./5^2 + (T - 0.5*L_time).^2./6.617.^2));
    
    % Linear propagation operator
    D_pump = -(0.5*(-1i*0*Kx.^2 - 1i*0*Ky.^2 - 20*1i*Kx + 1i*0*OMEGA.^2)*0.5*dZ);
    D_probe = -(0.5*(-1i*0*Kx.^2 - 1i*0*Ky.^2 + 1i*0*OMEGA.^2)*0.5*dZ);
    
    phi_out_pump = zeros(steps, N_time, N_space, N_space);
    phi_out_probe = zeros(steps, N_time, N_space, N_space);
    
    % Not symmetrized, will need to do that for better accuracy
    for ind = 1:steps
        
        %Pump
        
        lin_step_pump_1 = ifft2(fft2(phi_pump).*exp(D_pump));
        
        nonlin_step_pump = lin_step_pump_1.*exp((-2*gamma2_d*abs(lin_step_pump_1).^2)*dZ);
        
        lin_step_pump_2 = ifft2(fft2(nonlin_step_pump).*exp(D_pump));
        
        phi_pump = lin_step_pump_2;
        
        % Probe
        
        lin_step_probe_1 = ifft2(fft2(phi_probe).*exp(D_probe));
        
        nonlin_step_probe = lin_step_probe_1.*exp((-2*gamma2*abs(lin_step_pump_1).^2)*dZ);
        
        lin_step_probe_2 = ifft2(fft2(nonlin_step_probe).*exp(D_probe));
        
        phi_probe = lin_step_probe_2;
        
        phi_out_pump(ind, :, :, :) = phi_pump;
        phi_out_probe(ind, :, :, :) = phi_probe;
    end
    
    P_pump = zeros(1,steps);
    P_probe = zeros(1,steps);
    %Power
    dX = L_space/N_space;
    dY = L_space/N_space;
    dT = L_time/N_time;
    
    for ind = 1:steps
        P_pump(ind) = dX*dY*dT*sum(sum(sum(abs(phi_out_pump(ind,:,:,:)).^2)));
        P_probe(ind) = dX*dY*dT*sum(sum(sum(abs(phi_out_probe(ind,:,:,:)).^2)));
    end
    
    signal(jnd) = P_probe(end);
end

z = linspace(0,Z,steps);
figure(1)
% surf(x,z,abs(phi_out(:,:,N_time/2)).^2)
surf(t,z,abs(squeeze(phi_out_pump(:,:,N_space/2+1,N_space/2+1))).^2)
hold on
surf(t,z,abs(squeeze(phi_out_probe(:,:,N_space/2+1,N_space/2+1))).^2)
hold off
shading interp, colormap('jet')

figure(2)
surf(x,z,abs(squeeze(phi_out_pump(:,N_time/2+1,:,N_space/2+1))).^2)
shading interp
colormap('jet')
colorbar
xlabel('X (a.u.)')
ylabel('Distance (a.u.)')
%xlim([min(x) max(x)])


%Power conservation
figure(3)
plot(z,P_pump,'Linewidth',2)
xlabel('Distance (a.u.)')
ylabel('P')
%ylim([P(1)-1 P(1)+1])

figure(4)
plot(z,P_probe,'Linewidth',2)
xlabel('Distance (a.u.)')
ylabel('P')
%ylim([P(1)-1 P(1)+1])

figure(5)
plot(signal)

show_video = 0;
if show_video, for ind = 1:steps
        imshow(squeeze(abs(phi_out_pump(ind,:,:)).^2)./max(max(abs(phi_out_pump(ind,:,:)).^2)))
        pause(0.01)
    end
end

max(signal) - min(signal)
