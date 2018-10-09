function [normalized_probe_energy_out, varargout] = xpm_1d_fcn(t_fs, phi_pump, ... 
    phi_probe, sample_thickness_mm, num_z_steps, gamma_ep_minv_Winv, gamma_ee_minv_Winv, ... 
    gvm_fs_per_mm, ...
    pump_dispersion_fs2_per_mm, probe_dispersion_fs2_per_mm, output_normalized_energy_only)

c = 2.998e8;

tmax_fs = max(t_fs);
gamma_ep_mminv_Winv = 1e-3*gamma_ep_minv_Winv;
gamma_ee_mminv_Winv = 1e-3*gamma_ee_minv_Winv;

num_time_points = length(t_fs);

% Fourier space for x and y
omega = 2*pi/tmax_fs*[(0:num_time_points/2-1) (-num_time_points/2:-1)];

% Propagation distance
dz_mm = sample_thickness_mm/num_z_steps;
dt_s = 1e-15*tmax_fs/num_time_points;
initial_probe_energy_W =  dt_s*sum(sum(abs(phi_probe(:,:)).^2));


% Linear propagation operator
D_pump = -(1i*pump_dispersion_fs2_per_mm*omega.^2 - 1i*(gvm_fs_per_mm)*omega)*0.5*dz_mm;
D_probe = -(1i*probe_dispersion_fs2_per_mm*omega.^2)*0.5*dz_mm;


if (~output_normalized_energy_only)
    phi_out_probe = zeros(num_z_steps, num_time_points);
    phi_out_pump = zeros(num_z_steps, num_time_points);
end

z_mm = dz_mm:dz_mm:sample_thickness_mm;

debug = 1;
% Not symmetrized, will need to do that for better accuracy
for ind = 1:num_z_steps
    
    %Pump
    lin_step_pump_1 = ifft(fft(phi_pump).*exp(D_pump));
    nonlin_step_pump = lin_step_pump_1.*exp((1i*gamma_ee_mminv_Winv*abs(lin_step_pump_1).^2)*dz_mm);
    lin_step_pump_2 = ifft(fft(nonlin_step_pump).*exp(D_pump));
    phi_pump = lin_step_pump_2;
    


    
    % Probe
    lin_step_probe_1 = ifft(fft(phi_probe).*exp(D_probe));
    nonlin_step_probe = lin_step_probe_1.*exp((2*1i*gamma_ep_mminv_Winv*abs(lin_step_pump_1).^2)*dz_mm);
    lin_step_probe_2 = ifft(fft(nonlin_step_probe).*exp(D_probe));
    phi_probe = lin_step_probe_2;
    
    if (debug == 1)
        pause(0.1)
        figure(1)
        plot(abs(phi_probe).^2);
        figure(2)
        plot(abs(fftshift(fft(lin_step_probe_2))).^2)
    end
    
    if (~output_normalized_energy_only)
        phi_out_pump(ind, :) = phi_pump;
        phi_out_probe(ind,:) = phi_probe;
    end
end

final_probe_energy_W =  dt_s*sum(sum(abs(phi_probe(:,:)).^2));
normalized_probe_energy_out = (initial_probe_energy_W - final_probe_energy_W)/initial_probe_energy_W;

if ~output_normalized_energy_only
    varargout{1} = phi_out_probe;
    varargout{2} = phi_out_pump;
    varargout{3} = z_mm;
end

end

