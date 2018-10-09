q = 1.6e-19;

lambda_midIR =  1.9;
E_increment = 0.1;

d_array = [80e-10];
    
materials = [100 0 .32; d_array(1)*1e10 0 0; 100 0 0.32];
use_infinite_well = 0;


for ind = 1:length(d_array)

if (use_infinite_well)
    [wv_c1, E_c1] = GaAs_infinite_well(d_array(ind), 1, 'c', 500);
    [wv_c2, E_c2] = GaAs_infinite_well(d_array(ind), 2, 'c', 500);
    [wv_lh1, E_lh1] = GaAs_infinite_well(d_array(ind), 1, 'lh', 500);
    [wv_lh2, E_lh2] = GaAs_infinite_well(d_array(ind), 2, 'lh', 500);
    [wv_hh1, E_hh1] = GaAs_infinite_well(d_array(ind), 1, 'hh', 500);
    [wv_hh2, E_hh2] = GaAs_infinite_well(d_array(ind), 2, 'hh', 500);
    
    electron_states.energies = 1000/q*[E_c1; E_c2] + 1424;
    electron_states.wavefunctions = [wv_c1; wv_c2];
    
    light_hole_states.energies = -1000/q*[E_lh1; E_lh2];
    light_hole_states.wavefunctions = [wv_lh1; wv_lh2];
    
    heavy_hole_states.energies = -1000/q*[E_hh1; E_hh2];
    heavy_hole_states.wavefunctions = [wv_hh1; wv_hh2];
    
else
    [heavy_hole_states, light_hole_states, electron_states] = AlGaAs_states(materials, 1, E_increment);
end

[alpha_ND, Ep_1, Ep_2, absorption_edge] = ND2PA_TM_finite_parabolic(lambda_midIR, d_array(ind), q/1000*(electron_states.energies - 1424), electron_states.wavefunctions, -q/1000*light_hole_states.energies, light_hole_states.wavefunctions, -q/1000*heavy_hole_states.energies, heavy_hole_states.wavefunctions);

figure(3)
% plot((Ep_1 + Ep_2)/absorption_edge, alpha_ND)
plot(1.24./(Ep_2/q), alpha_ND);
hold on
end
ylabel('\alpha^{ND}_2 (cm/GW)'), xlabel('\lambda_{probe}'), grid on
title(strcat('\alpha^{ND}_2 for \lambda = ', sprintf(' %1.2f ', lambda_midIR), '\mum'))

hold off
