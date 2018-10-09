function [heavy_hole_states, light_hole_states, electron_states] = AlGaAs_states(materials, plot_wells, Eincrement)
% For a given matrix of materials, calculate the bound states in the
% conduction and valence bands.
% 
%   calls: schrodinger.m, solver for schrodinger equation
%   
%   materials: each row of the materials matrix is a different material
%   with thickness in Angstroms. The first column is the composition (Ex:
%   0.32) and the third column is the thickness. The second column is
%   reserved for the future when this supports more material systems
%   plot_wells: Display the wavefunctions for visualization
%   Eincrement: Passed through to schrodinger.m, energy resolution of
%   schrodinger solver (very important)

% Define constants
m0 = 9.109e-31;

% Go through the materials matrix and find the composition at every
% individual Angstrom step in z
x = [];
for ind = 1:size(materials, 1)
    x = [x materials(ind,3)*ones(1, materials(ind, 1))];  %#ok<AGROW>
end

% Bandgaps (meV)
T = 295;
Eg_GaAs = 1000*(1.519 - 5.405e-4*T^2/(T + 204));

% Effective masses
% Adachi properties of Aluminum Gallium Arsenide
% m_e_eff_AlGaAs_z = (0.0665 + 0.083.*x).*m0;
% m_hh_eff_AlGaAs_z = -(0.333 + 0.18*x).*m0;
% m_lh_eff_AlGaAs_z = -(0.090 + 0.09*x).*m0;

% m_e_eff_AlGaAs = (0.0665 + 0.083.*x).*m0;
% m_hh_eff_AlGaAs = -(0.333 + 0.42*x).*m0;
% m_lh_eff_AlGaAs = -(0.094 + 0.043*x).*m0;

m_e_eff_AlGaAs_z = (0.0665 + 0.083.*x).*m0;
m_hh_eff_AlGaAs_z = -(0.333 + 0.084*x).*m0;
m_lh_eff_AlGaAs_z = -(0.090 + 0.09*x).*m0;

% V (conduction band offset)
Ec_offset = 790.*x;
Ev_offset = -460.*x;

% Find the total thickness of the structure
z = 1e-10*(1:sum(materials(:,1)));

% Applied voltage, should be added as a parameter and not automatically set
% to 0
V_app= 0;

% Potential barriers with adjustments based on applied bias.
V_c = Ec_offset - V_app*z/max(z);
V_hh = Ev_offset - V_app*z/max(z);
V_lh = Ev_offset - V_app*z/max(z);

% Use the schrodinger solver to solve for electron, light hold and heavy
% hole states
[wavefunctions_hh, energies_hh, dz_hh] = schrodinger(V_hh, z, m_hh_eff_AlGaAs_z, Eincrement);
[wavefunctions_lh, energies_lh, dz_lh] = schrodinger(V_lh, z, m_lh_eff_AlGaAs_z, Eincrement);
[wavefunctions_c, energies_c, dz_c] = schrodinger(V_c, z, m_e_eff_AlGaAs_z, Eincrement);

% Shift the conduction band energies upward by the bandgap. Now the top of
% the valence bands are at 0 energy and the bottom of the conduction band
% is at 1424meV (GaAs bandgap). THIS SHOULD BE CHANGED IF THE WELLS AREN'T
% MADE OUT OF GAAS
energies_c = energies_c + Eg_GaAs;
V_c = V_c + Eg_GaAs;

% Save the heavy hole states into a structure for output
heavy_hole_states.wavefunctions = wavefunctions_hh;
heavy_hole_states.energies = energies_hh;
heavy_hole_states.dz = dz_hh;

% Save the light hole states into a structure for output
light_hole_states.wavefunctions = wavefunctions_lh;
light_hole_states.energies = energies_lh;
light_hole_states.dz = dz_lh;

% Save the electron states into a structure for output
electron_states.wavefunctions = wavefunctions_c;
electron_states.energies = energies_c;
electron_states.dz = dz_c;

% Plot the output if desired
if ~isempty(energies_hh) && (plot_wells == 1) 
    
    % Kind of a workaround to make sure we don't cover up this plot with
    % whavever we are plotting afterward
    figure(5)

    scale = .5*energies_hh(1)/max(wavefunctions_hh(1,:));
    for ind = 1:size(wavefunctions_hh, 1)
        %diagram = ones(size(z)).*energies_hh(ind);
        diagram = energies_hh(ind) + scale*wavefunctions_hh(ind, :);
        plot(z*1e9, diagram, 'b');
        hold on;
    end
    
    scale = .5*energies_lh(1)/max(wavefunctions_lh(1,:));
    for ind = 1:size(wavefunctions_lh, 1)
        %diagram = ones(size(z)).*energies_lh(ind);
        diagram = energies_lh(ind) + scale*wavefunctions_lh(ind, :);
        plot(z*1e9, diagram, 'r');
        hold on;
    end
    
    scale = .5*((energies_c(1) - Eg_GaAs)/max(wavefunctions_c(1,:)));
    for ind = 1:size(wavefunctions_c, 1)
        %diagram = ones(size(z)).*energies_c(ind);
        diagram = energies_c(ind) + scale*wavefunctions_c(ind, :);
        plot(z*1e9, diagram, 'g');
        hold on;
    end
    
    plot(z*1e9, V_c, 'k--');
    title('Solutions to Schrodinger equation for GaAs/AlGaAs structures'), xlabel('z (nm)'), ylabel('E (meV)');
    plot(z*1e9, V_hh, 'k--');
    hold off;
end


end

