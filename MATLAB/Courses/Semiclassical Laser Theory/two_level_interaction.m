
R0 = 1e8;
tpi2 = pi/(2*R0);


%%% PULSE SEQUENCE %%%
%%%%%%%%%%%%%%%%%%%%%%

T = 3e-6;

duration =[tpi2 T 2*tpi2 2*T 2*tpi2 2*T 2*tpi2 2*T 2*tpi2 2*T 2*tpi2 2*T 2*tpi2 2*T 2*tpi2 2*T 2*tpi2 2*T 2*tpi2];
R0_vec = R0*ones(1,length(duration));
R0_vec(2:2:end) = 0;
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%


% get an expression for the total number of "blocks" (pulse or FID)
num_blocks = length(R0_vec);



% Define relaxation times for level a, b, and dephasing
gamma_aa = 10;
gamma_bb = 10;
gamma_ab = 1;

% Define incoherent pump rates.
lambda_a = 0;
lambda_b = 1;

% w_in is the inhomogeneous width
w_in = 10e6;

% Define a range of detunings to apply inhomogeneous broadening. The center
% of delta is the detuning from the center (usually 0 for resonance
% excitement)
delta = 3*w_in*linspace(-1,1,200);
Ww = 1/(w_in*sqrt(pi))*exp(-delta.^2./w_in^2);

% Combine lambda and gamma into matrices for calling the equation of motion
% function
lambda_mat = [lambda_a lambda_b];
gamma_mat = [gamma_aa, gamma_bb, gamma_ab];


% Define the initial conditions so that all the poopulation is in state b
% (lower state). prev_rho is used because the initial conditions for each
% step will be the final rho of the previous
initial_conditions = zeros(1,3*length(delta));
initial_conditions(2:3:end) = 1;
prev_rho = initial_conditions;

% Define a cell array to store
t_cellarray = cell(num_blocks,1);
rho_cellarray = cell(num_blocks,1);
R_cellarray = cell(num_blocks,1);

% Start the pulse sequence at t = 0
start = 0;

% Iterate through each block and solve the density matrix differential
% equations, using the previous rho values as initial conditions. Do this
% for each detuning (delta) at the same time.
for ind = 1:num_blocks
    % The Runge Kutta solver can give variable time steps, so it's
    % important to store these as well as rho values
    [t_block,rho_block] = pulse_block(start, duration(ind), R0_vec(ind), delta, prev_rho, lambda_mat, gamma_mat);
    
    % Make sure the next block starts where this one left off. Set ICs as
    % mentioned before
    start = t_block(end);
    prev_rho = rho_block(end,:);
    
    % Store the block in the cell array
    t_cellarray{ind} = t_block;
    rho_cellarray{ind} = rho_block;
    R_cellarray{ind} = R0_vec(ind)*ones(size(t_cellarray{ind}));
end

% Group everything together into one matrix
t_total = cell2mat(t_cellarray);
rho_total = cell2mat(rho_cellarray);
R_total = cell2mat(R_cellarray);
R_total = R_total/max(R_total);

% Replicate the inhomogeneous spectrum to elementwise multiply at every
% time
Ww = repmat(Ww, length(t_total), 1);

% Pick out each density matrix component from the total matrix
rho_aa_total = rho_total(:,1:3:end);
rho_bb_total = rho_total(:,2:3:end);
rho_ab_total = rho_total(:,3:3:end);

% Integrated rho values (not used)
rho_aa_integrated = trapz(Ww.*rho_aa_total, 2)/length(delta);
rho_ab_integrated = trapz(Ww.*imag(rho_ab_total), 2)/length(delta);

% Option to plot a single response at a given detuning (detuning_index)
plot_single = 0;

if (plot_single) == 1
    
    % Choose which detuning index to look at
    detuning_index = 1;
    figure(1)
    subplot(3,1,1)
    plot(t_total, rho_aa_total(:,detuning_index)), ylabel('\rho_aa');
    subplot(3,1,2);
    plot(t_total, rho_bb_total(:,detuning_index));
    subplot(3,1,3)
    plot(t_total, imag(rho_ab_total(:,detuning_index)));
end

% Compute U,V,W from the population matrix
U = rho_ab_total + conj(rho_ab_total);
V = 1i*(rho_ab_total - conj(rho_ab_total));
W = rho_aa_total - rho_bb_total;

plot_individuals_together = 0;

if plot_individuals_together == 1
    figure(2)
    points_to_plot = 10;
    
    skip_indices = ceil(size(U,2)/points_to_plot);
    downsample_indices = 1:skip_indices:size(U,2);
    
    t_plot = t_total(downsample_indices);
    U_plot = U(:,downsample_indices);
    V_plot = V(:,downsample_indices);
    W_plot = W(:,downsample_indices);
    
    % Plot the (non-integrated) values separately and ignore that they have
    % different magnitudes
    subplot(3,1,1)
    plot(t_total, U_plot), ylabel('U');
    subplot(3,1,2)
    plot(t_total, V_plot), ylabel('V');
    subplot(3,1,3)
    plot(t_total, W_plot), ylabel('W');
end

plot_integrated = 1;

if plot_integrated == 1
    if length(delta) > 1
        
        d_delta = delta(2) - delta(1);
        W_integrated = trapz(Ww.*W,2)*d_delta;
        V_integrated = trapz(Ww.*V,2)*d_delta;
        U_integrated = trapz(Ww.*U,2)*d_delta;
        
    else
        % Only one detuning, no integration
        U_integrated = U;
        V_integrated = V;
        W_integrated = W;
    end
    
    figure(3)
    
    % Plot the Bloch vector components
    subplot(3,1,1)
    t_us = t_total*1e6;
    plot(t_us, U_integrated), ylim([-1 1]), ylabel('U')
    hold on
    plot(t_us, R_total, 'r')
    hold off
    subplot(3,1,2)
    plot(t_us, V_integrated), ylabel('V')
    hold on
    plot(t_us, R_total, 'r')
    hold off
    subplot(3,1,3)
    plot(t_us, W_integrated), ylabel('W'), xlabel('t (\mus)');
    hold on
    plot(t_us, R_total, 'r')
    hold off
    
    figure(4)
    plot(t_us, V_integrated), ylabel('V'), xlabel('t (\mus)');
    hold on
    plot(t_us, R_total, 'r')
    hold off
    
end

