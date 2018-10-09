R0 = 1e8;
tpi2 = pi/(2*R0);

T = 13.91e-6;

duration =[2*tpi2 T tpi2 T/2];
R0_vec = R0*ones(1,length(duration));
R0_vec(2:2:end) = 0;

% get an expression for the total number of "blocks" (pulse or FID)
num_blocks = length(R0_vec);

% Define relaxation times for level a, b, and dephasing
gamma_aa = 5e4;
gamma_bb = 5e4;
gamma_ab = 0e4;

% Define incoherent pump rates.
lambda_a = 0;
lambda_b = 1;

% w_in is the inhomogeneous width
w_in = 5e6;