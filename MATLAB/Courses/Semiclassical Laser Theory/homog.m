num_steps = 8;
R0_vec = R0*ones(1,num_steps);
R0_vec(2:2:end) = 0;

t_fid = 100;

duration = td*ones(1,num_steps);
duration(2:2:end) = td*0.5;
