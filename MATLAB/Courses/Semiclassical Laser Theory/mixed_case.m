rho_aa_0 = 1;
m = 1000;
rho_aa_mat = rho_aa_0*ones(1,m)
pop = zeros(1,100);
pop(1) = 1/m*sum(rho_aa_mat)

t_step = 0.01;

for ind = 2:length(pop)
    random_mat = round(0.52*rand(1,round(m*pop(ind-1))));
    rho_aa_mat(rho_aa_mat == 1) = xor(rho_aa_mat(rho_aa_mat == 1), random_mat);
    pop(ind) = 1/m*sum(rho_aa_mat);
end
plot(pop)