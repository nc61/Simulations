%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% T1 Measurement sequence %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 5e-6;

duration =[2*tpi2 T tpi2 T/2];
R0_vec = R0*ones(1,length(duration));
R0_vec(2:2:end) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% T2 Measurement sequence (CPMG) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 3e-6;

duration =[tpi2 T 2*tpi2 2*T 2*tpi2 2*T 2*tpi2 2*T 2*tpi2 2*T 2*tpi2 2*T 2*tpi2 2*T 2*tpi2 2*T 2*tpi2 2*T 2*tpi2];
R0_vec = R0*ones(1,length(duration));
R0_vec(2:2:end) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stimulated Echo sequence %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 5e-6;
T = 12e-6;

R0_vec = R0*[1 0 1 0 1 0];
duration = [tpi2 t tpi2 T tpi2 2*T];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%