%runs all the experiments and plots the figures

n = 20;
r = 10;
max_rom_dim = 12;

Sigma_heat = BQ_system(n,'heat');
Sigma_flow = BQ_system(n,'flow');


tic
%heat model plots
%outputs 
plot_outputs(Sigma_heat,r);


%precompute heat model reachability Gramian P_heat
fprintf(1, 'Precomputing full-order reachability Gramian P \n');
fprintf(1, '-------------------------------------------\n');
P_heat = gen_sylv_naive(Sigma_heat,Sigma_heat,100,1e-10);

%relative erroes
plot_relative_H2_error(Sigma_heat,max_rom_dim,P_heat);

%convergence properties
plot_convergence(Sigma_heat,r,P_heat);

%flow model plots
%outputs
plot_outputs(Sigma_flow,r);


%precompute flow model reachability Gramian P_flow
fprintf(1, 'Precomputing full-order reachability Gramian P \n');
fprintf(1, '-------------------------------------------\n');
P_flow = gen_sylv_naive(Sigma_flow,Sigma_flow,100,1e-10);
%relative errors
plot_relative_H2_error(Sigma_flow,max_rom_dim,P_flow);



%convergence properties
plot_convergence(Sigma_flow,r,P_flow);
toc