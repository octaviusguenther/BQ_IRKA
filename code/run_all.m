%runs all the experiments and plots the figures

n=12;
%heat model plots
%outputs 
plot_outputs('heat',n,8);

%relative erroes
plot_relative_H2_error('heat',n,12);

%convergence properties
plot_convergence('heat',n,5);

%flow model plots
%outputs
plot_outputs('flow',n,8);

%relative errors
plot_relative_H2_error('flow',n,12);

%convergence properties
plot_convergence('flow',n,5);

