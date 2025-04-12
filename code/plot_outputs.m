function plot_outputs(modelname,full_order,reduced_order)

% initiate Full ordermodel
FOM = BQ_system(full_order,modelname);
n = FOM.dim;

%specify input u
u = @(t) cos(2*pi*t);

%specify grid
t = linspace(0,2,2e+4);
%solve FOM
[t,Y] = get_output(FOM,u,t);

%calculate ROM by applying BQ_IRKA
%TO DO specify tolerance
%TO DO specify max_iterations
ROM = BQ_IRKA_v3(FOM,reduced_order,'-max_iter',500,'-tol',1e-14);
r = ROM.dim;

%solve ROM
[~,Y_r] = get_output(ROM,u,t);

%L_inf norm of error y - y_r
y_error = abs(Y - Y_r);

%plots
%fiure 1: outputs y(k) and y_r(k)
fig1=figure;
plot(t,Y,'-b',t,Y_r,'--r');

%plot properties
title_string = sprintf('Observed outputs, %s-model, FOM-dim=%d, ROM-dim=%d',modelname,n,r);
title(title_string,'Interpreter','latex');

legend('Full-order system','Reduced-order system','Location','northeast');
xlabel('time $t$','Interpreter','latex');
ylabel('outputs $y(t)$','Interpreter','latex');

%save as pdf
export_string1 = sprintf('./../plots/%s_outputs_FOM_dim_%d_ROM_dim_%d.pdf',modelname,n,r);
exportgraphics(fig1,export_string1,'Resolution',300);

%figure2: L_inf error between y and y_r
fig2 = figure;

%plot y_error
darkgreen = [0 0.5 0];
semilogy(t,y_error,'-','Color',darkgreen);
yticks([1e-8,1e-5,1e-2,1e-1]);

%plot properties

title_string = sprintf('Magnitutde of output error, %s-model, FOM-dim=%d, ROM-dim=%d',modelname,n,r);
title(title_string,'Interpreter','latex');
xlabel('time $t$','Interpreter','latex');
ylabel('$|y(t)-y_r(t)|$','Interpreter','latex');

%save as pdf
export_string2 = sprintf('./../plots/%s_y_error_FOM_dim_%d_ROM_dim_%d.pdf',modelname,n,r);
exportgraphics(fig2,export_string2,'Resolution',300);







