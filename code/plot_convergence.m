%plot the error evolution of the iteratives
function plot_convergence(model,ROM_dim,P)
%full-order system


Sigma= model;
n= Sigma.dim;
model_name = Sigma.name;
tol = 1e-10;

%initial guess for the algorithm
initial_guess1 = BQ_system(ROM_dim,'rand');
r1 = initial_guess1.dim;

%initial_guess2 = BQ_system(2,model);
%r2 = initial_guess2.dim;

%running the algorithm for only 10 iterations given a random initial guess
[Sigma_IRKA_randini,conv_evo_rand] = BQ_IRKA_v3(Sigma,r1,'-initial_guess',initial_guess1,'-error_evo',true,'-max_iter',20,'-P',P);

%examine the optimal interpolation points, which are given as
%-1*\sigma(A_opt)
p = optimal_points(Sigma_IRKA_randini);
p
%create a initial guess having the eigenvalues of the calculated model as
%poles and the rest is random
rng(42);
initial_guess2 = BQ_system(r1,0);
initial_guess2.A = diag(p);
initial_guess2.N = randi([-10 10],r1,r1);
initial_guess2.b = randi([-10 10],r1,1);
M = randi([-10 10],r1,r1);
                    
initial_guess2.M = M'*M;
initial_guess2.dim = size(initial_guess2.A,1);


%running the algorithm for 10 iteration given optimal initial guess
[Sigma_IRKA_optini,conv_evo_opt] = BQ_IRKA_v3(Sigma,r1,'-initial_guess',initial_guess2,'-error_evo',true,'-max_iter',20,'-tol',tol,'-P',P);

%running eigenvalue devations
[Sigma_IRKA,~,changein_errs] = BQ_IRKA_v3(Sigma,r1,'-initial_guess',initial_guess1,'-tol',tol,'-max_iter',1000);
[Sigma_IRKA_opt,~,changein_errs2] = BQ_IRKA_v3(Sigma,r1,'-initial_guess',initial_guess2,'-tol',tol,'-max_iter',1000);

%plotting the error of consecutives iterations
%define figure
fig1 = figure;
title_string=sprintf('convergence %s-model FOM-dim=%d ROM-dim=%d',model_name,n,r1);
sgtitle(title_string,'Interpreter','latex');
%define frist subplot
subplot(2,1,1);

%first subplot
x = 1:1:10;

semilogy(x,abs(conv_evo_rand(1:10)),'--om',x,abs(conv_evo_opt(1:10)),'-og');
yticks([10^-6,10^-3,10^0])

%title_s = sprintf('convergence history for model=%s, n=%d, r=%d',model,n,r1);
title_s = sprintf('Relative $H_2$-errors');
title(title_s,'Interpreter','latex');
legend('random initial guess','optimal interpolation points guess',Location='northeast')

xlabel('Number of iteration','Interpreter','latex');
ylabel('relative $$\mathcal{H}_2$$ error $$\frac{\|\Sigma_r^{(k)}-\Sigma_r^{(k-1)}\|_{H_2}}{\|\Sigma_r^{(k-1)}\|_{H_2}}$$','Interpreter','latex');

%define second subplot
subplot(2,1,2);
%plot 2

semilogy(changein_errs,'--m')
yticks([10^-10,10^0])
hold on
semilogy(changein_errs2,'-g')
%plot 2 attributes
title('Deviation of the eigenvalues','Interpreter','latex');
legend('random initial guess','optimal interpolation points guess','Location','northeast')
xlabel('Number of iterations','Interpreter','latex');
ylabel('$$\max|\mu_i - \mu_{i+1}|$$','interpreter','latex');


poles(Sigma_IRKA)
poles(Sigma_IRKA_opt)
%save as pdf
export_string = sprintf('./../plots/%s_FOM%d_ROM%d_convergence.pdf',model_name,n,r1);
exportgraphics(fig1,export_string,'Resolution',500);

