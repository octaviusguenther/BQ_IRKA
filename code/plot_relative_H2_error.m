%create large Full order System (k=10) and corresponding initial reduced
%order systems (r=4,5,...)
%run BQ_IRKA_v3
%calculate for each initial system the absolute H2 - error between the full
%order models and the reduced order models
%Plot the error in a line chart
function plot_relative_H2_error(model_name,FOM_dim,MAX_ROM_dim)
k = FOM_dim;
Sigma = BQ_system(k,model_name);
n=Sigma.dim;

%specify reduced orders and BQ-IRKA options
reduced_order = 1:MAX_ROM_dim;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                ;
tol = 1e-12;
max_iter = 100;
Norm_formula ='P';

%pre calculate FOM reachability Gramian to be pass into error H2-norm
%function

FOM_P_Gram = gen_sylv(Sigma,Sigma);

%precompute for the relative error
h = getH2norm(Sigma,Norm_formula,FOM_P_Gram);

%initialize H2_error array
H2_error = zeros(size(reduced_order,2),2);


%iterator over reduced dimensions
for iter = 1:size(reduced_order,2)
    
    inital_guess = BQ_system(reduced_order(iter),'rand');
    %apply BQ-IRKA

        try 
        [Sigma_IRKA,~] = BQ_IRKA_v3(Sigma,reduced_order(iter),'Norm_formula',Norm_formula,'-tol',tol,'-max_iter',max_iter,'-initial_guess',inital_guess);
        
        catch error
            H2_error(iter,1) = NaN;
        
        end


        %calculate H2 - Norm
        H2_error(iter,1) = getErrorH2norm(Sigma,Sigma_IRKA,Norm_formula,FOM_P_Gram)/h;
         
end

%plot
x = reduced_order;
semilogy(x,(H2_error(:,1)),'-o');

%plot attributes
title_string = sprintf('%s-model dim=%d',model_name,n);
title(title_string);
legend('random initial guess','Location','northeast')
xlabel('reduced order','Interpreter','latex');
ylabel('relative $$\mathcal{H}_2$$ error $$\frac{\|\Sigma_{ROM} - \Sigma_{FOM}\|_{H_2}}{\|\Sigma_{FOM}\|_{H_2}}$$','Interpreter','latex');

%save as pdf
ax = gca;
export_string = sprintf('./../plots/%s_dim_%d_rel_error.pdf',model_name,n);
exportgraphics(ax,export_string,'Resolution',300);



