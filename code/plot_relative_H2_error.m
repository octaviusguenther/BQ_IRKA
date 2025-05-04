%Create either the 'heat' or the 'flow' 
%full-order System with dimensions FOM_dim and initial random reduced
%order systems with dimension r=1,...,MAX_ROM_dim).

%Run BQ_IRKA_v3 to compute all the reduced-system systems ranging from
%dimension r=1,...,MAX_ROM_dim.

%Calculate for each resulting system the relative H2 - error between the full
%order models and the reduced order models.

%Plot the error in a line chart
function plot_relative_H2_error(model,MAX_ROM_dim,P) 

Sigma = model;
n=Sigma.dim;
model_name = Sigma.name;
%specify reduced orders and BQ-IRKA options
reduced_order = 1:MAX_ROM_dim;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                ;
tol = 1e-10;
max_iter = 100;
Norm_formula ='P';

%pre calculate FOM reachability Gramian to be pass into error H2-norm
%function

FOM_P_Gram = P;

%precompute the H_2 norm of the full-order model Sigma
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
        iter_rel_error = tic;
        H2_error(iter,1) = getErrorH2norm(Sigma,Sigma_IRKA,Norm_formula,FOM_P_Gram)/h;
        fprintf(1, 'Rel.Error computed in %.2f s\n',toc(iter_rel_error));
        fprintf(1, '---------------------------------------\n')
end

%plot
x = reduced_order;
semilogy(x,(H2_error(:,1)),'-o');

%plot attributes
title_string = sprintf('%s-model dim=%d',model_name,n);
title(title_string,'Interpreter','latex');
legend('BQ-IRKA','Location','northeast')
xlabel('reduced order','Interpreter','latex');
ylabel('relative $$\mathcal{H}_2$$ error $$\frac{\|\Sigma_{ROM} - \Sigma_{FOM}\|_{H_2}}{\|\Sigma_{FOM}\|_{H_2}}$$','Interpreter','latex');

%save as pdf
ax = gca;
export_string = sprintf('./../plots/%s_dim_%d_rel_error.pdf',model_name,n);
exportgraphics(ax,export_string,'Resolution',300);



