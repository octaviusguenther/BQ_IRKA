%solve the generalized Sylvester equation
%Sigma1.A*X + X*Sigma2.A' + Sigma1.N*X*Sigma2.N' + Sigma1.b*Sigma2.b' = 0
%iteratively via the method presented in Lemma 3.5
%the iteration stops if the maximal number of iterations 'maxiter' is
%reached or if the relative norm of two consecutive iteratives is smaller than 'tol'
function X = gen_sylv_naive(Sigma1,Sigma2,maxiter,tol)

%specifiy dimensions
n = size(Sigma1.A,1);
r = size(Sigma2.A,2);

%inital guess 
X=eye(n,r);
start_iter = tic;
%1. solving sylvester equations for X=lim Xk
    for k = 1:maxiter
        %specify rhs : -NX_(j-1)N^T -bb_r^T
        rhs1 =-Sigma1.N*X*(Sigma2.N') - Sigma1.b*Sigma2.b';
        % matlab sylvester solves AX + XB = C for X
        X_new = sylvester(full(Sigma1.A),full(Sigma2.A'),full(rhs1));
        
        rel_error=norm(X_new - X)/norm(X);
        fprintf(1,'rel. Error norm(X_(k+1) - X_(k))\\norm(X_(k))=%d\n',rel_error);
        fprintf(1, '-------------------------------------------\n');
        if (norm(X_new - X)/norm(X)) < tol
        fprintf(1,'Sylvester iteration converged in %.2f sec\n',toc(start_iter));
            break;
        end

        X = X_new;
    end