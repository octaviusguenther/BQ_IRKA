%solves the generalized Sylvester equation
%Sigma1.A*X + X*Sigma2.A' + Sigma1.N*X*Sigma2.N2' + Sigma1.b*Sigma2.b'=0
%by solving the equivalent vectorized linear system via MATLAB's
%'backslash' operator
function X = gen_sylv(Sigma_1,Sigma_2,method)
    if (~exist('method','var'))
        method = 'direct';
    end

    n = size(Sigma_1.A,1);
    r = size(Sigma_2.A,1);
    
    rhs = reshape(Sigma_1.b*Sigma_2.b',[n*r,1]);
    A= -(kron(speye(r),Sigma_1.A) + kron(Sigma_2.A,speye(n)) + kron(Sigma_2.N,Sigma_1.N));
    
    switch method
        case 'gmres'
            [vec_x,flag] = gmres(A,rhs,10,1e-16,10);
        case 'direct'
            vec_x = A\rhs;
    end
    
    X=sparse(reshape(vec_x,[n,r]));