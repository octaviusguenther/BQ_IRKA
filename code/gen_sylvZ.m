%solves Sigma1.A^T*Z + Z*Sigma2.A + Sigma1.N^T*Z*Sigma2.N + Sigma1.M*X*Sigma2.M = 0
%via solving the equivalent vectorized linear system either directly with
%'backslash'
%X is the solution to sylvester equation 
%Sigma1.A*X + X*Sigma2.A' + Sigma1.N*X*Sigma2.N2' + Sigma1.b*Sigma2.b'=0
function Z = gen_sylvZ(Sigma1,Sigma2,X,method)

    if (~exist('method','var'))
        method = 'direct';
    end

    n = size(Sigma1.A,1);
    r = size(Sigma2.A,1);
    
    %specify right-hand side 
    rhs = -1*kron(Sigma2.M,Sigma1.M) * reshape(X,[n*r,1]);
    A= (kron(speye(r),Sigma1.A') + kron(Sigma2.A',speye(n)) + kron(Sigma2.N',Sigma1.N'));

    switch method
        %gmres method only for experimantel purposes
        case 'gmres'
            [vec_z,flag] =gmres(A,rhs,10,1e-14,10);
        case  'direct'
            vec_z=A\rhs;
    end

    Z=reshape(vec_z,[n,r]);