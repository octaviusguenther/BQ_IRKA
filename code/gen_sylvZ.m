function Z = gen_sylvZ(Sigma1,Sigma2,X,method)
%solves Sigma1.A^T*Z + Z*A2 + Sigma1.N^T*Z*N2 + MXMr = 0
    if (~exist('method','var'))
        method = 'direct';
    end

    n = size(Sigma1.A,1);
    r = size(Sigma2.A,1);
    
    rhs = -1*kron(Sigma2.M,Sigma1.M) * reshape(X,[n*r,1]);
    A= (kron(speye(r),Sigma1.A') + kron(Sigma2.A',speye(n)) + kron(Sigma2.N',Sigma1.N'));

    switch method
        case 'gmres'
            [vec_z,flag] =gmres(A,rhs,10,1e-14,10);
        case  'direct'
            vec_z=A\rhs;
    end

    Z=reshape(vec_z,[n,r]);