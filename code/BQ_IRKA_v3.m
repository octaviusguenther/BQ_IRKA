function [Sigma_opt,H2_error_iteratives,changein_errs] = BQ_IRKA_v3(Sigma,r,varargin)


%Defaults
error_evo = false;
max_iter = 100;
tol = 1e-14;
Norm_formula = 'P';
Sigma_r = BQ_system(r,'rand');

%precomputing full order reachability Gramian
%P = gen_sylv_direct(Sigma,Sigma);

%optionals
for ii = 1:2:nargin-2
       
        
        if strcmp('-error_evo', varargin{ii})
            error_evo = varargin{ii+1};
        elseif strcmp('-Norm_formula', varargin{ii})
            Norm_formula = varargin{ii+1};
        elseif strcmp('-tol', varargin{ii})
            tol = varargin{ii+1};
        elseif strcmp('-max_iter',varargin{ii})
            max_iter = varargin{ii+1};
        elseif strcmp('-initial_guess',varargin{ii})
            Sigma_r = varargin{ii+1};
        end
end

% setting n and r
n = size(Sigma.A,1);
r = Sigma_r.dim;

%initialize array storing the errors |Sigma_(k-1) - Sigma_k|
H2_error_iteratives = zeros(max_iter,1);

%precompute the H2-norm of full-order model for the relative error
if error_evo == true

    h = getH2norm(Sigma,'P');
end


%defining Sigma_0
Ar = Sigma_r.A; Nr = Sigma_r.N; br = Sigma_r.b; Mr = Sigma_r.M;

%intialize iterator, used for c
iter = 1;

%compute the H2-error between the inital guess and the FOM
if error_evo == true
    H2_error_iteratives(iter)= getErrorH2norm(Sigma,Sigma_r,Norm_formula);
end

%set initial changeinerr to get into the while loop
changein_errs(iter)= 1 + tol;

while (changein_errs(iter) > tol && iter < max_iter)
    
    %start clock
    iter_start = tic;

    fprintf(1, 'Current iterate is k = %d\n', iter)
    fprintf(1, '---------------------------------------\n')

    %update matrix of eigenvalues of (prior) Sigma_k
    %used for convergence criterion
    [~,D] = eig(full(Ar));


    %Solve the sylvester equation V(-D) - AV - NVN_hat = b*b_hat'
    
    rhs1 = kron(br,Sigma.b);
    V = (kron(-Ar,speye(n)) - kron(speye(r),Sigma.A) - kron(Nr,Sigma.N));
    vect_V = V\rhs1;
    
    V = reshape(vect_V,[n,r]);

    %Solve the sylvester equation 
    % W*(-D) - Ar^T*W - Nr^T*W*N_hat = - Mr*V*M

    rhs2 = kron(Mr',Sigma.M)*vect_V;
    W = (kron(Ar',speye(n)) + kron(speye(r),Sigma.A') + kron(Nr',Sigma.N'));
    
    vect_W = W\rhs2; 
    W = reshape(vect_W,[n,r]);
    
    %orthogonalization
    [V,~] = qr(V,"econ");
    [W,~] = qr(W,"econ");


    %Petrov Galerkin Projection
    olg = W'*V;

    Ar = olg\W'*Sigma.A*V;
    Nr = olg\W'*Sigma.N*V;
    br= olg\W'*Sigma.b;
    Mr = V'*Sigma.M*V;
    
    
    %define the (posterior) Sigma_k
    Sigma_temp = BQ_system(r,0);
    Sigma_temp.A = Ar; Sigma_temp.N = Nr; Sigma_temp.b = br; Sigma_temp.M = Mr;

    

    % Calculating the new mirrorred eigenvalues of (posterior ) Sigma_k
    D_new = eig(full(Ar));
 

    %increase iterator
    iter = iter +1;
    
    if error_evo == true
        % calculate the relative error of |Sigma_(k-1) - Sigma_k|
        h = getH2norm(Sigma_r);
        dh = getErrorH2norm(Sigma_r,Sigma_temp);
        H2_error_iteratives(iter) = dh/h;
        
        %update the prior Sigma_(k-1)
        Sigma_r.A = Ar; Sigma_r.N = Nr; Sigma_r.b = br; Sigma_r.M = Mr;
    end

    
    fprintf(1, 'Computing maximal deviations of the poles from current to previous iteration \n')
    fprintf(1, '-------------------------------------------\n')
    
    %calculate convergence criterion
    %which is max_i |u_i^(k+1) - u_i^(k)| and where the eigenvalues are
    %sorted according to their magnitude
    changein_errs(iter) = max(abs(sort(eig(D)) - sort(D_new)))/max(abs(eig(D)));
    
    fprintf('change in poles is max_{1<=j<=r} |mu(j)^i - mu(j)^i-1|/|mu(j)^i-1 = %.12f \n', changein_errs(iter))
    fprintf(1, '----------------------------------------\n')

    %end clock
    fprintf(1, 'Current iterate finished in %.2f s\n', toc(iter_start))
    fprintf(1, 'End of current iterate k = %d\n', iter)
    fprintf(1, '---------------------------------------\n')
    
    % if convergence criterion is reached or maximal iterations have been
    % reached
    if changein_errs(iter) <= tol || iter == max_iter
        Sigma_opt = BQ_system(r,0);
        Sigma_opt.A = Ar;
        Sigma_opt.N = Nr;
        Sigma_opt.b = br;
        Sigma_opt.M = Mr; 
        return;
    end
    

end

