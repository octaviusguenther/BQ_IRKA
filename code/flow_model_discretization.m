%create flow model discretization
%SISO system
%The resulting model is of dimension k+k^2
function [A,N,b,M] = flow_model_discretization(k)

%viscosity coefficient
nu = 0.1;

%stepsize
h = 1/(k+1);

%constructing the matrices A_1, A_2, B_0,B_1 and M
%initialize Poisson Matrix and define A_1
P = diag(-2*ones(1,k))+diag(ones(1,k-1),1) + diag(ones(1,k-1),-1);

A_1 = nu/(h^2)*P;

% initialize dimension of A_2
A_2 = zeros(k,k^2);

%A_2 is given as A_2 = [A_2^1,..., A_2^k,...,A_2^N]
% fill in the entries
% A_2^1 
A_21=zeros(k,k);
A_21(1,2) = -1;
A_21(2,1) = -1;

%A_2^N
A_2N = zeros(k,k);
A_2N(k-1,k)=1;
A_2N(k,k-1)=1;

%A_2^l, 1<l<N
%!!! one additional matrix in this sequence, will be truncated at the end
A_2l = zeros(k,k,k-1);
for l=2:size(A_2l,3)
    for i=1:k
        for j=1:k
            if (i == l-1 && j == l) 
                A_2l(i,j,l)=1;
            elseif (i == l && j == l-1)
                A_2l(i,j,l)=1;
            elseif (i==l && j == l +1)
                A_2l(i,j,l)=-1;
            elseif (i==l+1 && j ==l)
                A_2l(i,j,l)=-1;
            else 
                A_2l(i,j,l)=0;
            end
        end
    end

end

%add A_21 and A_2N to the sequence of matrices 
A_2l(:,:,1) = A_21;
A_2l(:,:,k) = A_2N;

%reshape the sequence into a matrix of dimension N + N^2
A_2 = reshape(A_2l(:,:,:),[k,k^2]);

%scale with h_2
A_2 = (1/(2*h))*A_2;

%define B_0
B_0 = zeros(k,1);
B_0(1) = nu/(h^2);

%define B_1
B_1 = zeros(k);
B_1(1,1) = 1/(2*h);


A = [sparse(A_1),sparse(A_2);sparse(k^2,k),kron(A_1,speye(k))+kron(speye(k),A_1)];
N=[sparse(B_1),sparse(k,k^2);kron(B_0,speye(k))+kron(speye(k),B_0),sparse(k^2,k^2)];
b=[sparse(B_0);sparse(k^2,1)];
M=[speye(k),sparse(k,k^2);sparse(k^2,k),sparse(k^2,k^2)];

