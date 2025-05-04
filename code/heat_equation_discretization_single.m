%create heat equation discretization
%SISO system, X is positive semidefinite solution
function [A,N,b,M] = heat_equation_discretization_single(k)

%specify grid refinement
h=1\(k+1);

%set up diag(1,-2,1)
a=-2;
Tk = diag(a*ones(1,k))+diag(ones(1,k-1),1) + diag(ones(1,k-1),-1);

%set up poisson matrix P
P = kron(speye(k),Tk) + kron(Tk,speye(k));

%set up the e_1,e_k,E_1,E_k according to the construction in the thesis
e1=zeros(k,1);
e1(1)=1;

ek = zeros(k,1);
ek(k)=1;

E1=e1*e1';
Ek=ek*ek';

e = ones(k,1);

%specify system matrices
A = (1\h^2)*(P + kron(E1,speye(k)));

%
%N has to be scaled with 0.1
N = 0.1*(1\h)*(kron(E1,speye(k)));

b = 0.1*sparse((1\h)*(kron(e1,e)));

%as output we choose the identity
M= kron(speye(k),speye(k));





