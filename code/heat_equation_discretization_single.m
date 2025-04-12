%create heat equation discretization
%SISO system, X is positive semidefinite solution
function [A,N,b,M] = heat_equation_discretization_single(k)

h=1\(k+1);

a=-2;
Tk = diag(a*ones(1,k))+diag(ones(1,k-1),1) + diag(ones(1,k-1),-1);

P = kron(speye(k),Tk) + kron(Tk,speye(k));

e1=zeros(k,1);
e1(1)=1;

ek = zeros(k,1);
ek(k)=1;

E1=e1*e1';
Ek=ek*ek';

e = ones(k,1);

A = (1\h^2)*(P + kron(E1,speye(k)));


N = (1\h)*(kron(E1,speye(k)));


b = sparse((1\h)*(kron(e1,e)));
M= kron(speye(k),speye(k));





