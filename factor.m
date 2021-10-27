function [L U]=factor(A,rho)
[m,n]=size(A);
if (m>=n)% if skinny
    L=chol(A'*A+rho*speye(n),'lower');
else
    L=chol(speye(m)+1/rho*(A*A'),'lower');
end
% force matlab to recoginize the upper/lower triangular structure
L=sparse(L);
U=sparse(L');
end