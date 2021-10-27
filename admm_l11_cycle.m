% minimize_x (1/2)||Ax-b||_2^2+lambda||x||_1
% A=X;
% b=Y-X*P-E;
% z=ZQ;
% [L U]=factor(A,rho);
% [m,n]=size(A);
function [x,z,u]=admm_l11_cycle(A,b,x,z,u,L,U,task_num,rho,lambda,tolerance)
[m,n]=size(A);
Atb=A'*b;

fun_value=[];
fun_value_difference=[1];
while fun_value_difference(end)>tolerance
    %x-update
    q=Atb+rho*(z-u);
    if(m>=n)
        x=U\(L\q);
    else
        x=q/rho-(A'*(U\(L\(A*q))))/rho^2;
    end
%     x=Quadratic_Objective_Term(A,b,x,z,u,L,U,Atb,rho,m,n);
    %z-update
    z=shrinkage(x+u,lambda/rho);
    %u-update
    u=u+(x-z);
    temp=A*x-b;
    fun_value=cat(1,fun_value,0.5*norm(temp(:))^2+sum(sum(abs(x))));
    if length(fun_value)>=2
        fun_value_difference=cat(1,fun_value_difference,abs((fun_value(end)-fun_value(end-1))/fun_value(end-1)));
    end
end
end

