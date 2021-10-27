function [P,E,fun_value_difference,fun_value]=ADMM_l11_Huber(X,Y,pos)
% global I;
max_iter=pos.max_iter;
tolerance=pos.tolerance;
rho=pos.rho;
lambda1=pos.lambda1;
lambda2=pos.lambda2;
[~,task_num]=size(Y);
[m,n]=size(X);

P=zeros(n,task_num);
ZP=zeros(n,task_num);
UP=zeros(n,task_num);

E=zeros(m,task_num);
ZE=zeros(m,task_num);
UE=zeros(m,task_num);
[L U]=factor(X,rho);
I=zeros(m,m);
    
for iter=1:max_iter
    [L U]=factor(X,rho);
    [P,ZP,UP]=admm_l11_cycle(X,Y-E,P,ZP,UP,L,U,task_num,rho,lambda1,1e-1); 
    [E,ZE,UE]=admm_huber_cycle(I,Y-X*P,E,ZE,UE,task_num,rho,lambda2,1e-1);
    fun_value(iter)=objective_huber_l11(X,Y,P,E,lambda1,lambda2);
    if iter>=2
        fun_value_difference(iter-1)=abs((fun_value(end)-fun_value(end-1))/fun_value(end-1));
    end
    if tolerance~=0
        if iter>=2
            if fun_value_difference(iter-1)<tolerance
                break;
            end
        end
    end
end
end