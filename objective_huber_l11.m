function value=objective_huber_l11(X,Y,P,E,lambda1,lambda2)
temp=Y-X*P;
value=0.5*norm(temp(:))^2+lambda1*sum(abs(P(:)))+sum(huber(E(:),lambda2));
end