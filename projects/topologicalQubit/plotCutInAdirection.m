function f = plotCutInAdirection(F, a, b)   %F -- matrix; cut along line y=ax+b


[M,N] = size(F);
f=zeros(1,M);
x=1:M;
y=a*x+b;
for j=1:M
    n_min=floor(1+y(j));
    n_max=ceil(1+y(j));
    beta_=y(j)-n_min;
    alpha_=1-beta_;
    %n_min
    f(x(j)) = F(x(j),n_min)*alpha_+F(x(j),n_max)*beta_;
    
end

end
    

        