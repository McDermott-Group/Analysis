ff=zeros(51,50); 
for m=2:50 
    ff(:,m)=plotCutInAdirection(R,-0.02720,m);
end
 plot(ff+0.25*ones(51,50)*diag(1:50))