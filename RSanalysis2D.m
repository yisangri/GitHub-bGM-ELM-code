%% Use peak factor (Vanmarcke's solution) to compute ordinates of the response spectrum 
clearvars d y

for k=1:numel(keq)
    peakf=peakfactor(tk,keq(k),ceq(k),mu(k,1),mu(k,2));
    sigx=sqrt(cov(1,1,k));
    d(k)=sigx*peakf;
    d2(k)=alp(k)*d(k)^2*(1+mu(k,1)^2/sig(k)^2);
end
[mzk f1]=sort(d+abs(mu(:,1))');
alpk=alp(f1);
f2=mzk>=0.95*alp*(d+abs(mu(:,1))')';
mz=mzk(f2);
maxz=(alpk(f2)/sum(alpk(f2)))*mz';
maxz2=sqrt(sum(d2));
fprintf('\t Mean peak response:\t%1.4f  \t%1.4f \n',maxz,maxz2);

