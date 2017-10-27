function [peak] = peakfactor(tk,keq,ceq,mux,mud) 

    lam0=Moment(0,keq,ceq);
    lam1=Moment(1,keq,ceq);
    lam2=Moment(2,keq,ceq);
    delt=sqrt(1-lam1^2/lam0/lam2);

    nuo=1/(pi)*sqrt(lam2/lam0)*exp(-0.5*(0)^2/lam0); % Both sides
    nu_t=tk*nuo/(-log(0.5704));
    peak=sqrt(2*log(nu_t*(1-exp(-delt^1.2*sqrt(pi*log(nu_t))))));
    
end