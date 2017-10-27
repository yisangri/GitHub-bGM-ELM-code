
function [output] = Moment1D(n,k,c)
global m c sc

    kk=k;
    H=@(w) sc./(kk-m*w.^2+c*1i*w);
    func=@(w) autoPSD(w).*abs(H(w)).^2;
    output=2*integral(@(w)func(w).*w.^n,0,inf,'AbsTol',1.e-15,'RelTol',1.e-15);
