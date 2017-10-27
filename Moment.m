
function [output] = Moment(n,k,c)
global m sc

    kk=k;
    cc=c;
    H=@(w) sc./(kk-m*w.^2+cc*1i*w);
    func=@(w) autoPSD(w).*abs(H(w)).^2;
    output=2*integral(@(w)func(w).*w.^n,0,inf,'AbsTol',1.e-15,'RelTol',1.e-15);
