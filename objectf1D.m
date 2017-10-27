
function [out]=objectf1D(param)
    global sigk m sc c
    kk=param(1);
    H=@(w) sc./(kk-m*w.^2+c*1i*w);
    func=@(w) autoPSD(w).*abs(H(w)).^2;
    temp1=2*integral(func,0,inf,'AbsTol',1.e-15,'RelTol',1.e-15);
    out=abs(temp1-sigk^2)/sigk^2;

