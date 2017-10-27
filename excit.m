

function [out] = excit2(main)
dt=main.dt;
dw =main.dw;
t=dt:dt:main.td;
w=dw:dw:dw*main.nf;
w=w-dw/2;
sigma=sqrt(2*autoPSD(w)*dw);
s1X=cos(t'*w)*diag(sigma);
s2X=sin(t'*w)*diag(sigma);
out.s=[s1X,s2X];
end

