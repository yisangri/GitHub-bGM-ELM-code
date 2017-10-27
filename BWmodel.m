

function dydt = BWmodel(t,y,GM,tG,str)


    xddg=interp1(tG,GM,t);
    
    dydt = zeros(3,1);
    dydt(1,1) = y(2);
    dydt(2,1) = - (str.k0*str.alp*y(1) + str.c0*y(2) + str.k0*(1-str.alp)*y(3))/str.m0 - xddg;
    dydt(3,1) = - str.gam*abs(y(2))*y(3)*abs(y(3))^(str.n-1) - str.eta*y(2)*abs(y(3))^str.n + str.A*y(2);
    
%     GMs=[GMs [xddg;t]];
%     ys=[ys [y(1);t]];
end