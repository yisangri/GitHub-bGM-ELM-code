%% STEP4 : mean-crossing rate && first passage probability (old)


global sigk sigdk m sc

sc=str.m0; %% see Eq.(18) in the GM-ELM paper
m=str.m0;
% tk=(numel(main.t)-tsm)*main.dt;
tk=27;
dth=0.03/50;
thr=0:dth:0.03;

fp=zeros(1,numel(thr));
x0=[str.k0 str.c0];

options = optimset('TolX',1e-12,'TolFun',1e-12,'display','off');
for kk=1:Ko
    sigk=sqrt(cov(1,1,kk));
    sigdk=sqrt(cov(2,2,kk));     
    [parameq(kk,:) fval(kk) fl(kk)]=fminsearch(@objectf,x0,options); 
    disp(['finding ', num2str(kk) '-th 2D ELS, error: ' num2str(fval(kk))]);
end

keq=(parameq(:,1));
ceq=(parameq(:,2));

for tr=1:size(thr,2)
    z0=thr(tr);
    tem_cr=0;
    safe_start_eta=0;
    for kk=1:Ko        
        lam0=Moment(0,keq(kk),ceq(kk));
        lam2=Moment(2,keq(kk),ceq(kk));
        mux=mu(kk,1);
        mud=mu(kk,2);
        delt=mud/sqrt(lam2);
        lam2t=lam2*2*pi*(normpdf(delt)+delt-delt*normcdf(-delt)).^2;     
        tem_cr(kk)=1/(2*pi)*sqrt(lam2t/lam0)*exp(-0.5*(z0-mux)^2/lam0);
    end   
    cr1(tr)=alp*tem_cr'; 
end 

fpc_pois1=1-exp(-cr1*tk)';%% compute the first passage probability

for tr=1:size(thr,2)
    z0=thr(tr);
    tem_cr=0;
    safe_start_eta=0;
    for kk=1:Ko        
        lam0=Moment(0,keq(kk),ceq(kk));
        lam2=Moment(2,keq(kk),ceq(kk));
        mux=mu(kk,1);
        mud=mu(kk,2);
        delt=mud/sqrt(lam2);
        lam2t=lam2*2*pi*(normpdf(delt)+delt-delt*normcdf(-delt)).^2;     
        tem_cr(kk)=1/(2*pi)*sqrt(lam2t/lam0)*exp(-0.5*(z0-mux)^2/lam0);
    end   
    cr2(tr)=alp*tem_cr'; 
end 

fpc_pois2=1-exp(-cr2*tk)';%% compute the first passage probability



fpc_pois=(fpc_pois1+fpc_pois2)/2;
cr=(cr1+cr2)/2;




figure(2)  
semilogy(thr,cr2,'-','LineWidth',2,'DisplayName','GM-ELM (bivariate)')
hold on; grid on;
ylim([1.e-5 10]); xlim([0 0.02])
xlabel('Displacement'); ylabel('Crossing Rate');
legend off; legend show;
set(gcf,'color','w')

figure(3)
semilogy(thr,fpc_pois2,'-','LineWidth',2,'DisplayName','GM-ELM (bivariate)')
ylim([1.e-4 1]);  xlim([0 0.02])
xlabel('Displacement')
ylabel('First-passage Probability')
legend off; legend show;
grid on; hold on;
set(gcf,'color','w')

