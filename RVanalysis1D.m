
global sigk m c sc

%% STEP4 : mean-crossing rate && first passage probability (old)
sc=str.m0;
m=str.m0; 
c=str.c0; 
% tk=(numel(main.t)-main.ts)*main.dt;
tk=27;

dth=0.03/50;
thr=0:dth:0.03;
 
fp=zeros(1,numel(thr));
x0=[str.k0];
options = optimset('TolX',1e-12,'TolFun',1e-12,'display','off');
keq=[];
for l=1:Ko
    sigk=sqrt(covU(l));    
    [keq(l) fval(l) fl(l)]=fminsearch(@objectf1D,x0,options); % When it does not converge, check autopsd func
    disp(['finding ', num2str(l) '-th 1D ELS, error: ' num2str(fval(l))]);
end
keq=(keq);
 
for tr=1:size(thr,2)
    z0=thr(tr);
    tem=0;
    for kk=1:Ko
        lam0=Moment1D(0,keq(kk));%% Compute the spectral moment check autopsd func
        lam2=Moment1D(2,keq(kk));
        tem(kk)=1/(2*pi)*sqrt(lam2/lam0)*exp(-0.5*(z0-muU(kk))^2/lam0); %% Compute the crossing rate for each sublinear system
    end
    cr1D(tr)=alpU*tem'; %% compute the crossing rate  
end 

fpc1D=1-exp(-cr1D*tk);%% compute the first passage probability

Mcr = load(['CR_MCS.txt']);%%load the MCS solution for crossing rate
Ecr = load(['CR_ELM.txt']);%%load the ELM solution for crossing rate
Tcr = load(['CR_TELM.txt']);%%load the TELM solution for crossing rate
Mfp = load(['FP_MCS.txt']);%%load the MCS solution for FP probability
Efp = load(['FP_ELM.txt']);%%load the ELM solution for FP probability
Tfp = load(['FP_TELM.txt']);%%load the TELM solution for FP probability


figure(2)  
semilogy(thr,cr1D,'-.','linewidth',2,'DisplayName','GM-ELM (univariate)');
grid on; hold on;
semilogy(Ecr(:,1),Ecr(:,2),':','linewidth',2,'DisplayName','ELM');
semilogy(Tcr(:,1),Tcr(:,2),'--','linewidth',2,'DisplayName','TELM');
semilogy(Mcr(:,1),Mcr(:,2),'ko','linewidth',2,'DisplayName','MCS');
legend off; legend show;

figure(3)
semilogy(thr,fpc1D,'-.','linewidth',2,'DisplayName','GM-ELM (univariate)');
grid on; hold on;
semilogy(Efp(:,1),Efp(:,2),':','linewidth',2,'DisplayName','ELM');
semilogy(Tfp(:,1),Tfp(:,2),'--','linewidth',2,'DisplayName','TELM');
semilogy(Mfp(:,1),Mfp(:,2),'ko','linewidth',2,'DisplayName','MCS');
legend off; legend show;






