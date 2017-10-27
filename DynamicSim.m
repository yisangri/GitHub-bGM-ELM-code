disp(['Generating ',num2str(NS),' samples']);  

global scal
scal=0.004;%% scaling factor

main.td=35;          %[s]
main.cut=50*pi;       %cut-off frequeny of the process
main.dw=main.cut/1000;
main.nf=round(main.cut/main.dw);
main.dt=0.02;
main.t =main.dt:main.dt:main.td;

%% STEP1 Generate sample response process==================================
[out] = excit(main);  
tspan=linspace(main.dt,main.td,main.td/main.dt);
countd=0;
propBW
while countd<=NS

    countd=countd+1;        
    u=normrnd(0,1,[1 2*main.nf])';
    xo(countd,:)=out.s*u;

    main.Fe=[-xo(countd,:)];
    
    [tt uu]=ode45(@(t,y) BWmodel(t,y,main.Fe,main.t,str),[main.dt main.td],zeros(1,3));

    zo(countd,:)=interp1(tt,uu(:,1),main.t);
    zod(countd,:)=interp1(tt,uu(:,2),main.t);
    ssig(countd)=std(zo(countd,:));

    conv_check(countd)=std(ssig)/mean(ssig)/sqrt(countd);
    if countd>NS && conv_check(countd)<1.0e-1 
        break
    end
    disp(['Dynamic analysis ... ', num2str(countd)])

    %{
    close all     
    figure(); hold on; plot(HistVar.eps(1,:), HistVar.Q(1,:)*1.e-3,'k')
    grid on; box on; set(gcf,'color','w')
    xlabel('Deformation (m)'); ylabel('Restoring force (kN)');    
    title('Hysteresis loop'); set(gca,'fontsize',12);
    xlim([-0.012 0.012])
    ylim([-800 800])
    %}

end
    
clear conv ssig

% Determine how long it takes to achieve stationary
tem=std(zo)/mean(std(zo));
gfit = fittype( '1/(1+exp(a*t+b))','independent','t','coefficients',{'a','b'});
fitobject = fit((main.t/main.dt)',tem',gfit,'StartPoint',[-1,0] );
for jj=1:numel(main.t)
    if abs(1-fitobject(jj))<=1.e-2
        break
    end
end
main.ts=round(main.t(jj)/main.dt);


zf1=reshape(zo(:,main.ts:end)',1,[]);
zf2=reshape(zod(:,main.ts:end)',1,[]);

% From symmetricity
zf=[zf1 zf1 -zf1 -zf1; zf2 -zf2 zf2 -zf2]';


