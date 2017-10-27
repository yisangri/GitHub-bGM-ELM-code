mphist=[];
upchists=[];
fpcmhists=[];
fpcmhists2=[];
zn=[];
ts=main.ts;
trialss=1;

zz=zo';
zzd=zod';

dth=0.03/50;
thres=0:dth:0.03;
    
for k=1:length(thres)
    count=0;
    samps1=[];
    samps2=[];
    samp1only=[];
    samp2only=[];
    for nex=2:size(zz,2)
        samps=[];           
        for nt=2:size(zz,1)
            theory=((zz(nt-1,nex))<thres(k) && (zz(nt,nex))>thres(k));
            if theory==1
                count=count+1;
            end
        end
    end
    upc(k)=count/(length(zz(:)))/main.dt;
     %disp([trialss k])        
end 

    
%% F-P

    fpcm=[];
    fpcm2=[];
    for k=1:length(thres) 
        count=0;
        count2=0;
        for nex=1:size(zz,2)
            if max(zz(:,nex))>thres(k)
                count=count+1;
            end
            if max(abs(zz(:,nex)))>thres(k)
                count2=count2+1;
            end
        end
        fpcm(k)=count/size(zz,2);
        fpcm2(k)=count2/size(zz,2);
    end

    
    for nex=1:size(zz,2)
        p_resp(nex)=max(abs(zz(:,nex)));  
    end
    
    upchists=[upchists upc'];
    fpcmhists=[fpcmhists fpcm'];
    fpcmhists2=[fpcmhists2 fpcm2'];
    mphist=[mphist mean(p_resp)];
    
%% Peak
    
%     countp=0;
%     countv=0;
%     for nex=2:size(zz,2)
%         peakss=[];
%         vallys=[];
% 
%         for nt=2:size(zz,1)       
%             if zzd(nt-1,nex)<0 && zzd(nt,nex)>0 
%                 countp=countp+1;
%                 peakss=[peakss ; (zz(nt-1,nex)+zz(nt,nex))/2];
%             end      
%             if zzd(nt-1,nex)>0 && zzd(nt,nex)<0 
%                 countv=countv+1;
%                 vallys=[vallys ; (zz(nt-1,nex)+zz(nt,nex))/2];
%             end        
%         end
%         msamps(nex,:)=[mean(abs(peakss)), mean(abs(vallys))];
%     end   
%     
%     mphists=[mphists mean(msamps(:))];
% 


mp_resp=mean(mphist)
save(['miniMCS']);

upc_mean=mean(upchists,2);
fpc_mean=mean(fpcmhists,2);
fpc_mean2=mean(fpcmhists2,2);

%% Plot

figure(2)
semilogy(thres(1:1:end),upc_mean(1:end)','o','MarkerSize',5,'MarkerFaceColor','w','DisplayName','MCS'); hold on

legend off; legend show;
set(gcf,'color','w')
% ylim([1.e-4 1])

figure(3)
semilogy(thres(1:1:end),(fpc_mean(1:end)),'o','MarkerSize',5,'MarkerFaceColor','w','DisplayName','MCS'); hold on
legend off; legend show;
set(gcf,'color','w')

