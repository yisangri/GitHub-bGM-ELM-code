clear all
tic  

%% Generate sample response processes
NS=100;
run DynamicSim
zf_org=zf;

%% GM parameter estimation for bivariate GM-ELM
disp('Expectation Maximization algorithm (2D):')
Ko=64;
input.z=zf;
input.K=Ko;

% Define initial value
bd=std(zf);
Kh=sqrt(Ko);
bdx=bd(1):(-bd(1)-bd(1))/(Kh-1):-bd(1);
bdy=bd(2):(-bd(2)-bd(2))/(Kh-1):-bd(2);

count=0;
for i=1:Kh
    for j=1:Kh
        count=count+1;
        S.mu(count,:)=3*[bdx(j)+normrnd(0,bd(1)*0.1) bdy(i)+normrnd(0,bd(2)*0.1)];
    end
end
for i=1:Ko
    S.Sigma(:,:,i)= mean(var(zf,1))*eye(2);
    S.ComponentProportion(i)=1/Ko;
end

% Iteration

options=statset('Display','iter','MaxIter',10000,'TolFun',1e-7);
GMmodel=fitgmdist(input.z,input.K,'Options',options,'Replicates',1,'start',S);
cov=GMmodel.Sigma;
alp=GMmodel.ComponentProportion;
mu=GMmodel.mu;
fprintf('\t number of mixture component :\t%1.0f\n', Ko);
fprintf('\t total number of dyanmic analysis : \t%1.0f\n', countd);
run RVanalysis2D

%% GM parameter estimation for univariate GM-ELM
bdxu=bd(1):(-bd(1)-bd(1))/(Ko-1):-bd(1);
SU.mu=zeros(Ko,1);
for i=1:Ko
    SU.mu(i)= 3*(bdxu(j)+normrnd(0,bd(1)*0.1));
    SU.Sigma(:,:,i)= mean(var(zf(:,1),1));
    SU.ComponentProportion(i)=1/Ko;
end

disp('Expectation Maximization algorithm (1D):')
options=statset('Display','iter','MaxIter',10000,'TolFun',1e-7);
GMmodelU=fitgmdist(input.z(:,1),input.K,'Options',options, 'Replicates',1 ,'start',SU);
covU=GMmodelU.Sigma;
alpU=GMmodelU.ComponentProportion;
muU=GMmodelU.mu;

%% FIGURE of inverse cdf
vdim=1;
if vdim==1
    xx=-0.02:0.001:0.02;
elseif vdim==2
    xx=0:0.001:0.8;
end

sig=reshape(sqrt(cov(vdim,vdim,:)),[],1);
cdf=0; pdf=0;
for j=1:numel(xx)
    cdf(j)=sum(alp'.*normcdf(xx(j),mu(:,vdim),sig));
    pdf(j)=sum(alp'.*normpdf(xx(j),mu(:,vdim),sig));
end

figure(1)
subplot(1,2,1)
semilogy(xx,pdf,'LineWidth',2); title('PDF'); hold on
subplot(1,2,2)
semilogy(xx,1-cdf); title('Complementary Cumulative Distribution'); hold on;

%% Random Vibration Anlaysis
run RVanalysis2D
run RSanalysis2D
run RVanalysis1D

minu=round(toc/60);
disp(num2str(minu));

figure(3)
legend off
figure(2)
legend off

