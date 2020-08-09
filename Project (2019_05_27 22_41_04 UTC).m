clear all,clc
load dbTreasuryCurve.mat % ParFitted = Yield curve
N=size(ParFitted,2); 
YieldChanges=ParFitted(2:end,:)-ParFitted(1:end-1,:); % "returns"/ yield change
CurrentYields=ParFitted(end,:); % last yield

%% Our data
plot(ParFitted) 
%montrer inversion yield curve ici :
surf(ParFitted)
mesh(ParFitted)
%% linear dependance of the data : the more the yield are far appart from each other, the less than seem to behave LINEARLY.
Corr=corr(YieldChanges)
imagesc(Corr);
colormap(jet);
colorbar;
clf;
%% using only 2 maturities ( 6 and 2 years)
%paramhat = copulafit('t',u) %  returns an estimate, rhohat, of the matrix of linear correlation parameters for a t copula, and an estimate of the degrees of freedom parameter, nuhat, given the data in u.
%yield 6mois et 2 ans
scatterhist(YieldChanges(:,1),YieldChanges(:,2))

%Transform the data to the copula scale (unit square) using a kernel
%estimator of the cumulative distribution function
u = ksdensity(YieldChanges(:,1),YieldChanges(:,1),'function','cdf');
v = ksdensity(YieldChanges(:,2),YieldChanges(:,2),'function','cdf');

% Get the margins
figure;
scatterhist(u,v)
xlabel('u')
ylabel('v')
hold on % t-distrib could be fitted => dependances on both sides !
clf;
%% compute empyrical CDF (unsmoothed)
nobs = size(YieldChanges,1);

invCDF1 = sort(YieldChanges(:,1));
n1 = length(YieldChanges(:,1));
invCDF2 = sort(YieldChanges(:,2));
invCDF3 = sort(YieldChanges(:,3));
invCDF4 = sort(YieldChanges(:,4));
invCDF5 = sort(YieldChanges(:,5));
invCDF6 = sort(YieldChanges(:,6));

stairs(invCDF1,(1:nobs)/nobs)
hold on;
stairs(invCDF2,(1:nobs)/nobs)
stairs(invCDF3,(1:nobs)/nobs)
stairs(invCDF4,(1:nobs)/nobs)
stairs(invCDF5,(1:nobs)/nobs)
stairs(invCDF6,(1:nobs)/nobs)


legend('3M','2Y','5Y','10Y','20Y','30Y');
xlabel('X');
ylabel('Cumulative Probability');
clf;
%% Alternatively

% empirical cdf,
% for a continuous distribution, smoother model than the step function computed by ecdf. 
% Another way is to use kernel smoothing with ksdensity. 

[f1,x1] = ecdf(YieldChanges(:,1),'Bounds','on')
[f2,x2] = ecdf(YieldChanges(:,2),'Bounds','on')
[f3,x3] = ecdf(YieldChanges(:,3),'Bounds','on')
[f4,x4] = ecdf(YieldChanges(:,4),'Bounds','on')
[f5,x5] = ecdf(YieldChanges(:,5),'Bounds','on')
[f6,x6] = ecdf(YieldChanges(:,6),'Bounds','on')


% Graphs with bounds :
close all
ecdf(YieldChanges(:,1),'Bounds','on')
%stairs(x1,f1,'b','LineWidth',2) same graph
hold on
%ecdf(YieldChanges(:,1),'Bounds','on') with confidence interval

%plot(x1,f)
hold on
ecdf(YieldChanges(:,2),'Bounds','on')
ecdf(YieldChanges(:,3),'Bounds','on')
ecdf(YieldChanges(:,4),'Bounds','on')
ecdf(YieldChanges(:,5),'Bounds','on')
ecdf(YieldChanges(:,6),'Bounds','on')

% F1 = ksdensity(YieldChanges(:,1),x1,'function','cdf','width',.15)
% F2 = ksdensity(YieldChanges(:,2),x2,'function','cdf','width',.15)
% F3 = ksdensity(YieldChanges(:,3),x3,'function','cdf','width',.15)
% F4 = ksdensity(YieldChanges(:,4),x4,'function','cdf','width',.15)
% F5 = ksdensity(YieldChanges(:,5),x5,'function','cdf','width',.15)
% F6 = ksdensity(YieldChanges(:,6),x6,'function','cdf','width',.15)
clf;
%% Smoothed CDF
plot(x1,f1) % smoothed
hold on 
plot(x2,f2)
plot(x3,f3)
plot(x4,f4)
plot(x5,f5)
plot(x6,f6)
%integrate and get the PDF : not smooth enough  => use kernel.
clf;
%% Comparaison Gaussian and Real 3M(using gaussian kernel) PDF :

for i=1:1015
     r(i)=random('normal',0,0.05)
end
r=r/80
[fn,xn]=ksdensity(r)
plot(xn,fn,'r')
hold on
ksdensity(YieldChanges(:,1))

clf;
%% Density of the change in yield depending on maturities using non param' gaussian kernel
subplot(2,3,1)
%[pdf1]=gaussiankernel(ParFitted(:,1),1017)
%plot(pdf1(1,:),pdf1(2,:))
[f11,x11] = ksdensity(YieldChanges(:,1))
plot(x11,f11)
title('3months change of  Yield');
xlabel('Values');
ylabel('Density');


subplot(2,3,2)
% [pdf2]=gaussiankernel(ParFitted(:,2),1017)
% plot(pdf2(1,:),pdf2(2,:));
[f22,x22] = ksdensity(YieldChanges(:,2))
plot(x22,f22)
title('2year Yield ');
xlabel('Values');
ylabel('Density');


subplot(2,3,3)
% [pdf3]=gaussiankernel(ParFitted(:,3),1017)
% plot(pdf3(1,:),pdf3(2,:));
[f33,x33] = ksdensity(YieldChanges(:,3))
plot(x33,f33)
title('5 years Yield');
xlabel('Values');
ylabel('Density');

subplot(2,3,4)
% [pdf4]=gaussiankernel(ParFitted(:,4),1017)
% plot(pdf4(1,:),pdf4(2,:));
[f44,x44] = ksdensity(YieldChanges(:,4))
plot(x44,f44)

title('10 years Yield ');
xlabel('Values');
ylabel('Density');

subplot(2,3,5)
% [pdf5]=gaussiankernel(ParFitted(:,5),1017)
% plot(pdf5(1,:),pdf5(2,:));
[f55,x55] = ksdensity(YieldChanges(:,5))
plot(x55,f55)
title('20 years Yield ');
xlabel('Values');
ylabel('Density');

subplot(2,3,6)
% [pdf6]=gaussiankernel(ParFitted(:,6),1017)
% plot(pdf6(1,:),pdf6(2,:));
[f66,x66] = ksdensity(YieldChanges(:,6))
plot(x66,f66)

title('30 years Yield ');
xlabel('Values');
ylabel('Density');

close all;
%% test of normality 
jbtest(YieldChanges(:,1)) %reject normality 
jbtest(YieldChanges(:,2)) %reject normality 
jbtest(YieldChanges(:,3)) %reject normality
jbtest(YieldChanges(:,4)) %reject normality
jbtest(YieldChanges(:,6)) %reject normality

%multivariate test 
H= mardiatest(YieldChanges,0.05); % reject normality as well
%% Kolmogorov test to see if cdf are differents or not
%to do
%%
% fit t coupla to yield changes (bivariate here) 
[Rho,nu] = copulafit('t',[u v],'Method','ApproximateML')
%% Fit Copula on whole data set
[Nu,Mu,Sigma]=FitT(YieldChanges)
[s,C]=cov2corr(Sigma);

%% Estimate (smoothed) marginals and store result as quantiles ***adapted from Meucci(see references)***
Xs=[];
Fs=[];


for n=1:N
    % truncate empirical distribution where yield changes give rise to
    % negative yields
    DY=YieldChanges(:,n);
    NonFeasible=find(CurrentYields(n)+DY<0);
    DY(NonFeasible)=[];
    L=length(DY);

    % kernel smoothing of the truncated empirical cdf ( as done before)
    X_Low=median(DY)-1.2*(median(DY)-min(DY));
    X_High=median(DY)+1.2*(max(DY)-median(DY));

    Step=(X_High-X_Low)/5000;
    X=[X_Low : Step : X_High]';

    eps=.3*std(DY);
    F=0*X;
    for l=1:L
        F=F+normcdf(X,DY(l),eps);
    end
    F=F/(L+1);
    F(end)=1.001;

    Xs=[Xs X];
    Fs=[Fs F];
end

plot(Xs,Fs)
ylim([0 1])

close all;
%% generate Monte Carlo simulations of the copula
J=50000;
tSample = mvtrnd(C,Nu,J/2);
tSample = [tSample
    -tSample];

U=tcdf(tSample,Nu);

% map Monte Carlo simulations of the copula into Monte Carlo
% simulations for the joint distribution
MPrior=[];
for n=1:N
    yi=.001+.998*U(:,n);
    y=Fs(:,n);
    x=Xs(:,n);
    xi = interp1(y,x,yi);
    MPrior=[MPrior xi];
end
%% Enter views and confidence ***from Meucci(2006)***

[J,N]=size(MPrior) % J=50k simulations of 6 assets(=N)

% initialize pick matrix
K=2;
P=zeros(K,N); 

% view on the 2-5-10 butterfly 
P(1,2)=-.5; P(1,3)=1; P(1,4)=-.5;
a_b(1,:)=.0001*[0 5];
c(1)=.25;

% view on slope
P(2,2)=-1; P(2,5)=1;
a_b(2,:)=.0001*[0 10];
c(2)=.25;


% compute posterior

% Step 1: rotate the market into the views' coordinates
P_bar=[P; null(P)'];
V=MPrior*P_bar';

% Step 2: sort the panel of Monte Carlo views scenarios
[W,C]=sort(V(:,[1:K]));
for k=1:K
    x=C(:,k);
    y=[1:J];
    xi=[1:J];
    yi = interp1(x,y,xi);
    C(:,k)=yi/(J+1);
end
F=zeros(50000,2);
for k=1:K
    % Step 3: compute the marginal posterior cdf of each view
    F(:,k)=[1:J]'/(J+1);
    F_hat(:,k)=unifcdf(W(:,k),a_b(k,1),a_b(k,2));
    F_tilda(:,k)=(1-c(k))*F(:,k)+c(k)*F_hat(:,k);
    
    % Step 4: compute joint posterior realizations of the views
    dummy = interp1(F_tilda(:,k),W(:,k),C(:,k),'linear','extrap');
    V_tilda(:,k) = dummy;
end

% Step 5: compute joint posterior realizations of the market distribution
V_tilda=[V_tilda V(:,[K+1:end])];
MPost=V_tilda*inv(P_bar');

% map posterior into returns of securities

load dbPricingFeatures
J=size(MPrior,1); % number of Monte Carlo simulations
N=size(MPrior,2); % number of market factors
L=length(WeeklyCarries);  % number of securities ( 2,5,10 and 30 only)

% mapping the prior
LinearPr = -MPrior*KRDs';
ConvexityFactorPr = (mean(MPrior,2)).^2;
QuadraticPr = ConvexityFactorPr*OACs';
RetsPr = ones(J,1)*WeeklyCarries' + LinearPr + QuadraticPr;

% mapping the posterior
LinearPost = -MPost*KRDs';
ConvexityFactorPost = (mean(MPost,2)).^2;
QuadraticPost = ConvexityFactorPost*OACs';
RetsPost = ones(J,1)*WeeklyCarries' + LinearPost + QuadraticPost;

%% Summary statistics
mean_prior= mean(RetsPost); 
std_prior= std(RetsPost);
kurtosis_prior= kurtosis(RetsPost);
skewness_prior= skewness(RetsPost);


mean_posterior= mean(RetsPr);
std_posterior= std(RetsPr);
kurtosis_posterior= kurtosis(RetsPr);
skewness_posterior= skewness(RetsPr);

%% tables 
% table(mean_prior ,mean_posterior);
% table(std_prior ,std_posterior);
% table(kurtosis_prior ,kurtosis_posterior);
% table(skewness_prior ,skewness_posterior);
%% Using posterior yield change( not working )

% Tau=0.05;
% Sigma= cov(YieldChanges);
% S = size(P_bar*Sigma*P_bar',1);
% BSB         = (Tau*P_bar*Sigma*P_bar') \ eye(S,S)
% 
% BL_mu = RetsPost + Tau*Sigma*P_bar'*BSB*(YieldChanges - P_bar*RetsPost);
% BL_sigma = (1 + Tau).*Sigma - (Tau^2).*Sigma*B'*BSB*(em_r - B*Mu_Eq);
% BL_sigma_inv = BL_sigma\eye(23);
% BL_weights = 1/gamma*BL_sigma_inv*BL_mu; 
% 
% 
% display('Black Litterman weights with COP views, tau = 1/T')
% array2table(BL_weights,'RowNames',Names)
% 
% figure(1)
% x =(1:N)';
% 
% bar(x,[BL_weights RetsPost])
% XTickLabel = Names;
% XTick      = 1:N;
% XTickLabelRotation = 60;
% legend('Black Litterman weights','MC weigths')
% title ('Black Littermann COP approach ')
