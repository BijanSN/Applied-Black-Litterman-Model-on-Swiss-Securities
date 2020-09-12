% Applied Black-Litterman on Swiss portfolios
clear all, clc

% Load the data extracted from GetData.R and Data preprocessing

prices = readtable('Prices_SwissPortfolio.xlsx');
smi    = readtable('Benchmark.xlsx');
dates  = readtable('Dates.xlsx');
dates  = table2array(dates)

assetNames= prices.Properties.VariableNames;

prices.Dates = dates
prices= movevars(prices, 'Dates', 'Before', 'ZURN_SW')
head(prices)

%% Summary statistics / Additional Data

ret = tick2ret(prices(:, 2:end)); % All but Dates (non-numerical)
assetRetns = ret(:, assetNames);
benchRet = tick2ret(smi);
numAssets = size(assetRetns, 2);
rf=0;

%%
%The inputs for P,q,and omagea depends on the manager's views.

v = 3;  % total 3 views
P = zeros(v, numAssets);
q = zeros(v, 1);
Omega = zeros(v);

% View 1 : Nestlé SA is going to have 2% annual return with uncertainty 1e-3.
%          This is a absolute view, with high uncertainty.

P(1, assetNames=="NESN_SW") = 1; 
q(1) = 0.02;
Omega(1, 1) = 1e-3;

% View 2 : Givaudan SA is going to have 4% annual return with uncertainty 1e-3.
%          This is a absolute view, with high uncertainty.

P(2, assetNames=="GIVN_SW") = 1; 
q(2) = 0.04;
Omega(2, 2) = 1e-3;

% View 3 : Novartis SA is going to outperform Credit suisse Group by 1% annual return with uncertainty 1e-5.
%          This is a relative view ,with low uncertainty.
P(3, assetNames=="NOVN_SW") = 1; 
P(3, assetNames=="CSGN_SW") = -1; 
q(3) = 0.01;
Omega(3, 3) = 1e-5;

% Summary : 
T = array2table([P q diag(Omega)])
T.Properties.VariableNames = {'ZURN_SW' 'UHR_SW' 'SGSN_SW' 'ROG_SW' 'CFR_SW' 'NOVN_SW' 'NESN_SW' 'LONN_SW' 'LHN_SW' 'GIVN_SW' 'GEBN_SW' 'CSGN_SW' 'Returns' ,'Uncertainty'}

% convert anuualised views into daily
convertBusinessdays = 1/252;
q = q*convertBusinessdays; 
Omega = Omega*convertBusinessdays;

%covariance matrix of assets returns :
Sigma = cov(assetRetns.Variables);
tau = 1/size(assetRetns.Variables, 1);

%C= uncertainty*covar
C = tau*Sigma; 

%% Pi is the implied equilibrium return and C is the uncertainty in our prior.

[wtsMarket, PI] = MarketPortfolioAndImpliedReturn(assetRetns.Variables, benchRet.Variables);

mu_bl = (P'*(Omega\P) + inv(C)) \ ( C\PI + P'*(Omega\q));
cov_mu = inv(P'*(Omega\P) + inv(C));

varNames = {'Asset_Ticker','Prior_Belief_of_Expected_Return', 'Black_Litterman_Blended_Expected_Return'}
T= table(assetNames', PI*252, mu_bl*252,'VariableNames',varNames)%,{'Variable Names',})

%% 
port = Portfolio('NumAssets', numAssets, 'lb', 0, 'budget', 1, 'Name', 'MPT-Efficient Allocation');
port = setAssetMoments(port, mean(assetRetns.Variables), Sigma);
wts = estimateMaxSharpeRatio(port);

portBL = Portfolio('NumAssets', numAssets, 'lb', 0, 'budget', 1, 'Name', 'Black-Litterman-Modified Allocation');
portBL = setAssetMoments(portBL, mu_bl, Sigma + cov_mu);  
wtsBL = estimateMaxSharpeRatio(portBL);

ax1 = subplot(1,2,1);
idx = wts>0.001;
pie(ax1, wts(idx), assetNames(idx));
title(ax1, port.Name ,'Position', [-0.05, 1.6, 0]);

ax2 = subplot(1,2,2);
idx_BL = wtsBL>0.001;
pie(ax2, wtsBL(idx_BL), assetNames(idx_BL));
title(ax2, portBL.Name ,'Position', [-0.05, 1.6, 0]);

variableNames = {'Asset_Ticker','MPT-Efficient Allocation', 'Black-Litterman-Modified Allocation'}
table(assetNames', wts, wtsBL,'VariableNames',variableNames)
