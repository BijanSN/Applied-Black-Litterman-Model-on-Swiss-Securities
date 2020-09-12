% Applied Black-Litterman on Swiss portfolios
clear all, clc

% Load the data extracted from GetData.R and Data preprocessing

prices = readtable('Prices_SwissPortfolio.xlsx');
smi    = readtable('Benchmark.xlsx');
dates  = readtable('Dates.xlsx');
dates  = table2array(dates);

assetNames= prices.Properties.VariableNames;

prices.Dates = dates;
prices= movevars(prices, 'Dates', 'Before', 'ZURN_SW');

%% Summary statistics / Additional Data

ret = tick2ret(prices(:, 2:end)); % All but Dates (non-numerical)
assetRetns = ret(:, assetNames);
benchRet = tick2ret(smi);
numAssets = size(assetRetns, 2);
rf=0;

%%
%The inputs for P,q,and omega depends on the manager's views :

vprompt = 'How many views do you have ?\n';
v = input(vprompt);

%initialisation of P,Q and Omega :
P = zeros(v, numAssets);
q = zeros(v, 1);
Omega = zeros(v);

for i=1:v
    
    AssetViewed=0;
    AssetViewed2=0;
    absoluteStr=0;
    conf=0;
    expReturns=0;
    
    Display = [' Please indicate the instructions for the view number : ',num2str(v)];
    disp(Display)
    
    absolutePrompt='Is it an absolute view? Y/N [Y] \n';
    absoluteStr = input(absolutePrompt,'s');
    if isempty(absoluteStr)
        absoluteStr = 'Y';
    end
    
    % Asset ( Pick Matrix)
    if absoluteStr== 'Y' 
        AssetViewedPrompt='Which Asset is concerned? \n For example: ZURN_SW, NOVN_SW, NESN_SW, GIVN_SW \n';
        AssetViewed= input(AssetViewedPrompt,'s');
        
    else    
        AssetViewedPrompt='Which Asset is concerned first? \n For example: ZURN_SW, NOVN_SW, NESN_SW, GIVN_SW \n';
        AssetViewed= input(AssetViewedPrompt,'s');
         
        AssetViewedPrompt2='Which Asset is concerned second? \n For example: ZURN_SW, NOVN_SW, NESN_SW, GIVN_SW \n';
        AssetViewed2= input(AssetViewedPrompt2,'s');
    end

    if absoluteStr== 'Y'   
            P(i, assetNames==string(AssetViewed)) = 1; 
    else  
            P(i, assetNames==string(AssetViewed)) = 1; 
            P(i, assetNames==string(AssetViewed2)) = -1;
    end
   
    
    % View on annualised returns(Q)
   
    if absoluteStr== 'Y'   
        expReturnsprompt = 'What are your expectations regarding annual returns ? \n For example : 0.02 \n ';
        expReturns = input(expReturnsprompt);
    else
        RelexpReturnspromptdisp=['By how much will ', num2str(AssetViewed), ' outperform ', num2str(AssetViewed2), ' annual returns ? \nFor example : 0.02 \n '];
        disp(RelexpReturnspromptdisp)
        expReturnsprompt = '';
        expReturns = input(expReturnsprompt);
    end
       
    if isempty(expReturns)
         expReturns = '0.02';
    end
    q(i)=expReturns;
    
    % Confidence (Omega)
    confPrompt = 'How confident are you ?[Certain/Unsure] \n';
    conf = input(confPrompt,'s');
    
    if (conf==string('Certain'))
        Omega(i,i)=1e-5;
    else
        Omega(i,i)=1e-3;
    end
    
    fprintf('\n \n \n')
 end

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

variableNames = {'Asset_Ticker','MPT_Efficient_Allocation', 'Black_Litterman_Modified_Allocation'}
table(assetNames', wts, wtsBL,'VariableNames',variableNames)
