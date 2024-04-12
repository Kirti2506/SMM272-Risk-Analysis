clear all; close all; clc
rng(123)

%Log-Returns
data = xlsread("CW_Q4_new.xlsm","Log-returns");
AA_return = data(:,2);
INTC_return = data(:,3);
JPM_return = data(:,4);
PG_return = data(:,5);
historical_returns = data(:,2:5);

% Mean returns
AA_mean = mean(AA_return);
INTC_mean = mean(INTC_return);
JPM_mean = mean(JPM_return);
PG_mean = mean(PG_return);

% Spot Prices
AA_price = 27.96;
INTC_price = 43.47;
JPM_price = 180.9;
PG_price = 160.4;

% Strike Prices
AA_K = AA_price*1.05;
INTC_K = INTC_price*0.9;
JPM_K = JPM_price*1.0;
PG_K = PG_price*1.10;

% Option Prices (see Excel tab "2. Black-Scholes")
Call_AA = 5.31;
Put_AA = 5.55;
Call_INTC = 7.48;
Put_INTC = 1.98;
Call_JPM = 12.58;
Put_JPM = 9.00;
Call_PG = 4.85;
Put_PG = 15.67;

alpha = 0.99;

% 1. Covariance Matrix
covariance_matrix = cov(historical_returns);

% 2. Generate Random Numbers
num_simulations = 10000;
mean_returns = [AA_mean,INTC_mean,JPM_mean,PG_mean]; 
stock_returns = mvnrnd(mean_returns, covariance_matrix, num_simulations);

% 3. Scaling of returns
time_horizon = 10;
scaled_returns = stock_returns * sqrt(time_horizon);

% 4. Simulated prices
AA_sim_prices = AA_price*exp(stock_returns(:,1));
INTC_sim_prices = INTC_price*exp(stock_returns(:,2));
JPM_sim_prices = JPM_price*exp(stock_returns(:,3));
PG_sim_prices = PG_price*exp(stock_returns(:,4));

% 5. Calculation of Option Prices on Day 10
[C_AA,P_AA] = blsprice(AA_sim_prices,AA_K,0.04,1-10/252,0.4893);
[C_INTC,P_INTC] = blsprice(INTC_sim_prices,INTC_K,0.04,(3/4)-10/252,0.2975);
[C_JPM,P_JPM] = blsprice(JPM_sim_prices,JPM_K,0.04,0.5-(10/252),0.2118);
[C_PG,P_PG] = blsprice(PG_sim_prices,PG_K,0.04,(3/4)-10/252,0.1621);

% 6. Simulate Portfolios

PL_INTC = (C_INTC - Call_INTC);
PL_JPM = (P_JPM - Put_JPM);
PL_AA = (C_AA - Call_AA);
PL_PG = (P_PG - Put_PG);

PL_port = -3*PL_INTC + 6*PL_JPM + 6*PL_AA - 2*PL_PG;

% 7. VaR and ES Computation

[VaR_INTC  ES_INTC]= get_risk_sample(PL_INTC, alpha);
[VaR_JPM  ES_JPM]= get_risk_sample(PL_JPM, alpha);
[VaR_AA  ES_AA]= get_risk_sample(PL_AA, alpha);
[VaR_PG  ES_PG]= get_risk_sample(PL_PG, alpha);
[VaR_port  ES_port]= get_risk_sample(PL_port , alpha);

VaRanalysis=table(VaR_INTC, VaR_JPM,VaR_AA,VaR_PG, VaR_port)
ESanalysis=table(ES_INTC, ES_JPM,ES_AA,ES_PG, ES_port)

%-----graphs-----%

h=figure('Color',[1 1 1])
subplot(5,1,1)
histogram(-3*PL_INTC, 'normalization','pdf')
xlabel('Simulated P\&L (INTC)','interpreter','latex')
ylabel('PDF','interpreter','latex')
xlim([min(PL_port) max(PL_port)]) 
subplot(5,1,2)
histogram(6*PL_JPM, 'normalization','pdf')
xlabel('Simulated P\&L (JPM)','interpreter','latex')
ylabel('PDF','interpreter','latex')
xlim([min(PL_port) max(PL_port)])
subplot(5,1,3)
histogram(6*PL_AA, 'normalization','pdf')
xlabel('Simulated P\&L (AA)','interpreter','latex')
ylabel('PDF','interpreter','latex')
xlim([min(PL_port) max(PL_port)])
subplot(5,1,4)
histogram(-2*PL_PG, 'normalization','pdf')
xlabel('Simulated P\&L (PG)','interpreter','latex')
ylabel('PDF','interpreter','latex')
xlim([min(PL_port) max(PL_port)])
subplot(5,1,5)
histogram(PL_port, 'normalization','pdf')
xlabel('Simulated P\&L (portfolio)','interpreter','latex')
ylabel('PDF','interpreter','latex')
xlim([min(PL_port) max(PL_port)]) 

h=figure('Color',[1 1 1])
subplot(2,1,1)
X = categorical({'INTC','JPM','AA','PG','Portfolio'});
bar(X, [VaR_INTC, VaR_JPM, VaR_AA, VaR_PG, VaR_port],0.4)
title('VaR','interpreter','latex')
subplot(2,1,2)
X = categorical({'INTC','JPM','AA','PG','Portfolio'});
bar(X, [ES_INTC, ES_JPM, ES_AA, ES_PG, ES_port],0.4)
title('Expected Shortfall','interpreter','latex')

% 8. Marginal and Component VaR 

eps = VaR_port*0.01;
position = find((PL_port<-VaR_port+eps).*(PL_port>-VaR_port-eps));
MVaR =  -[mean(PL_INTC(position)) mean(mean(PL_JPM(position))) mean(mean(PL_AA(position))) mean(mean(PL_PG(position)))]
CVaR =  [-3*MVaR(1) 6*MVaR(2) 6*MVaR(3) -2*MVaR(4)]
%check
[sum(CVaR) VaR_port]
CVaRp = CVaR/sum(CVaR)

% 9. Compute Marginal Expected Shortfall
position = find((PL_port<-VaR_port));
MES =  -[mean(PL_INTC(position)) mean(mean(PL_JPM(position))) mean(mean(PL_AA(position))) mean(mean(PL_PG(position)))]
CES =  [-3*MES(1) 6*MES(2) 6*MES(3) -2*MES(4)]
[sum(CES) ES_port]
CESp = CES/sum(CES)

h=figure('Color',[1 1 1])
subplot(2,1,1)
X = categorical({'INTC','JPM','AA','PG'});
bar(X, CVaRp,0.4 )
title('Risk Contribution (VaR)','interpreter','latex')
subplot(2,1,2)
X = categorical({'INTC','JPM','AA','PG'});
bar(X, CESp,0.4 )
title('Risk Contribution (ES)','interpreter','latex')

%-----function-----%

function [VaR ES] = get_risk_sample(pandl, alpha)

[Nobs Nseries] = size(pandl);
NConf = length(alpha);

for i=1:NConf
    for j=1:Nseries
    VaR(i,j) = -quantile(pandl(:,j), 1-alpha(i));
    
    ES(i,j) = - mean(pandl(:,j).*(pandl(:,j)<=-VaR(i,j)))/(1-alpha(i)) ...
            + ((1-alpha(i))-mean((pandl(:,j)<=-VaR(i,j))))*VaR(i,j)/(1-alpha(i))  ;  
    end
end
end
