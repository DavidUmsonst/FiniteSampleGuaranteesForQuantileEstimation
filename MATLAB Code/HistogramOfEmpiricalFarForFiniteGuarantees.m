% Testing the samples sizes to get certain confidence intervals for
% different distributions
gamma_test=0.95; % quantile we look for
rho_guarantee=0.05; % confidence interval with 1-rho_guarantee
epsilon_diff=0.01; % allowed deviation from desired quantile
[~,den_test]=rat(gamma_test);
N_DKW=ceil(log(2/rho_guarantee)/(2*epsilon_diff^2*den_test))*den_test;
N_Beta=FiniteSampleBoundBinomConfInt1(gamma_test,1-rho_guarantee,epsilon_diff);
n_r=4; % degrees of freedom for chi squared distribution
% obtain probability distirbution of Levy distribution
pd_levy=makedist('Stable','alpha',0.5,'beta',1,'gam',1,'delta',0);

%% Determine thresholds for different sets of samples with same size
N_iter=1000;

JD_ApproxChi2DKW=zeros(1,N_iter);
JD_ApproxCUSUMDKW=zeros(1,N_iter);
JD_ApproxLevyDKW=zeros(1,N_iter);
JD_ApproxChi2Beta=zeros(1,N_iter);
JD_ApproxCUSUMBeta=zeros(1,N_iter);
JD_ApproxLevyBeta=zeros(1,N_iter);

yD_cusumDKW=zeros(1,N_DKW+1);
yD_cusumBeta=zeros(1,N_Beta+1);
delta_cusum=6;

disp('Starting to calculate thresholds')
for n=1:N_iter
    % draw samples from chi2 dist and determine thresholds
    yD_chi2DKW=chi2rnd(n_r,1,N_DKW); % draw samples of the input
    yD_chi2Binom=yD_chi2DKW(1:N_Beta);
    for i=1:N_DKW
       yD_cusumDKW(i+1)=max(0,yD_cusumDKW(i)+yD_chi2DKW(i)-delta_cusum);
    end
    yD_cusumBeta=yD_cusumDKW(1:N_Beta+1);
    
    % determine approximated threshold for DKW inequality
    
    yD_CusumDKWordered=sort(yD_cusumDKW); % order samples ascending
    
    yD_ChiDKWordered=sort(yD_chi2DKW); % order samples ascending

    JD_ApproxChi2DKW(n)=yD_ChiDKWordered(N_DKW*gamma_test);
        
    JD_ApproxCUSUMDKW(n)=yD_CusumDKWordered(N_DKW*gamma_test);
    
    
    
    % determine approximated threshold for Binom inequality
    
    yD_CusumBinomordered=sort(yD_cusumBeta); % order samples ascending
    
    yD_ChiBinomordered=sort(yD_chi2Binom); % order samples ascending

    JD_ApproxChi2Beta(n)=yD_ChiBinomordered(floor(N_Beta*gamma_test));
        
    JD_ApproxCUSUMBeta(n)=yD_CusumBinomordered(floor(N_Beta*gamma_test));

    
    % draw samples from Levy dist and determine thresholds for DKW
    yD_LevyDKW=random(pd_levy,1,N_DKW); % draw samples of the input
    
    yD_LevyDKWordered=sort(yD_LevyDKW); % order samples ascending

    JD_ApproxLevyDKW(n)=yD_LevyDKWordered(N_DKW*gamma_test);
    
    % Determine threshold for Levy for Beta
    
    yD_LevyBinom=yD_LevyDKW(1:N_Beta); % draw samples of the input
    
    yD_LevyBinomordered=sort(yD_LevyBinom); % order samples ascending

    JD_ApproxLevyBeta(n)=yD_LevyBinomordered(floor(N_Beta*gamma_test));
end
println('Finished calculating thresholds')

%% Evaluate determined thresholds

N_test=1000000;
yDtest_Chi2=chi2rnd(n_r,1,N_test);
yDtest_cusum=zeros(1,N_test+1);
falseAlarmRate_Chi2DKW=zeros(1,N_iter);
falseAlarmRate_CUSUMDKW=zeros(1,N_iter);
falseAlarmRate_Chi2Beta=zeros(1,N_iter);
falseAlarmRate_CUSUMBeta=zeros(1,N_iter);
for i=1:N_test
       yDtest_cusum(i+1)=max(0,yDtest_cusum(i)+yDtest_Chi2(i)-delta_cusum);
end

yDtest_Levy=random(pd_levy,1,N_test);
falseAlarmRate_LevyDKW=zeros(1,N_iter);
falseAlarmRate_LevyBeta=zeros(1,N_iter);

disp('Start obtaining empirical false alarm rate')
for j=1:N_iter
    falseAlarmRate_Chi2DKW(j)=sum(yDtest_Chi2>JD_ApproxChi2DKW(j))/N_test;
    
    falseAlarmRate_CUSUMDKW(j)=sum(yDtest_cusum>JD_ApproxCUSUMDKW(j))/N_test;
    
    falseAlarmRate_LevyDKW(j)=sum(yDtest_Levy>JD_ApproxLevyDKW(j))/N_test;
     
    
    
    falseAlarmRate_Chi2Beta(j)=sum(yDtest_Chi2>JD_ApproxChi2Beta(j))/N_test;
    
    falseAlarmRate_CUSUMBeta(j)=sum(yDtest_cusum>JD_ApproxCUSUMBeta(j))/N_test;
    
    falseAlarmRate_LevyBeta(j)=sum(yDtest_Levy>JD_ApproxLevyBeta(j))/N_test;
end

disp('Finished obtaining empirical false alarm rates')
%%
disp('Percentage of false alarm rates outside of desired interval for chi2 and DKW')
100/N_iter*(sum(falseAlarmRate_Chi2DKW<1-gamma_test-epsilon_diff)+sum(falseAlarmRate_Chi2DKW>1-gamma_test+epsilon_diff))

disp('Percentage of false alarm rates outside of desired interval for chi2 and Beta')
100/N_iter*(sum(falseAlarmRate_Chi2Beta<1-gamma_test-epsilon_diff)+sum(falseAlarmRate_Chi2Beta>1-gamma_test+epsilon_diff))

disp('Percentage of false alarm rates outside of desired interval for CUSUM and DKW')
100/N_iter*(sum(falseAlarmRate_CUSUMDKW<1-gamma_test-epsilon_diff)+sum(falseAlarmRate_CUSUMDKW>1-gamma_test+epsilon_diff))

disp('Percentage of false alarm rates outside of desired interval for CUSUM and Beta')
100/N_iter*(sum(falseAlarmRate_CUSUMBeta<1-gamma_test-epsilon_diff)+sum(falseAlarmRate_CUSUMBeta>1-gamma_test+epsilon_diff))

disp('Percentage of false alarm rates outside of desired interval for Levy and DKW')
100/N_iter*(sum(falseAlarmRate_LevyDKW<1-gamma_test-epsilon_diff)+sum(falseAlarmRate_LevyDKW>1-gamma_test+epsilon_diff))

disp('Percentage of false alarm rates outside of desired interval for Levy and Beta')
100/N_iter*(sum(falseAlarmRate_LevyBeta<1-gamma_test-epsilon_diff)+sum(falseAlarmRate_LevyBeta>1-gamma_test+epsilon_diff))
%%
figure
subplot(3,1,1)
histogram(falseAlarmRate_Chi2DKW,0.03:0.001:0.07);
hold on
histogram(falseAlarmRate_Chi2Beta,0.03:0.001:0.07);
legend({['$N_{\mathrm{DKW}}=' num2str(N_DKW),'$'],['$N_{\mathrm{beta}}=', num2str(N_Beta), '$']},'AutoUpdate','off','Interpreter','Latex')
xline(1-gamma_test-epsilon_diff,'k-.');
xline(1-gamma_test+epsilon_diff,'k-.');
title('Histogram of the empirical false alarm rates for $\chi^2(4)$ distribution','Interpreter','Latex')
subplot(3,1,2)
histogram(falseAlarmRate_LevyDKW,0.03:0.001:0.07);
hold on
histogram(falseAlarmRate_LevyBeta,0.03:0.001:0.07);
xline(1-gamma_test-epsilon_diff,'k-.');
xline(1-gamma_test+epsilon_diff,'k-.');
title('Histogram of the empirical false alarm rates for L\''evy distribution','Interpreter','Latex')
subplot(3,1,3)
histogram(falseAlarmRate_CUSUMDKW,0.03:0.001:0.07);
xlabel('Empirical false alarm rate','Interpreter','Latex')
hold on
histogram(falseAlarmRate_CUSUMBeta,0.03:0.001:0.07);
xline(1-gamma_test-epsilon_diff,'k-.');
xline(1-gamma_test+epsilon_diff,'k-.');
title('Histogram of the empirical false alarm rates for CUSUM detector','Interpreter','Latex')
