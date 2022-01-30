% In this script we utilize Repeated Training/Test Split as a
% crossvalidation method to estimate the threshold for a CUSUM detector
% without resetting, where the residual data comes from measurements of a
% TCLab process

load('DataOfTCLabForDetectorTuning.mat') % load TCLab data
y=TempMeasured; % extract temperature measurements from TCLab
deltaY=y-Tss+273.15; % determine variation around desired steady state
deltaYHat=[0     0     1     0
           0     0     0     1]*deltaX_controller; % extract estimates of output variation from controller state
residual=deltaY-deltaYHat; % determine residual signals


N_start=781; % set time for when the TCLab has reached its steady state
r_val=residual(:,N_start+1:end); % extract residuals used for mean, covariance, and threshold estimation
N_muSig=1000; % set how how many residuals are used for the estimation of the mean and covariance
r_muSig=r_val(:,1:N_muSig); % extract residual which are used for estimating the mean and covariance

mu_r=mean(r_muSig,2); % estimate mean
Sigma_r=cov(r_muSig'); % estimate covariance
Sigma_r_sqrt=sqrtm(Sigma_r); % calculate square root of covariance matrix for normalization

r_Normalized=Sigma_r_sqrt\(r_val(:,N_muSig+1:end)-mu_r); % normalize the remaining residuals with the estimated mean and covariance
N_normalized=length(r_Normalized(1,:)); % number of normalized residuals available for testing and training


N_T=10000; % number of threshold estimates per gamma
rho_guarantee=0.05; % confidence interval with 1-rho_guarantee
epsilon_diff=0.01; % allowed deviation from desired quantile
gamma_test=0.95:0.005:0.99; % quantiles we want to investigate
N_gamma=length(gamma_test);


delta_cusum=3; % forgetting factor of the CUSUM detector

yD_cusum=zeros(1,N_normalized+1); % initialize CUSUM detector output
N_cusum=length(yD_cusum)-1;

% calculate CUSUM detector output based on the normalized residuals
for i=1:N_normalized
   yD_cusum(i+1)=max(0,yD_cusum(i)+r_Normalized(:,i)'*r_Normalized(:,i)-delta_cusum);
end

% Plot autocorrelation of the CUSUM detector output
figure
autocorr(yD_cusum,'NumLags',100)
ylabel('Sample auto correlation','Interpreter','Latex')
xlabel('Lag','Interpreter','Latex')
title('Sample autocorrelation of $\mathcal{D}$','Interpreter','Latex')

% initialize empirical false alarm rate matrix
EmpFAR_cusum=zeros(N_T,N_gamma);

% estimate thresholds and calculate empirical false alarm rates with
% Repeated Training/Test Spli
for i=1:N_gamma
    % determine samples needed for given gamma
    N_binom=FiniteSampleBoundBinomConfInt2(gamma_test(i),1-rho_guarantee,epsilon_diff);
    N_test=N_cusum+1-N_binom;
    for j=1:N_T
        %randomize CUSUM outpurs
        index_rand=randperm(N_cusum+1);
        yD_rand=yD_cusum(index_rand);
        % take the first N_binom of the randomized outpurs for estimating
        % the threshold
        yD_train=yD_rand(1:N_binom);
        yD_sort=sort(yD_train);
        JD_approx=yD_sort(floor(N_binom*gamma_test(i)));
        % take the remaining outputs for testing the threshold
        yD_test=yD_rand(N_binom+1:end);
        % calculate empirical false alam rate from test set with
        % approximated threshold
        EmpFAR_cusum(j,i)=sum(yD_test>JD_approx)/N_test;
    end
end

% determine the percentage of empirical false alarm rates outside of the
% desired interval
PercOutsideOfConvInterval=zeros(1,N_gamma);
for i=1:N_gamma
    PercOutsideOfConvInterval(i)=100/N_T*(sum(EmpFAR_cusum(:,i)<1-gamma_test(i)-epsilon_diff)+sum(EmpFAR_cusum(:,i)>1-gamma_test(i)+epsilon_diff));
end
disp('Maximum percentage outside of the desired confidence interval for the investigated false alarm rates is')
max(PercOutsideOfConvInterval)


% plot box plot
EmpFAR_cusum_flipped=flip(EmpFAR_cusum,2);
AccFAR=flip(1-gamma_test);
figure
boxplot(EmpFAR_cusum_flipped,AccFAR, 'positions', AccFAR)
xlabel('Acceptable false alarm rate $1-\gamma$','Interpreter','Latex')
ylabel('Empirical false alarm rate','Interpreter','Latex')
hold on
x_shade=[AccFAR(1),AccFAR(end),AccFAR(end),AccFAR(1)];
y_shade=[AccFAR(1)-epsilon_diff,AccFAR(end)-epsilon_diff,AccFAR(end)+epsilon_diff,AccFAR(1)+epsilon_diff];
fill(x_shade,y_shade,'k-.','FaceAlpha',0.1)
xlabel('Acceptable false alarm rate $1-\gamma$','Interpreter','Latex')
ylabel('Empirical false alarm rate','Interpreter','Latex')
legend('$\pm\epsilon$-band','Interpreter','Latex')
 

%% Plot histogram for empirical false alarm rate for gamma=0.95
figure
histogram(EmpFAR_cusum(:,1),0.03:0.001:0.07);
xlabel('Empirical false alarm rate','Interpreter','Latex')
hold on
xline(1-gamma_test(1)-epsilon_diff,'k-.');
xline(1-gamma_test(1)+epsilon_diff,'k-.');
title('Histogram of the empirical false alarm rates for CUSUM detector','Interpreter','Latex')
