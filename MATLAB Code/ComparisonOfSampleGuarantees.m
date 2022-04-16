gamma_iter=0.01:0.01:0.99;
n_gamma=length(gamma_iter);
gamma_test=0.95; % quantile we look for
rho_guarantee=0.05; % confidence interval with 1-rho_guarantee
epsilon_diff=0.01; % allowed deviation from desired quantile
N_beta=zeros(1,n_gamma);
N_DKW=zeros(1,n_gamma);
N_VP=zeros(1,n_gamma);

% calculate the sample guarantees for several values of gamma
for i=1:n_gamma
    [~,den_gamma]=rat(gamma_iter(i));
    N_beta(i)=FiniteSampleBoundBetaConfInt(gamma_iter(i),1-rho_guarantee,epsilon_diff);
    N_DKW(i)=ceil(log(2/rho_guarantee)/(2*epsilon_diff^2*den_gamma))*den_gamma;
    N_VP(i)=ceil(1/den_gamma*((4/9)*gamma_iter(i)*(1-gamma_iter(i))/(epsilon_diff^2*rho_guarantee)-1))*den_gamma-1;
end
% plot finite sample guarantees
figure
plot(gamma_iter,N_DKW,'-.',gamma_iter,N_VP,'r--',gamma_iter,N_beta,'k')
xlabel('$\gamma$','Interpreter','latex')
ylabel('Sample size','Interpreter','latex')
legend('$N_{\mathrm{DKW}}$','$N_{\mathrm{VP}}$','$N_{\mathrm{beta}}$','Interpreter','latex')

