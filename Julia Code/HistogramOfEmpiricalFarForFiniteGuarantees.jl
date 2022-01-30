using Random, Distributions, Plots, LaTeXStrings

# Function for determining sample size from confidence interval of beta distribution
function FiniteSampleBoundBetaConfInt(gamma,desired_conf,desired_epsilon)
    if desired_epsilon>min(gamma,1-gamma)      # ideally epsilon should be smaller than or equal to min(gamma,1-gamma)
        println("epsilon is chosen too large!") 
    end

    if gamma<0.5 
        gamma=1-gamma
    end

    fract=rationalize(gamma)
    n2=denominator(fract)
    q=1-gamma
    QuantileOfGaussDist=quantile.(Normal(),(1-(1-desired_conf)/2))
    a_low=-QuantileOfGaussDist*sqrt(gamma*q)/(desired_epsilon*sqrt(n2))
    b_low=1/(3*desired_epsilon*n2)*(2*(0.5-gamma)*QuantileOfGaussDist^2-(1+gamma)) 
    k_sqrt=-a_low/2+sqrt(a_low^2/4-b_low)
    k_low=k_sqrt^2
    return ceil(Int,k_low)*n2
end

# Testing the samples sizes to get certain confidence intervals for
# different distributions
gamma_test=0.95 # quantile we look for
rho_guarantee=0.05 # confidence interval with 1-rho_guarantee
epsilon_diff=0.01 # allowed deviation from desired quantile
fract=rationalize(gamma_test)
den_test=denominator(fract)
N_DKW=ceil(Int,log(2/rho_guarantee)/(2*epsilon_diff^2*den_test))*den_test
N_beta=FiniteSampleBoundBetaConfInt(gamma_test,1-rho_guarantee,epsilon_diff)
n_r=4 # degrees of freedom for chi squared distribution



## Determine thresholds for different sets of samples with same size
N_iter=1000

JD_ApproxChi2DKW=zeros(N_iter)
JD_ApproxChi2Beta=zeros(N_iter)
JD_ApproxCUSUMDKW=zeros(N_iter)
JD_ApproxCUSUMBeta=zeros(N_iter)
JD_ApproxLevyDKW=zeros(N_iter)
JD_ApproxLevyBeta=zeros(N_iter)

yD_cusumDKW=zeros(N_DKW+1)
yD_cusumBeta=zeros(N_beta+1)
delta_cusum=6 # forgetting factor of cusum detector
chi2=Chisq(n_r)
lev=Levy()
println("Starting to calculate thresholds")
for n=1:N_iter
    # draw samples from chi2 dist and determine thresholds
    yD_chi2DKW=rand(chi2,N_DKW) # draw samples of the input
    yD_chi2Beta=yD_chi2DKW[1:N_beta]
    for i=1:N_DKW
        yD_cusumDKW[i+1]=max(0,yD_cusumDKW[i]+yD_chi2DKW[i]-delta_cusum)
    end
    yD_cusumBeta[1:(N_beta+1)]=yD_cusumDKW[1:(N_beta+1)]
    
    # determine approximated threshold for DKW inequality
    
    yD_cusumDKW_Ordered=sort(yD_cusumDKW) # order samples ascending
    
    yD_chi2DKW_Ordered=sort(yD_chi2DKW) # order samples ascending
    
    JD_ApproxChi2DKW[n]=yD_chi2DKW_Ordered[floor(Int,N_DKW*gamma_test)]
        
    JD_ApproxCUSUMDKW[n]=yD_cusumDKW_Ordered[floor(Int,N_DKW*gamma_test)]
    
    
    
    # determine approximated threshold for Beta inequality
    
    yD_cusumBeta_Ordered=sort(yD_cusumBeta) # order samples ascending
    
    yD_chi2Beta_Ordered=sort(yD_chi2Beta) # order samples ascending
   
    JD_ApproxChi2Beta[n]=yD_chi2Beta_Ordered[floor(Int,N_beta*gamma_test)]
        
    JD_ApproxCUSUMBeta[n]=yD_cusumBeta_Ordered[floor(Int,N_beta*gamma_test)]
    
    
    # draw samples from Levy distribution and determine thresholds for DKW
    yD_LevyDKW=rand(lev,N_DKW) # draw samples of the input
    
    yD_LevyDKWordered=sort(yD_LevyDKW) # order samples ascending

    JD_ApproxLevyDKW[n]=yD_LevyDKWordered[floor(Int,N_DKW*gamma_test)]
    
    # Determine threshold for Levy distribution for Beta
    
    yD_LevyBeta=yD_LevyDKW[1:N_beta] # draw samples of the input
    
    yD_LevyBetaordered=sort(yD_LevyBeta) # order samples ascending

    JD_ApproxLevyBeta[n]=yD_LevyBetaordered[floor(Int,N_beta*gamma_test)]
    
end
println("Finished calculating thresholds")


# Evaluate determined thresholds

N_test=1000000;
yDtest_Chi2=rand(chi2,N_test);
yDtest_cusum=zeros(N_test+1)
falseAlarmRate_Chi2DKW=zeros(N_iter)
falseAlarmRate_CUSUMDKW=zeros(N_iter)
falseAlarmRate_Chi2Beta=zeros(N_iter)
falseAlarmRate_CUSUMBeta=zeros(N_iter)
for i=1:N_test
    yDtest_cusum[i+1]=max(0,yDtest_cusum[i]+yDtest_Chi2[i]-delta_cusum)
end

yDtest_Levy=rand(lev,N_test);
falseAlarmRate_LevyDKW=zeros(N_iter);
falseAlarmRate_LevyBeta=zeros(N_iter);

println("Start obtaining empirical false alarm rate")

for j=1:N_iter
        
    falseAlarmRate_Chi2DKW[j]=sum(yDtest_Chi2.>JD_ApproxChi2DKW[j])/N_test
    
    falseAlarmRate_CUSUMDKW[j]=sum(yDtest_cusum.>JD_ApproxCUSUMDKW[j])/N_test
    
    falseAlarmRate_LevyDKW[j]=sum(yDtest_Levy.>JD_ApproxLevyDKW[j])/N_test
    
    
    
    falseAlarmRate_Chi2Beta[j]=sum(yDtest_Chi2.>JD_ApproxChi2Beta[j])/N_test
    
    falseAlarmRate_CUSUMBeta[j]=sum(yDtest_cusum.>JD_ApproxCUSUMBeta[j])/N_test
    
    falseAlarmRate_LevyBeta[j]=sum(yDtest_Levy.>JD_ApproxLevyBeta[j])/N_test
end    

println("Finished obtaining empirical false alarm rates")


println("Percentage of false alarm rates outside of desired interval for chi2 and DKW:")
println(100/N_iter*(sum(falseAlarmRate_Chi2DKW.<1-gamma_test-epsilon_diff)+sum(falseAlarmRate_Chi2DKW.>1-gamma_test+epsilon_diff)))

println("Percentage of false alarm rates outside of desired interval for chi2 and beta")
println(100/N_iter*(sum(falseAlarmRate_Chi2Beta.<1-gamma_test-epsilon_diff)+sum(falseAlarmRate_Chi2Beta.>1-gamma_test+epsilon_diff)))

println("Percentage of false alarm rates outside of desired interval for CUSUM and DKW")
println(100/N_iter*(sum(falseAlarmRate_CUSUMDKW.<1-gamma_test-epsilon_diff)+sum(falseAlarmRate_CUSUMDKW.>1-gamma_test+epsilon_diff)))

println("Percentage of false alarm rates outside of desired interval for CUSUM and beta")
println(100/N_iter*(sum(falseAlarmRate_CUSUMBeta.<1-gamma_test-epsilon_diff)+sum(falseAlarmRate_CUSUMBeta.>1-gamma_test+epsilon_diff)))

println("Percentage of false alarm rates outside of desired interval for Levy and DKW")
println(100/N_iter*(sum(falseAlarmRate_LevyDKW.<1-gamma_test-epsilon_diff)+sum(falseAlarmRate_LevyDKW.>1-gamma_test+epsilon_diff)))

println("Percentage of false alarm rates outside of desired interval for Levy and beta")
println(100/N_iter*(sum(falseAlarmRate_LevyBeta.<1-gamma_test-epsilon_diff)+sum(falseAlarmRate_LevyBeta.>1-gamma_test+epsilon_diff)))


# plot empirical false alarm rate for chi-squared distribution
p1=histogram([falseAlarmRate_Chi2DKW falseAlarmRate_Chi2Beta], fillalpha=[1 0.4],xlims=(0.03, 0.07), label=[L"N_{\mathrm{DKW}}=18460" L"N_{\mathrm{beta}}=2180"],bins=[:sturges 50])
vline!([0.04 0.06],linestyle= [:dashdot :dashdot], color =[:black :black],label=nothing)

# plot empirical false alarm rate for Levy distribution
p2=histogram([falseAlarmRate_LevyDKW falseAlarmRate_LevyBeta],fillalpha=[1 0.4],xlims=(0.03, 0.07),legend=false,bins=[:sturges 50])
vline!([0.04 0.06],linestyle= [:dashdot :dashdot], color =[:black :black],legend=false)

# plot empirical false alarm rate for distribution from CUSUM detector
p3=histogram([falseAlarmRate_CUSUMDKW falseAlarmRate_CUSUMBeta],xlims=(0.03, 0.07),fillalpha=[1 0.4],legend=false,bins=[25 75])
vline!([0.04 0.06],linestyle= [:dashdot :dashdot], color =[:black :black],legend=false)

#plot histogram of empirical false alarm rate
plot(p1,p2,p3,layout=(3,1))
