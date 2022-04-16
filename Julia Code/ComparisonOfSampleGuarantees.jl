
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
    return ceil(k_low)*n2
end

gamma_iter=0.01:0.01:0.99
n_gamma=length(gamma_iter)
rho_guarantee=0.05 # confidence interval with 1-rho_guarantee
epsilon_diff=0.01 # allowed deviation from desired quantile
N_beta=zeros(n_gamma)
N_DKW=zeros(n_gamma)
N_VP=zeros(n_gamma)

# calculate the sample guarantees for several values of gamma
for i=1:n_gamma
    irreducibleFraction=rationalize(gamma_iter[i])
    den_gamma=denominator(irreducibleFraction)
    N_beta[i]=FiniteSampleBoundBetaConfInt(gamma_iter[i],1-rho_guarantee,epsilon_diff)
    N_DKW[i]=ceil(log(2/rho_guarantee)/(2*epsilon_diff^2*den_gamma))*den_gamma
    N_VP[i]=ceil((1/den_gamma)*(4/9*gamma_iter[i]*(1-gamma_iter[i])/(epsilon_diff^2*rho_guarantee)-1))*den_gamma-1
end
# plot finite sample guarantees
pyplot()
plot(gamma_iter,[N_beta,N_DKW,N_VP],linestyle= [:solid :dashdot :dash], color =[:black :blue :red], label=[L"N_{\mathrm{beta}}" L"N_{\mathrm{DKW}}" L"N_{\mathrm{VP}}"] )#,:black)#,label="N_{\mathrm{beta}}")
xlabel!(L"$\gamma$")
ylabel!("Sample size")
