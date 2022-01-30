function [N_sample] = FiniteSampleBoundBetaConfInt(gamma_quantile,desired_conf,epsilon_desired)
%FINITESAMPLEBOUNDBETACONFINT This function takes desired quantile, p, the 
%desired confidence, desired_conf, the desired deviation from the quantile,
%epsilon_deired and calculates the number of samples needed to estimate the
%quantile given the desired confidence and deviation

if epsilon_desired>min(gamma_quantile,1-gamma_quantile)  % ideally epsilon should be smaller than or equal to min(gamma,1-gamma)
    disp('\epsilon is chosen too large!')
end

if gamma_quantile<0.5
    gamma_quantile=1-gamma_quantile;
end

[~,n2]=rat(gamma_quantile);  
q=1-gamma_quantile;
QuantileOfGaussDist=norminv(1-(1-desired_conf)/2);
a_low=-QuantileOfGaussDist*sqrt(gamma_quantile*q)/(epsilon_desired*sqrt(n2));
b_low=1/(3*epsilon_desired*n2)*(2*(0.5-gamma_quantile)*QuantileOfGaussDist^2-(1+gamma_quantile));    
k1low_sqrt=-a_low/2+sqrt(a_low^2/4-b_low);
k_low=k1low_sqrt^2;
N_sample=ceil(k_low)*n2


end



