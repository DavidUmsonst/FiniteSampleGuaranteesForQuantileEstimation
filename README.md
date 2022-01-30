# Finite Sample Guarantees For QuantileEstimation - An Application to Detector Threshold Tuning

In this repository you can find the code used in the article **Finite Sample Guarantees For Quantile Estimation - An Application to Detector Threshold Tuning**. You can find the arXiv version of the article [here](https://arxiv.org/abs/2105.12239).

In this article, we determine how many independent and identically distributed samples of a certain distribution are needed to estimate the distribution's quantile with a certain confidence

This repository provides the code to reproduce the result presented in the paper. The code was originally coded in MATLAB, but we also provide two Julia scripts to reproduce the simulated results.

## How to use the code:
The code can be divided into the following three separate parts:
1. **Comparing the finite guarantees**

    To reproduce Fig. 2 that compares the three finite guarantees we propose in our article one can...
    - **USING MATLAB**
    
      ...run the file `ComparisonOfSampleGuarantees.m`, which uses `FiniteSampleBoundBetaConfInt.m`.  
    - **USING JULIA**
    
      ...run the file  `ComparisonOfSampleGuarantees.jl`.
2. **Creating histograms for different distributions**

    To reproduce Fig. 3 that compares the performance of two guarantees for three different distributions one can...
    - **USING MATLAB**
    
      ...run the file `HistogramOfEmpiricalFarForFiniteGuarantees.m`, which uses the file `FiniteSampleBoundBetaConfInt.m`.  
    - **USING JULIA**
    
      ...run the file  `HistogramOfEmpiricalFarForFiniteGuarantees.jl`.
3. **Applying the results to real data**

    To reproduce Fig. 4, Fig. 5, and Fig. 6 one can run the file `HistogramOfEmpiricalFarForFiniteGuarantees.m`, which uses the real data saved in `DataOfTCLabForDetectorTuning.mat` and the file`FiniteSampleBoundBetaConfInt.m`.  *Note for reproducing these result only MATLAB code is available!*
    
