# Generating Evidence Based Guidance for Analysts

## Summary:

The goal of this project was to conduct a study explicating the conditions under which a statistical analysis of interest produced accurate results. The study demonstrated that the statistical analysis of interest would lead to inaccurate conclusions as currently implemented. 

## Application of Results:

We wrote a report detailing our study and providing evidence-based guidance to analysts wishing to implement the statistical analysis of interest in the future. 

## Details of the Analysis:

The statistical analysis of interest is a missing data technique for augmenting missing values in such a way that all important features of the dataset are retained. To evaluate the performance of this technique, we conducted a Monte Carlo data simulation, induced missing values, and applied the missing data technique. Since we simulated the data and induced the missing values ourselves, we could compare characteristics of the augmented values to characteristics of the simulated values, thus evaluating performance of the technique. 

This simulation has three components:

(1) 0% Missing Data
(2) MCAR Simulation with 10% Missing Data
(3) MAR Simulation with 10% Missing Data

* All components are available at https://osf.io/uqdg5/?view_only=4bb97ffdd542422a9513ee5d8100a60f

* R and SAS scripts are for (3)

As an example, we provide R and SAS scripts for (3)
* 
* 
* R and SAS scripts are available for all three components for the condition where alpha_0 = 0.2, the scale variance = 1.1, and the location-scale covariance = 0.3 (sim. This is the condition for which figures were presented in the body of the paper.

The R file stored in the 1st component contains code for simulating the data. The R files stored in the 2nd and 3rd components contain code for simulating the data, amputing the data, and imputing the data.

For the data simulation, researchers can easily make alterations to the population parameter values, number of individuals, number of repeated measures per individual, and number of replications.

For the data amputation, researchers can easily make alterations to the proportion of missing data.

Additional modifications to the data simulation, amputation, and imputation processes will be more difficult (e.g. modifications to the type/number of simulated covariates, relationships between covariates, missing data pattern, imputation model).

The SAS files stored in all three components contain code for fitting the multilevel location scale model to each replication for the condition where the number of individuals is 50 and the number of repeated measures per individual is 10.

Researchers can easily make alterations to the SAS scripts such that the model is being fitted to alternative sample size conditions by editing the file path for the SAS library and the simulated data set they desire to import and fit the model to.
