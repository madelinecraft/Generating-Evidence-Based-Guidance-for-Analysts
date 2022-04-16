# Generating Evidence Based Guidance for Analysts


















This simulation has three components:

0% Missing Data
MCAR Simulation with 10% Missing Data
MAR Simulation with 10% Missing Data

R and SAS scripts are available for all three components for the condition where alpha_0 = 0.2, the scale variance = 1.1, and the location-scale covariance = 0.3. This is the condition for which figures were presented in the body of the paper.

The R file stored in the 1st component contains code for simulating the data. The R files stored in the 2nd and 3rd components contain code for simulating the data, amputing the data, and imputing the data.

For the data simulation, researchers can easily make alterations to the population parameter values, number of individuals, number of repeated measures per individual, and number of replications.

For the data amputation, researchers can easily make alterations to the proportion of missing data.

Additional modifications to the data simulation, amputation, and imputation processes will be more difficult (e.g. modifications to the type/number of simulated covariates, relationships between covariates, missing data pattern, imputation model).

The SAS files stored in all three components contain code for fitting the multilevel location scale model to each replication for the condition where the number of individuals is 50 and the number of repeated measures per individual is 10.

Researchers can easily make alterations to the SAS scripts such that the model is being fitted to alternative sample size conditions by editing the file path for the SAS library and the simulated data set they desire to import and fit the model to.
