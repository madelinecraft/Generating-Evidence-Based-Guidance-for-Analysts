# Generating Evidence Based Guidance for Analysts

## Summary:

The goal of this project was to conduct a study explicating the conditions under which a statistical analysis of interest produced accurate results. The study demonstrated that the statistical analysis of interest would lead to inaccurate conclusions as currently implemented. 

## Application of Results:

We wrote a report detailing our study and providing evidence-based guidance to analysts wishing to implement the statistical analysis of interest. 

## Details of the Analysis:

The statistical analysis of interest is a missing data technique for augmenting missing values in such a way that all important features of the dataset are retained. To evaluate the performance of this technique, we conducted a Monte Carlo data simulation, induced missing values, and applied the missing data technique. Since we simulated the data and induced the missing values ourselves, we could compare characteristics of the augmented values to characteristics of the simulated values, thus evaluating performance of the technique. 

A more extensive compilation of the R and SAS scripts written for this project is available at https://osf.io/uqdg5/?view_only=4bb97ffdd542422a9513ee5d8100a60f

The written report is stored above as a Word document ("MI for MLSM Manuscript 2.0.docx").

A small compilation of the R and SAS scripts written for this project is stored above. 

* An R script for simulating the data, induced the missing values, and applying the missing data technique is stored above for one of the simulation conditions ("simampimpUNIX.R). The simulation condition presented is the one for which figures were presented in the body of the written report.

* Alterations to the conditions of the data simulation (population parameter values, number of individuals, number of repeated measures per individual, and number of replications) can be made. Alterations to the proportion of missing data can also be made.

* Additional modifications will be more difficult (e.g. modifications to the type/number of simulated covariates, relationships between covariates, missing data pattern, imputation model).

* A SAS script for fitting the statistical model of interest to the simulated, augmented data for component (3) is stored above for the condition where the number of indivdiuals is 50 and the number of repeated measures per individual is 10 ("s10_50.sas"). 

* Alterations to the SAS scripts can be made such that the statistical model of interest is being fitted to alternative sample size conditions by editing the file path for the SAS library and the simulated data set they desire to import and fit the model to.
