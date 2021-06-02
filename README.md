# Mammal-responses-logging-concessions

Data and codes to run the Bayesian Multi-species Occupancy Models (MSOM) used in [Carvalho Jr et al. 2021. Mammal responses to reduced-impact logging in Amazonian forest concessions. *Forest Ecolgy and Management*. 496, 119401](https://doi.org/10.1016/j.foreco.2021.119401)

### Description

The *data_jamari.rds* file contains the data for the analysis. This include camera-trap data (species detection histories) and site and sampling covariates in the format required by the analysis code (model 1 and 2).

*model1.r* and *model2.r* have the R code to run the Bayesian multi-species occupancy model. Model 1 is for avaluating logging status (unlogged or logged), and model 2 for evaluation of logging intensity and road density effects on mammal community.

Codes were adapted from original codes written by Elise Zipkin and available at the Github page of the [Zipkin Quantitative Ecology Lab](https://github.com/zipkinlab/Community_model_examples-covariate_model/blob/master/covariate%20model%20code.r), with some chunks taken from supplementary material from [Yamaura et al. 2011.](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2664.2010.01922.x)


# Contact Us
If you have any questions please contact <elildojr@gmail.com>
