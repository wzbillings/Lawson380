# R Codes for MATH 380 (Math Research) with Dr. Lawson; Fall 2019.

Project goal: develop least squares methods to obtain parameters from data which is best modelled by the logistic ODE.

Description of contents:

## Scripts
Contains five major scripts with the bulk of the R code for the project.
* `Dimensionaless_exploration.R`: provides code for generating dimensionless logistic growth time-series data, and estimating the growth parameter from the data. Unsure if this should remain independent or not.
* `Helpers.R`: provides functions that needed to be in multiple places at one time. Currently, this is `validate_inputs()` and `prep_data()`. The other scripts likely will not run without sourcing this file first.
* `Least_squares_methods.R`: provides various methods for fitting models to get the parameters back from time series logistic growth data.
* `Old_stuff.R`: functions and stuff that were necessary at one point but have now been put here because they were taking up space and not doing much.
* `Sample_data_generation.R`: provides several methods for generating logistic growth time series data in multiple ways. Needs to be cleaned up as many of the functions could be condensed into a few functions with some default arguments.

## Tests
Contains R-markdown documents and R-scripts used for monkeying around.

## Manuscript
This folder will eventually contain everything needed for my final manuscript for this project.

## Presentation
Yep, this will also exist soon.

