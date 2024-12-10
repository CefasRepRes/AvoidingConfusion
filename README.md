# AvoidingConfusion

R and stan scripts to produce the outputs and figures of Spence *et al*. (2024) "Avoiding Confusion: Modelling Image Identification Surveys with Classification Errors", *Methods in Ecology and Evolution*.

`ML functions.R` contains the main functions required for the other scripts to run.<br>
`Poisson generating process (3.1).R` runs the Poisson example in Section 3.1.<br>
`Poisson and NB generating processes (3.2).R` runs the mixed Poisson - Negative Binomial example in Section 3.2.<br>
`Sensitivity_to_C (3.3).R` contains the code for the comparison study in Section 3.3.<br>
`Zooplankton survey (4).R` runs the analysis of the Zooplankton survey data in Section 4. This file also contains the code to produce Figures S7 and S8 in Section S5 of the Supporting Information.<br>
`Simulation-based calibration.R` runs the simulation-based calibration in Supporting Information Section S3.<br>
`Zero-inflated Poisson with covariates.R` runs the additional example with a covariate influencing the count generating process, see Section S6 of the Supporting Information.

`poisson_model.stan` contains the stan model used in `Poisson generating process (3.1).R` and `Simulation-based calibration.R`.<br>
`plankton_model.stan` contains the stan model used in `Zooplankton survey (4).R`.<br>
`zeroinf_poisson_model.stan` contains the stan model used in `Zero-inflated Poisson with covariates.R`.

`plank.csv` is the data used by `Zooplankton survey (4).R`.

## Requirements

R version 4.3.2 or higher is required. The following R packages need to be installed to run the R scripts:

```R
install.packages( c("rstan","ggplot2","reshape2",
                    "posterior","bayesplot"
                    "patchwork","MCMCprecision",
                    "VGAM","parallel", "MASS"))
```

## Citation
Spence M., Barry J., Bartos T., Blackwell R., Scott J., Pitois S. 2024. “Avoiding Confusion: Modelling Image Identification Surveys with Classification Errors.” *Methods in Ecology and Evolution*. In press.
