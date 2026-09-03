# Spence et al. Section 4 NumPyro transcription

Python 3.13.12 

This bundle implements the published `plankton_model.stan` posterior structure:

- 20 parallel chains, 5,000 warm-up draws, and 5,000 retained draws by default
- each row of `P` has a `Dirichlet(1)` prior
- raw validation rows follow `Multinomial(P_i)`
- `alpha_b` has the uniform simplex prior, represented as `Dirichlet(1)`
- `b_par` and `k_nb` use `HalfNormal(1000)`, equivalent to the published positive-truncated `Normal(0, 1000)` priors
- `log_mu` uses `Normal(0, 1000)` and `mu = exp(log_mu)`
- each field prevalence `q_l` uses `Dirichlet(b_par * alpha_b)`
- field class counts use `Multinomial(q_l @ P)`
- observed totals use `NegativeBinomial2(mu, k_nb)`

The reprex in this branch is built around a timeseries of predictions data that resembles what we might expect to collect, but the model has fitted poorly to it, so at present should be considered unfinished.
The model, model priors and even architecture does not appear to be sensible for the present dataset. Our validation dataframe / confusion matrix is extremely small, and a dirchlet of 1 is able to inject a disproportionately large sum of confusion into the confusion matrix.
To my best understanding the full original model is replicated, the field series and validation counts jointly update `P`. Timeseries plots are post-processing layers and do not alter the fitted posterior.

## Setup (miniforge prompt instructions for WINDOWS). First clone the repo and cd into it:
```bash
cd %USERPROFILE%\Documents
git clone --branch pythonavoidingconfusion https://github.com/CefasRepRes/AvoidingConfusion.git
cd AvoidingConfusion
conda create -n pythonavoidingconfusion python=3.13
conda activate pythonavoidingconfusion
pip install -r requirements.txt
```

## Run 
Plot_class_counts_timeseries.py is a wrapper around the bayesian functions and should be the entry point.

```bash
python plot_class_counts_timeseries.py --input-json timeseries.json --validation-json validation.json --target-class fish_larvae --mc-samples 5000 --mc-seed 42 --outdir summary_timeseries_out
```
