source("ML functions.R")

library(rstan)
library(MCMCprecision)
library(ggplot2)
library(reshape2)
library(posterior)

rstan_options(auto_write = T)
ncores <- parallel::detectCores() - 1  #Set this to less than 4 if there are less than 4 cores available on machine
options(mc.cores = ncores)

set.seed(5)

##Simulation-based calibration function (used in loop further down)

mcmc_calib <- function(mu, P, C, l, st_model, nsamps, thin_period, seed) {
  
  P_real <- P
  mu_real <- mu
  #Simulate data according to above parameters
  distribution <- rep("Poisson", 3)
  Xobs <- observedW.gen.multi(n = l, mu = mu_real, pmat = P_real,
                              distribution = distribution)
  xmeans <- colMeans(Xobs)
  
  alpha_P <- matrix(1, nrow = 3, ncol = 3) #Specify Dirichlet hyperparameters for each category
  for (i in 1:3) {
    alpha_P[i,i] <- 15
  }
  
  w_dat <- list(L = l, M = 3, w = Xobs, C = C, alpha_P = alpha_P)
  mcmc.fit <- sampling(st_model, data = w_dat, iter = nsamps, thin = thin_period, seed = seed) #Fit stan model
  arr.fit <- as.array(mcmc.fit)
  Rh <- max(apply(arr.fit, 3, Rhat))
  ess <- min(apply(arr.fit, c(2, 3), ess_basic))
  posterior <- as.matrix(mcmc.fit)
  
  return(list(samples = posterior, Rh = Rh, ess = ess, xmeans = xmeans))
  
}

#Setting initial variables and loading model

nsites <- 15  #Number of locations
s <- c(50, 200, 100) #Totals of each category used to generate confusion matrix

st_model <- stan_model(file = "poisson_model.stan")

#Estimating thinning period to achieve desired effective sample size

target_ess <- 2000
actual_ess <- 0
nIter <- 500

mu <- rgamma(3, shape = 1, rate = 0.001)
nmu <- length(mu)

alpha_P <- matrix(1, nrow = nmu, ncol = nmu)
for (i in 1:nmu) {
  alpha_P[i,i] <- nsites
}
P <- t(sapply(1:nmu, function(i) {rdirichlet(1, alpha_P[i,])}))
C <- Cgen(s, P)

while (actual_ess < target_ess) {
  
  posterior <- mcmc_calib(mu, P, C, nsites, st_model, nIter, thin = 1, seed = 2)$samples
  npars <- dim(posterior)[2]
  actual_ess <- min(sapply(1:npars, function(i) {ess_basic(posterior[,i])}))
  nIter <- nIter*2 #Double number of iterations if effective sample size is too low
  
}

nsamps <- nIter #Number of non-burn-in samples to draw during each MCMC run

#Initialising variables for main calibration routine
#Users can ignore WARNING bits about effective sample sizes - because
#the calibration doesn't actually happen until the main loop below.

## Used nsamples = 1000 for actual run 
nsamples <- 100  #Number of prior draws to make (also the number of independent MCMC runs)
thin <- ceiling(nsamps*2 / target_ess) #Keep 1 out of every 'thin' samples - multiplying nsamps by 2 in expression due to burn-in taking up the first half of all MCMC samples (stan default)
draws <- data.frame(matrix(0, nrow = nsamples, ncol = 3)) #Data frame to store rank statistics of category mean parameters (the mus) calculated from each posterior sample
names(draws) <- c("mu[1]", "mu[2]", "mu[3]")

P_draws <- data.frame(matrix(0, nrow = nsamples, ncol = 9)) #Data frame to store rank statistics of P entries calculated from each posterior sample
P_labs <- expand.grid(1:3, 1:3)
names(P_draws) <- sapply(1:9, function(x) {paste0("P[", P_labs[x,1],",", P_labs[x,2],"]")})

emp_quants <- data.frame(matrix(0, nrow = nsamples, ncol = 3)) #Data frame to store the empirical quantiles of the column means of the simulated data
names(emp_quants) <- c("mu[1]", "mu[2]", "mu[3]")

nchains <- 4 #Number of chains per MCMC run
target_chain_ess <- target_ess/nchains #Target per-chain effective sample size
stan_seed <- 5

#Main simulation-based calibration loop
# This takes a while to run

for (sample in 1:nsamples) {
  #Draw parameter values from priors
  mu <- rgamma(3, shape = 1, rate = 0.001)
  P <- t(sapply(1:3, function(i) {rdirichlet(1, alpha_P[i,])}))
  C <- Cgen(s, P)
  
  #Run calibration 
  n_mixed <- 0
  while (n_mixed < 1) {
    fit_obj <- mcmc_calib(mu, P, C, nsites, st_model, nsamps, thin, stan_seed)
    if (fit_obj$Rh < 1.05) {
      n_mixed = n_mixed + 1 #Keep sample only if Rhat is acceptable
    }
  }
  ess <- fit_obj$ess
  #Run for more samples (nsamps_new) if effective sample size is too low  
  if (ess < target_chain_ess) {
    thin_new <- ceiling(nsamps*2/ess)
    nsamps_new <- target_chain_ess*thin_new
    while (n_mixed < 1) {
      fit_obj <- mcmc_calib(mu, P, C, nsites, st_model, nsamps_new, thin_new, stan_seed)
      if (fit_obj$Rh < 1.05) { #Keep sample only if Rhat is acceptable
        n_mixed = n_mixed + 1
      }
    }
  }
  posterior <- fit_obj$samples
  wmeans <- fit_obj$xmeans
  #Calculate rank statistics for the mus and update draws
  draws[sample,] <- sapply(1:3, function(x){sum(posterior[,x] < mu[x])})
  #Calculate rank statistics for P entries and update P_draws
  P_draws[sample,] <- sapply(1:9, function(x){sum(posterior[,x+3] < P[as.numeric(P_labs[x,1]),as.numeric(P_labs[x,2])])})
  #Calculate quantile of w mean
  emp_quants[sample,] <- sapply(1:3, function(i) {mean(wmeans[i] < posterior[,i])})
}

##Quantiles of uncorrected counts in mean posteriors

#Code to create histogram of empirical quantiles of the w means in the posterior draws of the mus

quantpltdat <- melt(emp_quants)
names(quantpltdat)[2] <- "Quantile"

# Figure S2 in supporting information
emp_quant_plt <- ggplot(quantpltdat) + geom_histogram(aes(x = Quantile), binwidth = 0.05, color = "#005b96", fill = '#6497b1', linewidth = 0.25) + facet_wrap(vars(variable), nrow = 1, labeller = label_parsed) + ylab("Count")
#jpeg("empquantplot1.jpg", height=10, width=30, res=600, units="cm")
emp_quant_plt
#dev.off()

##Diagnosing uniformity

#Code to create empirical cdf plot for ranks of posterior draws of the mus

## Ignore R comment
draws_stacked <- melt(draws)
names(draws_stacked)[2] <- "Rank"
npars <- ncol(draws)
draws_prob <- data.frame(matrix(0, nrow = target_ess*npars, ncol = 2))
draws_prob[,1] <- melt(sapply(1:npars, function(x) {rep(paste0("mu[",x,"]"), target_ess)}))[,3]
draws_prob[,2] <- melt(sapply(1:npars, function(x) {sapply(1:target_ess, function(y) {length(which(draws_stacked[((x-1)*nsamples + 1):(x*nsamples),2] <= y))/nsamples})}))[,3]
true_cdf <- melt(replicate(npars, sapply(1:target_ess, function(x) {x/target_ess})))[,3]
upper_ecdf <- sapply(true_cdf, function(x) {qbinom(0.95, nsamples, x)})/nsamples
lower_ecdf <- sapply(true_cdf, function(x) {qbinom(0.05, nsamples, x)})/nsamples
row_inds <- melt(replicate(npars, seq(1, target_ess)))[,3]
draws_with_actual <- cbind(draws_prob, true_cdf, upper_ecdf, lower_ecdf, row_inds)
names(draws_with_actual) <- c("variable", "Value", "true", "upper", "lower", "Rank")
ecdf_mu_general <- ggplot(draws_with_actual, aes(x = Rank)) + geom_line(aes(y = Value), color = "black") + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, fill = "green", linetype = 0) + facet_wrap(vars(variable), nrow = 1, labeller = label_parsed)

# Figure S1 in supporting information
#jpeg("sbc_mu_ecdf90_5.jpg", height=10, width=30, res=600, pointsize=1, units="cm")
ecdf_mu_general
#dev.off()

##Bias plots

set.seed(5)

nsamples <- 1000
draws <- data.frame(matrix(0, nrow = nsamples, ncol = 3))
names(draws) <- c("mu[1]", "mu[2]", "mu[3]")

alpha_P <- matrix(1, nrow = 3, ncol = 3)
for (i in 1:3) {
  alpha_P[i,i] <- 15
}

for (i in 1:nsamples){
  mu <- rgamma(3, shape = 1, rate = 0.001)
  Pt <- sapply(1:3, function(i) {rdirichlet(1, alpha_P[i,])})
  biases <- mu - Pt%*%mu
  draws[i,] <- biases
}

biasdat <- melt(draws)
names(biasdat)[2] <- "Bias"

# Figure S3 in supporting information
bias_plt <- ggplot(biasdat) + geom_histogram(aes(x = Bias), binwidth = 100, color = "#005b96", fill = '#6497b1', linewidth = 0.25) + facet_wrap(vars(variable), nrow = 1, labeller = label_parsed) + ylab("Count")
#jpeg("biasplot3.jpg", height=10, width=30, res=600, units="cm")
bias_plt
#dev.off()

#Trimmed bias plot
bias_plt_tr <- ggplot(biasdat[which(abs(biasdat$Bias) < 100),]) + geom_histogram(aes(x = Bias), binwidth = 10, color = "#005b96", fill = '#6497b1', linewidth = 0.25) + facet_wrap(vars(variable), nrow = 1, labeller = label_parsed) + ylab("Count")
#jpeg("biasplottr2.jpg", height=10, width=30, res=600, units="cm")
bias_plt_tr
#dev.off()
