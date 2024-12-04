#Functions used in code
source("ML functions.R")

library(MCMCprecision)
library(parallel)
library(bayesplot)
library(rstan)
library(posterior)
library(ggplot2)
library(patchwork)
library(reshape2)
library(MASS)

#Specify required parameters/arguments

distributions = c(rep("Negbin", 2), rep("Poisson", 3))
#Means for each of the five generating distributions
mu = c(10, 100, 50, 200, 30)
nmu = length(mu)
#Size parameter for the two Negative Binomial distributions
k = c(2, 5)
#Initial probability matrix
Pinit = matrix(0, ncol=5, nrow=5)
Pinit[1,] = c(0.8, 0.05, 0.05, 0.05, 0.05) 
Pinit[2,] = c(0.02, 0.9, 0.03, 0.01, 0.04)
Pinit[3,] = c(0.05, 0.03, 0.85, 0.05, 0.02)
Pinit[4,] = c(0.01, 0.01, 0.01, 0.95, 0.02)
Pinit[5,] = c(0.03, 0.01, 0.05, 0.01, 0.9)

#Simulate data W according to above choices
#The matrix Wobs has n rows and nmu columns, where nmu is the
#length of the vector mu.

set.seed(5)
Wobs = observedW.gen.multi(n=15, mu, Pinit, distributions, k)
Wmeans <- colMeans(Wobs)

# Create an initial confusion matrix

CM = matrix(0, ncol=nmu, nrow=nmu)
for (j in 1:nmu) {
  CM[j,] = rmultinom(1, size=100, prob=Pinit[j,])
}
CM

# Choose geometric jump parameters for MCMC sampling from the X matrix

ggprobs <- geomprobs.gen(Wobs, delta = 0.9, Unif = F, sig = 4)
ggprobs

# Extra arguments for main MCMC function
# These commented-out ones are the ones used in the paper.
# Lesser values used for testing

#thin = 100
#Nits = 2000000
#Nthin = Nits / thin   #Make sure this is an integer
#burnin = 10000     #Must be less than Nthin

thin = 10
Nits = 20000
Nthin = Nits / thin   #Make sure this is an integer
burnin = 1000     #Must be less than Nthin

# Choose standard deviations of proposal distributions for each of mu and k.

cand_sd_mu = c(8.0, 0.5)
cand_sd_k = c(8.0, 1.5)

# Is proposal on log scale or not?

log_prop_mu = c(F, T)
log_prop_k = c(F, T)

#Multiple chain setup

M <- dim(Wobs)[2]
nreps <- dim(Wobs)[1]
nchains <- 6
chain_names <- sapply(1:nchains, function(x) {paste("chain", x, sep = ":")})

trials <- array(0, dim = c(Nthin-burnin, nchains, 5+2+25+(15*25)))
names_mu <- c("mu[1]", "mu[2]", "mu[3]", "mu[4]", "mu[5]")
names_k <- c("k[1]", "k[2]")
matlabs <- expand.grid(1:5, 1:5)
names_P <- sapply(1:25, function(x) {paste0("P[", matlabs[x,1],",", matlabs[x,2],"]")})
names_X <- unlist(lapply(1:15, function(y) {sapply(1:25, function(x){paste0("X", y, "[", matlabs[x,1],",", matlabs[x,2],"]")})}))
par_names <- c(names_mu, names_k, names_P, names_X)
dimnames(trials) <- list(Iteration = NULL, Chain = chain_names, Parameter = par_names)

#Fitting multiple chains in parallel

n_cores <- detectCores()
cl <- makeCluster(min(n_cores-1, nchains))

clusterEvalQ(cl, library(MCMCprecision))
clusterExport(cl, c("dlike_log_NB_Po", "dlike_log_NB_Po_rep", "Xinitial", "knewprior", "updateX_NB_Po", "updatemuk_NB_Po"))
trial_list <- clusterCall(cl, mcmcML_NB_Po2, Wdata=Wobs, C_matrix=CM, distributions = distributions, Nits=Nits, thin=thin,
                          geom.probs=ggprobs, cand_sd_mu=cand_sd_mu,
                          cand_sd_k=cand_sd_k, log_prop_mu=log_prop_mu,
                          log_prop_k=log_prop_k, fit_init=F, C_sam=T)

stopCluster(cl)

for (chain in 1:nchains) {
  
  P_reduced  <- t(apply(trial_list[[chain]]$P[,,(burnin+1):Nthin], 3, as.vector))
  X_reduced <- t(apply(trial_list[[chain]]$X[,,(burnin+1):Nthin,], 3, as.vector))
  trials[,chain,] <- cbind(trial_list[[chain]]$mu[(burnin+1):Nthin,], trial_list[[chain]]$k[(burnin+1):Nthin,], P_reduced, X_reduced)
}

###Bayesplot and Rstan diagnostics

#Acceptance probabilities for each X matrix and each joint (mu,k) proposed jump.

lapply(trial_list, function(chain) {chain$Xprob}) #Between 15% and 25%
lapply(trial_list, function(chain) {chain$mukprob}) #Between 3% and 5%

#Effective sample size and Rhat

npars <- dim(trials)[3]
max(sapply(1:npars, function(i) {rstan::Rhat(trials[,,i])})) #1.001
min(sapply(1:npars, function(i) {ess_basic(trials[,,i])})) #5360

##Plots

#Preparing MCMC output for plotting

posterior <- do.call(rbind, lapply(1:dim(trials)[2], function(j) {trials[,j,]}))
dfmu.fit <- data.frame(posterior[,1:5])
dfk.fit <- data.frame(posterior[,6:7])
mu_names <- c("mu[1]", "mu[2]", "mu[3]", "mu[4]", "mu[5]")
k_names <- c("k[1]", "k[2]")
dfmu.plt <- melt(dfmu.fit)
dfk.plt <- melt(dfk.fit)
dfmu.plt$true <- unlist(lapply(1:5, function(x) {rep(mu[x], nrow(dfmu.plt)/5)}))
dfk.plt$true <- unlist(lapply(1:2, function(x) {rep(k[x], nrow(dfk.plt)/2)}))
dfmu.plt$wmeans <- unlist(lapply(1:5, function(x) {rep(Wmeans[x], nrow(dfmu.plt)/5)}))

#Histograms of posterior draws of mu components with true values plotted on
#Not included in paper

mu_hist_means <- mcmc_recover_hist(trials[,,1:5], mu[1:5], facet_args = list(labeller = ggplot2::label_parsed)) + yaxis_text(on = T)
mu_hist_means

#Histograms of posterior draws of k components with true values plotted on
#Not included in paper

k_hist_means <- mcmc_recover_hist(trials[,,6:7], k, facet_args = list(labeller = ggplot2::label_parsed)) + yaxis_text(on = T)
k_hist_means

#Histograms of posterior draws of P entries with true values plotted on
#Not included in paper

true_P <-  as.vector(Pinit)
P_hist_means <- mcmc_recover_hist(trials[,,8:32], true_P)
P_hist_means

#Trace plots of posterior draws of mu components
#S3 in supporting information

mu_trace <- mcmc_trace(trials, regex_pars = "mu\\[[1,2,3,4,5]\\]", facet_args = list(labeller = ggplot2::label_parsed, ncol = 1, strip.position = "left"))
#jpeg("toy_PNB_trace_mu_FF2.jpg", height=12, width=15, res=600, pointsize=1, units="cm")
mu_trace
#dev.off()

#Trace plots of posterior draws of k components
# S4 in supporting in formation

k_trace <- mcmc_trace(trials, regex_pars = "k\\[[1,2]\\]", facet_args = list(labeller = ggplot2::label_parsed, ncol = 1, strip.position = "left"))
k_trace

#Violin plots of posterior draws of mu components with true values plotted (white
# dots) and mean of observed data (red dots).

par_viol_mu_withpoints <- ggplot(dfmu.plt) + geom_violin(aes(variable, value), fill = "#6497b1", scale = "width") + geom_point(aes(variable, true), color = "white", size = 2) + geom_point(aes(variable, wmeans), color = "#6497b1", fill = "red", size = 2.5, shape = 23) + scale_x_discrete(labels = expression(mu[1], mu[2], mu[3], mu[4], mu[5])) + xlab("Parameter") + ylab("Value")
par_viol_mu_withpoints

#Violin plots of posterior draws of k components with true values plotted (white
# dots)

par_viol_k_withpoints <- ggplot(dfk.plt) + geom_violin(aes(variable, value), fill = "#6497b1", scale = "width") + geom_point(aes(variable, true), color = "white", size = 2) + scale_x_discrete(labels = expression(k[1], k[2])) + xlab("Parameter") + ylab("Value")
par_viol_k_withpoints

#Combined violin plots of mu and k.
# Figure 3 in paper

combplt_withpoints <- par_viol_mu_withpoints + par_viol_k_withpoints + plot_layout(widths = c(2, 1), axes = "collect")
#jpeg("viol_muk_withpoints3.jpg", height=16, width=24, res=600, pointsize=1, units="cm")
combplt_withpoints
#dev.off()

##Credible intervals

## 90% credible interval for each component of mu
quantile(trials[,,1], c(0.05,0.95))
quantile(trials[,,2], c(0.05,0.95))
quantile(trials[,,3], c(0.05,0.95))
quantile(trials[,,4], c(0.05,0.95))
quantile(trials[,,5], c(0.05,0.95))

## 90% credible interval for each component of k
quantile(trials[,,6], c(0.05,0.95))  
quantile(trials[,,7], c(0.05,0.95))

