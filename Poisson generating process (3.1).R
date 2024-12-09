source("ML functions.R")
# You may need to install the package "vectrs" to get bayesplot to work
library(rstan)
library(bayesplot)
library(posterior)
library(ggplot2)
library(reshape2)
library(MCMCprecision)

rstan_options(auto_write = T)
ncores <- 4  #Change to less than 4 if there are less than 4 cores available on machine
options(mc.cores = ncores)

##Main run

#Specify the underlying probability matrix and vector of row means

P_real <- matrix(c(0.9, 0.05, 0.05, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8), nrow = 3, ncol = 3,
                 byrow = T)
mu_real <- c(10, 100, 50)

set.seed(5)

#Generate C (original confusion matrix estimate) from s (true category counts
#in training set)

s <- c(50, 200, 100)
C <- Cgen(s, P_real)
C

#Generate prior matrix parameters

alpha_P <- matrix(1, nrow = 3, ncol = 3)

#Simulate N = 15 observations according to above parameters

N <- 15
distribution = rep("Poisson", length(mu_real))
Wobs <- observedW.gen.multi(n = N, mu = mu_real, pmat = P_real, distribution = distribution)
Wobs

#Calculate column means

Wmeans <- colMeans(Wobs)
Wmeans

#Running stan model

w_dat <- list(N = s, L = N, M = 3, P = P_real, w = Wobs, C = C, alpha_P = alpha_P)

st_model <- stan_model(file = "poisson_model.stan")

mcmc.fit <- sampling(st_model, data = w_dat, iter=4000, seed = 5)

##Model checking

posterior <- as.matrix(mcmc.fit)
q1 = round(quantile(posterior[,1], probs=c(0.05, 0.95)), 1) #90% credible interval for category 1 mean
q2 = round(quantile(posterior[,2], probs=c(0.05, 0.95)), 1) #90% credible interval for category 2 mean
q3 = round(quantile(posterior[,3], probs=c(0.05, 0.95)), 1) #90% credible interval for category 3 mean

q1; q2; q3

##Bayesplot diagnostics

#Preparing data for plotting

arr.fit <- as.array(mcmc.fit)
mat.fit <- as.matrix(mcmc.fit)
dfmu.fit <- data.frame(mat.fit[,1:3])
mu_names <- c("mu[1]", "mu[2]", "mu[3]")
names(dfmu.fit) <- mu_names
dfmu.plt <- melt(dfmu.fit) # ignore R message about "No id variables ..."
dfmu.plt$true <- unlist(lapply(1:3, function(x) {rep(mu_real[x], nrow(dfmu.plt)/3)}))
dfmu.plt$wmeans <- unlist(lapply(1:3, function(x) {rep(Wmeans[x], nrow(dfmu.plt)/3)}))

#Effective sample size and Rhat

npars <- dim(arr.fit)[3]
max(sapply(1:npars, function(i) {rstan::Rhat(arr.fit[,,i])})) #3277 (R seed = 5, stan seed =  5)
min(sapply(1:npars, function(i) {ess_basic(arr.fit[,,i])})) #1.00186 (R seed = 5, stan seed = 5)

#Trace plots of (post burn-in) MCMC draws for the category mean parameters

par_trace <- mcmc_trace(arr.fit, regex_pars = "mu\\[[1,2,3]\\]", facet_args = list(labeller = ggplot2::label_parsed, ncol = 1, strip.position = "left"))
par_trace

#Histogram of posterior draws of each of the category mean parameters with true values plotted on

compare_hist <- mcmc_recover_hist(arr.fit[,,c(1, 2, 3)], mu_real, facet_args = list(labeller = ggplot2::label_parsed)) + xlab("Parameter value") + ylab("Number of samples")
compare_hist

# Violin plots for each of the category mean parameters with true values plotted on
# (white dot) and also the mean of the w data (red dot)
# Included in paper as Figure 2

par_viol_withpoints <- ggplot(dfmu.plt) + geom_violin(aes(variable, value), fill = "#6497b1", scale = "width") + geom_point(aes(variable, true), color = "white", size = 2) + geom_point(aes(variable, wmeans), color = "#6497b1", fill = "red", shape = 23, size = 2.5) + scale_x_discrete(labels = expression(mu[1], mu[2], mu[3])) + xlab("Parameter") + ylab("Value")
#jpeg("viol_mu_Pois_withpoints4.jpg", height=15, width=15, res=600, pointsize=1, units="cm")
par_viol_withpoints
#dev.off()

