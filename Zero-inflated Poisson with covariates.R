source("ML functions.R")
# You may need to install the package "vectrs" to get bayesplot to work
library(rstan)
library(bayesplot)
library(posterior)
library(ggplot2)
library(reshape2)
library(MCMCprecision)
library(patchwork)
library(VGAM)

rstan_options(auto_write = T)
ncores <- 4  #Change to less than 4 if there are less than 4 cores available on machine
options(mc.cores = ncores)

##Main run

#Specify the underlying probability matrix and vector of row means

P_real <- matrix(c(0.7, 0.3, 0.15, 0.85), nrow = 2, ncol = 2,
                 byrow = T)
mu_real <- c(100, 200)
q_real <- c(0.1, 0.4)

set.seed(5)

#Generate C (original confusion matrix estimate) from s (true category counts
#in training set)

s <- c(5000, 1500)
C <- Cgen(s, P_real)
C

#Set prior matrix parameters

alpha_P <- matrix(1, nrow = 2, ncol = 2)
alpha_q <- matrix(1, nrow = 2, ncol = 2)

#Specify true covariate relationship

intercepts_real <- c(3,4)
betas_real <- c(2, 1.6)

#Simulate data

## Perhaps wrong number of mus being passed to observed.gen.
## Need to sort that out. See the file in ML functions gull.R in the root
## directory. This may well be the correct one. Transfer to ML functions.R
## in the code directory

N <- 15
covs <- rnorm(N, 0, 1)
distribution = rep("ZIPoisson", N)

Wobs2 <- observedW.gen.multi.covar(n = N, pmat = P_real, distribution = distribution, q = q_real, covs = covs, intercepts = intercepts_real, betas = betas_real)
Wobs2

#Running stan model

w_dat <- list(L = N, w = Wobs2, C = C, alpha_P = alpha_P, alpha_q = alpha_q, X = covs)

st_model <- stan_model(file = "zeroinf_poisson_model.stan")

mcmc.fit <- sampling(st_model, data = w_dat, iter=4000,chains=4, seed = 5)

##Model checking

#Preparing data for plotting

arr.fit <- as.array(mcmc.fit)
mat.fit <- as.matrix(mcmc.fit)
df.plt.list <- list()
coef_names <- c("p[1]", "p[2]", "q[1]", "q[2]", "alpha[1]", "alpha[2]", "beta[1]", "beta[2]")
for (i in 1:4){
  df.fit <- data.frame(mat.fit[,2*(i-1) + (1:2)])
  names(df.fit) <- coef_names[(2*i-1):(2*i)]
  df.plt.list[[i]] <- melt(df.fit) # ignore R message about "No id variables ..."
}

df.plt.list[[1]]$true <- unlist(lapply(1:2, function(x) {rep(P_real[x,x], nrow(df.plt.list[[1]])/2)}))
df.plt.list[[2]]$true <- unlist(lapply(1:2, function(x) {rep(q_real[x], nrow(df.plt.list[[2]])/2)}))
df.plt.list[[3]]$true <- unlist(lapply(1:2, function(x) {rep(intercepts_real[x], nrow(df.plt.list[[3]])/2)}))
df.plt.list[[4]]$true <- unlist(lapply(1:2, function(x) {rep(betas_real[x], nrow(df.plt.list[[4]])/2)}))

zeroinfplot <- ggplot(df.plt.list[[1]]) + geom_violin(aes(variable, value), fill = "#6497b1", scale = "width") + geom_point(aes(variable, true), color = "white", size = 2)+ scale_x_discrete(labels = expression(P[list(1,1)], P[list(2,2)])) + xlab("Parameter") + ylab("Value")
zeroinfplot <- zeroinfplot + ggplot(df.plt.list[[2]]) + geom_violin(aes(variable, value), fill = "#6497b1", scale = "width") + geom_point(aes(variable, true), color = "white", size = 2)+ scale_x_discrete(labels = expression(q[1], q[2])) + xlab("Parameter") + ylab("Value")
zeroinfplot <- zeroinfplot + ggplot(df.plt.list[[3]]) + geom_violin(aes(variable, value), fill = "#6497b1", scale = "width") + geom_point(aes(variable, true), color = "white", size = 2)+ scale_x_discrete(labels = expression(alpha[1], alpha[2])) + xlab("Parameter") + ylab("Value")
zeroinfplot <- zeroinfplot + ggplot(df.plt.list[[4]]) + geom_violin(aes(variable, value), fill = "#6497b1", scale = "width") + geom_point(aes(variable, true), color = "white", size = 2)+ scale_x_discrete(labels = expression(beta[1], beta[2])) + xlab("Parameter") + ylab("Value")

#Violin plots
zeroinfplot <- zeroinfplot + plot_layout(axis_titles = "collect")
#jpeg("ziplot_withP.jpg", height = 15, width = 15, res = 600, units = "cm")
zeroinfplot
#dev.off()

#Effective sample size and Rhat

npars <- dim(arr.fit)[3]
max(sapply(1:npars, function(i) {rstan::Rhat(arr.fit[,,i])})) #1.002166 (R seed = 2, stan seed =  5, N = 100) 1.001077(R seed = 5, stan seed = 5, N = 15)
min(sapply(1:npars, function(i) {ess_basic(arr.fit[,,i])})) #3542.353 (R seed = 2, stan seed = 5) 3667.139 (R seed = 5, stan seed = 5, N = 15)

