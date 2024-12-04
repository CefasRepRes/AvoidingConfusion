library(rstan)
library(bayesplot)
library(posterior)
library(abind)
library(ggplot2)
library(patchwork)

rstan_options(auto_write = T)
ncores <- parallel::detectCores() - 1   #Set this to number of cores available on machine
options(mc.cores = ncores)

#Reading plankton data

plank <- read.csv("plank_data.csv", header=TRUE, na.strings="")
plank_dat <- plank[c(1,2,5:19), c(3,4,5)] #Removing points with very low counts and discarding irrelevant columns
parnames <- c("Copepods", "Detritus", "Non-copepods")

#Setting parameters for MCMC

L <- dim(plank_dat)[1] #Number of locations
M <- dim(plank_dat)[2] #Number of categories

C <- matrix(c(909, 7, 39, 71, 3977, 74, 42, 16, 547), M, M, byrow=T) #Reported confusion matrix

w_dat_plank <- list(L = L, M = M, w = plank_dat, C = C)

#Running MCMC

set.seed(5)

st_model <- stan_model(file = "plankton_model.stan")

#Ignore warnings about R-hat and effective sample size which are due to two badly mixing chains (this issue is addressed below)
mcmc.fit <- sampling(st_model, data = w_dat_plank, iter = 10000, chains = 20, control=list(max_treedepth=12,adapt_delta=0.95), warmup = 5000, seed = 2)

#Extracting inconsistent chains

ex.fit <- extract(mcmc.fit, permuted = FALSE)
## Chains 1 and 15 are bad chains
bad_chains <- as.vector(which(max(ex.fit[5000,,152]) - ex.fit[5000,,152]>1000))

####Bayesplot diagnostics (with only good chains kept)

posterior <- as.matrix(mcmc.fit)
posterior.reduced <- as.array(mcmc.fit)[,-bad_chains,] #Removing bad chains
npars <- dim(posterior)[2]

#Effective sample size and Rhat

min(sapply(1:npars, function(i) {ess_basic(posterior.reduced[,,i])})) #25624 for stan seed = 2 (R seed = 5)
max(sapply(1:npars, function(i) {Rhat(posterior.reduced[,,i])})) #1.00078 for stan seed = 2 (R seed = 5)

###Plots

dimnames(posterior.reduced)$parameters[2] <- "k"
dimnames(posterior.reduced)$parameters[71] <- "Copepods"
dimnames(posterior.reduced)$parameters[72] <- "Detritus"
dimnames(posterior.reduced)$parameters[73] <- "Non-copepods"

##Trace plots
#In supporting information

#jpeg("plank_NBD_muk_trace1.jpg", height=8, width=15, res=600, pointsize=3, units="cm")
mcmc_trace(posterior.reduced, pars = c("mu", "k"), transformations = list(mu = "log"), facet_args = list(labeller = ggplot2::label_parsed, ncol = 1, strip.position= "left"))
#dev.off()

#jpeg("plank_NBD_P_trace1.jpg", height=8, width=15, res=600, pointsize=6, units="cm")
mcmc_trace(posterior.reduced, pars = c("P[1,1]", "P[2,2]", "P[3,3]"), facet_args = list(labeller = ggplot2::label_parsed, ncol = 1, strip.position= "left"),)
#dev.off()

##Histograms

#Histogram of posterior distributions of mu and k
#Not included in paper

#jpeg("plank_NBD_muk_hist6.jpg", height=8, width=15, res=600, pointsize=2, units="cm")
mcmc_hist(posterior.reduced, pars = c("mu", "k"), facet_args = list(labeller = ggplot2::label_parsed), freq = T)
#dev.off()

#Histograms of posterior category means
#Not included in paper

#jpeg("plank_NBD_y_hist4.jpg", height=8, width=15, res=600, pointsize=4, units="cm")
mcmc_hist(posterior.reduced, pars = c("Copepods", "Detritus", "Non-copepods"), freq = T) #+ labs(title = "Category means")
#dev.off()

#Histograms of posterior correlations between categories
# Not im paper

dimnames(posterior.reduced)$parameters[93] <- "Copepods and detritus"
dimnames(posterior.reduced)$parameters[94] <- "Copepods and non-copepods"
dimnames(posterior.reduced)$parameters[97] <- "Non-copepods and detritus"

#jpeg("plank_NBD_corrhist5a.jpg", height=10, width=30, res=600, pointsize=4, units="cm")
mcmc_hist(posterior.reduced, pars = c("Copepods and detritus", "Copepods and non-copepods", "Non-copepods and detritus"), freq = T)
#dev.off()

#Histogram of posterior predictive totals (z-values)
# Not included in the paper

zs <- apply(posterior.reduced, c(1,2), function(x) {rnbinom(1, size = x[2], mu = exp(x[1]))})
zs.df <- data.frame(Total = as.numeric(zs))
zplt <- ggplot(zs.df) + geom_histogram(aes(x = Total, y = ..density..), color = "#005b96", fill = '#6497b1', linewidth = 0.25) + scale_x_log10(labels = function(x) format(x, scientific = F)) + ylab("Density")

#jpeg("ztotalwide.jpg", height = 15, width = 25, res = 600, units = "cm")
zplt
#dev.off()

##Proportion and ratio plots with kernel density estimates

nsamps <- dim(posterior.reduced)[1]
nchains <- dim(posterior.reduced)[2]
drawlength <- nsamps * nchains
ratiodat <- data.frame(matrix(0, nrow = drawlength * (M^2), ncol = 4))
datmeans <- colMeans(plank_dat)
indices <- c(68, 69, 70)
ratios <- array(0, dim = c(M, M, drawlength))

for (i in 1:M){
  for (j in 1:M){
    if (i > j) {
      ratios[i,j,] <- apply(posterior.reduced, c(1,2), function(x) {x[indices[j]]/(x[indices[i]] + x[indices[j]])})
    }
    if (i < j){
      ratios[i,j,] <- apply(posterior.reduced, c(1,2), function(x) {x[indices[j]]/x[indices[i]]})
    }
    else {
      ratios[i,i,] <- apply(posterior.reduced, c(1,2), function(x) {x[indices[i]]/sum(x[indices])})
    }
  }
}

for (i in 1:M){
  for (j in 1:M){
    ratiodat[(((i-1)*M + j - 1)*drawlength + 1):(((i-1)*M + j)*drawlength),1] <- ratios[i,j,]
    ratiodat[(((i-1)*M + j - 1)*drawlength + 1):(((i-1)*M + j)*drawlength),2] <- rep(parnames[i], drawlength)
    ratiodat[(((i-1)*M + j - 1)*drawlength + 1):(((i-1)*M + j)*drawlength),3] <- rep(parnames[j], drawlength)
    if (i == j) {
      ratiodat[(((i-1)*M + j - 1)*drawlength + 1):(((i-1)*M + j)*drawlength),4] <- rep(datmeans[i]/sum(datmeans), drawlength)
    }
    if (i > j) {
      ratiodat[(((i-1)*M + j - 1)*drawlength + 1):(((i-1)*M + j)*drawlength),4] <- rep(datmeans[j]/(datmeans[i] + datmeans[j]), drawlength)
    }
    if (i < j) {
      ratiodat[(((i-1)*M + j - 1)*drawlength + 1):(((i-1)*M + j)*drawlength),4] <- rep(datmeans[j]/datmeans[i], drawlength)
    }
  }
}

ratiopltlist <- list(list(),list())
counter <- rep(0,2)
for (i in 1:M) {
  for (j in 1:M) {
    if (i == j) {
      counter[1] <- counter[1] + 1
      tmp <- ggplot(ratiodat) + geom_density(data = ratiodat[intersect(which(ratiodat$X2 == parnames[i]), which(ratiodat$X3 == parnames[j])),-2], aes(x = X1, y = ..density..), color = "#005b96", fill = '#6497b1', linewidth = 0.25) + geom_vline(data = ratiodat[intersect(which(ratiodat$X2 == parnames[i]), which(ratiodat$X3 == parnames[j])),-2], aes(xintercept = X4), color = "red", linewidth = 0.75) + xlab(parnames[i]) + ylab("Density")
      ratiopltlist[[1]][[counter[1]]] <- tmp
    }
    if (i > j) {
      counter[2] <- counter[2] + 1
      tmp <- ggplot(ratiodat) + geom_density(data = ratiodat[intersect(which(ratiodat$X2 == parnames[i]), which(ratiodat$X3 == parnames[j])),], aes(x = X1, y = ..density..), color = "#005b96", fill = '#6497b1', linewidth = 0.25) + geom_vline(data = ratiodat[intersect(which(ratiodat$X2 == parnames[i]), which(ratiodat$X3 == parnames[j])),], aes(xintercept = X4), color = "red", linewidth = 0.75) + xlab(parnames[j]) + ylab(parnames[i]) + scale_x_continuous(limits = c(0,1))
      ratiopltlist[[2]][[counter[2]]] <- tmp
    }
  }
}

propplt <- ratiopltlist[[1]][[1]]
for (p in 2:length(ratiopltlist[[1]])) {
  propplt <- propplt + ratiopltlist[[1]][[p]]
}

#3-panel proportion plot
## Figure 7 in paper
propplt <- propplt + plot_layout(axis_titles = "collect")
#jpeg("meanpltwide.jpg", height = 10, width = 30, res = 600, units = "cm")
propplt
#dev.off()

ratioplt <- ratiopltlist[[2]][[1]]
for (p in 2:length(ratiopltlist[[2]])) {
  ratioplt <- ratioplt + ratiopltlist[[2]][[p]]
}

des <- "
  1#
  32
"

#Triangular ratio plot
ratioplt <- ratioplt + plot_layout(design = des, axis_titles = "collect")
#jpeg("ratiopltwide.jpg", height = 15, width = 25, res = 600, units = "cm")
ratioplt
#dev.off()

## Figure 8 in paper
#Bottom left corner
ratioplt2 <- ratiopltlist[[2]][[2]] + xlab("Proportion") + ylab("Density")
#jpeg("ratioplt3wide.jpg", height = 15, width = 25, res = 600, units = "cm")
ratioplt2
#dev.off()


####Old 9-panel ratio plot

ratiopltlist <- list()
counter = 0
for (i in 1:M) {
  for (j in 1:M) {
    counter = counter + 1
    if (i == j) {
      tmp <- ggplot(ratiodat) + geom_density(data = ratiodat[intersect(which(ratiodat$X2 == parnames[i]), which(ratiodat$X3 == parnames[j])),], aes(x = X1, y = ..density..), color = "#005b96", fill = '#6497b1', linewidth = 0.25) + geom_vline(data = ratiodat[intersect(which(ratiodat$X2 == parnames[i]), which(ratiodat$X3 == parnames[j])),], aes(xintercept = X4), color = "red", linewidth = 0.75) + xlab(parnames[j]) + ylab(parnames[i])
      ratiopltlist[[counter]] <- tmp
    }
    if (i < j) {
      tmp <- ggplot(ratiodat) + geom_density(data = ratiodat[intersect(which(ratiodat$X2 == parnames[i]), which(ratiodat$X3 == parnames[j])),], aes(x = X1, y = ..density..), color = "#005b96", fill = '#6497b1', linewidth = 0.25) + geom_vline(xintercept = 1, color = "#03396c", linewidth = 0.5, linetype = "dashed") + geom_vline(data = ratiodat[intersect(which(ratiodat$X2 == parnames[i]), which(ratiodat$X3 == parnames[j])),], aes(xintercept = X4), color = "red", linewidth = 0.75) + scale_x_log10(limits = c(0.01,100), labels = function(x) format(x, scientific = F)) + xlab(parnames[j]) + ylab(parnames[i])
      ratiopltlist[[counter]] <- tmp
    }
    else {
      tmp <- ggplot(ratiodat) + geom_density(data = ratiodat[intersect(which(ratiodat$X2 == parnames[i]), which(ratiodat$X3 == parnames[j])),], aes(x = X1, y = ..density..), color = "#005b96", fill = '#6497b1', linewidth = 0.25) + geom_vline(data = ratiodat[intersect(which(ratiodat$X2 == parnames[i]), which(ratiodat$X3 == parnames[j])),], aes(xintercept = X4), color = "red", linewidth = 0.75) + xlab(parnames[j]) + ylab(parnames[i]) + scale_x_continuous(limits = c(0,1))
      ratiopltlist[[counter]] <- tmp
    }
  }
}

ratioplt <- ratiopltlist[[1]]
for (p in 2:length(ratiopltlist)) {
  ratioplt <- ratioplt + ratiopltlist[[p]]
}

ratioplt <- ratioplt + plot_layout(nrow = 3, axis_titles = "collect")

#jpeg("ratioplot_withlines_dens.jpg", height=15, width=25, res=600, units="cm")
ratioplt
#dev.off()

###90% HPD credible intervals and mean values

q = round(quantile(exp(posterior.reduced[,,1]), probs = c(0.05, 0.95)), 1) #credible interval for mu
qk = round(quantile(posterior.reduced[,,2], probs=c(0.05, 0.95)), 1) #credible interval for k
q1 = round(quantile(posterior.reduced[,,71], probs=c(0.05, 0.95)), 1) #credible interval for copepods
q2 = round(quantile(posterior.reduced[,,72], probs=c(0.05, 0.95)), 1) #credible interval for detritus
q3 = round(quantile(posterior.reduced[,,73], probs=c(0.05, 0.95)), 1) #credible interval for non-copepods
qcorr1 = round(quantile(posterior.reduced[,,93], probs=c(0.05, 0.95)), 1) #credible interval for correlation between copepods and detritus
qcorr2 = round(quantile(posterior.reduced[,,94], probs=c(0.05, 0.95)), 1) #credible interval for correlation between copepods and non-copepods
qcorr3 = round(quantile(posterior.reduced[,,97], probs=c(0.05, 0.95)), 1) #credible interval for correlation between non-copepods and detritus

m1 <- round(mean(posterior.reduced[,,71]), 1) #Posterior mean of copepods expected value
m2 <- round(mean(posterior.reduced[,,72]), 1) #Posterior mean of detritus expected value
m3 <- round(mean(posterior.reduced[,,73]), 1) #Posterior mean of non-copepods expected value

q 
qk 
q1 
q2 
q3
qcorr1 
qcorr2
qcorr3 

m1 #453.7
m2 #7812.4
m3 #1273.1








