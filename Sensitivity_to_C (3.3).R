#Software to investigate the effect of the number of observations (in s)
#used to generate the confusion matrix (Section 3.3 in paper)
source("ML functions.R")
library(rstan)
rstan_options(auto_write = T)
options(mc.cores = 4)
library(ggplot2)
library(patchwork)
library(MCMCprecision)
library(reshape2)

####Required functions

inv_logit <- function(x){
  1 / (1 + exp(-x))
}

trans_funct <- function(theta){
  c(exp(theta[1:2]),inv_logit(theta[3:4]))
}

## IS THIS THE SAME AS ON OF THE LIKELIHOOD FUNCTIONS IN "ML FUNCTIONS.R"
## DO WE NEED IT?

log_lik <- function(theta,w,C,n_1,n_2){
  mu <- exp(theta[1:2])
  p <- inv_logit(theta[3:4])
  sum(dpois(w,p*mu + (1-p[2:1])*mu[2:1],log=T)) + dbinom(C[1,1],n_1,p[1],log=T) + dbinom(C[2,2],n_2,p[2],log=T)
}

priors <- function(theta){
  mu <- exp(theta[1:2])
  p <- inv_logit(theta[3:4])
  sum(dgamma(mu,1,1e-2,log=T)) + sum(dbeta(p,1.0,1.0,log=T)) 
}

###For rejection sampling

get_mu <- function(lambda,ps){
  mu <- rep(NA,2)
  mu[1] <- (ps[2] * lambda[1] - (1-ps[2]) * lambda[2])/ (ps[1] + ps[2] - 1)
  mu[2] <- (1 / ps[2]) * (lambda[2] - (1 - ps[1]) * mu[1])
  return(mu)
}

rel_prior <- function(mu){
  prod(dgamma(mu,1,1e-4)/dgamma(0,1,1e-4))
}

rej_sample <- function(w,C,ns,N_iter){
  lambda <- matrix(rgamma(2*N_iter,1 + rowSums(w),ncol(w)),nrow=2)
  ps_ <- matrix(rbeta(2*N_iter,1+C,1+ns-C),nrow=2)
  mus_ <- sapply(1:N_iter,function(x,lambda,ps){
    get_mu(lambda[,x],ps[,x])
  },ps=ps_,lambda=lambda)
  
  probs <- apply(mus_,2,rel_prior)
  accept <- which(runif(N_iter) < probs)
  
  ps <- t(ps_[,accept])
  mus <- t(mus_[,accept])
  
  return(cbind(mus,ps))
}

###Main run setup

mu_real <- c(20, 80)
P_real <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, ncol = 2,
                 byrow = T)
N <- 5

set.seed(10)
distribution = rep("Poisson", length(mu_real))
Wobs2 <- t(observedW.gen.multi(n = N, mu = mu_real, pmat = P_real, distribution = distribution))
Wobs2

mu_draws_list <- list()
mu_draws_hex_list <- list()
P_draws_list <- list()
P_draws_hex_list <- list()
outputs <- list()
coef_names_mu <- c("mu1", "mu2")
coef_names_P <- c("p1","p2")

#Storing possible C matrices

Cs <- list()
rowtotals <- c(0, 5, 10, 100, 1000)
ntest <- length(rowtotals)
for (i in 1:ntest){
  Cs[[i]] <- list()
  for (j in 1:ntest){
    if (i != 1){
      if (j != 1){
        s <- c(rowtotals[i]-rowtotals[i-1], rowtotals[j]-rowtotals[j-1])
        C_tmp <- Cgen(s, P_real)
        Cs[[i]][[j]] <- Cs[[i-1]][[j-1]] + C_tmp 
      }
      else{
        if (rowtotals[1] != 0){
          s <- c(rowtotals[i]-rowtotals[i-1], rowtotals[1])
          C_tmp <- Cgen(s, P_real)
          C_tmp[1,] <- C_tmp[1,] + Cs[[i-1]][[1]][1,]
          Cs[[i]][[1]] <- C_tmp
        }
        else {
          C_tmp <- rbind(t(rmultinom(1, rowtotals[i]-rowtotals[i-1],P_real[1,])), rep(0,nrow(P_real))) 
          C_tmp[1,] <- C_tmp[1,] + Cs[[i-1]][[1]][1,]
          Cs[[i]][[1]] <- C_tmp
        }
      }
    }
    else {
      if (j != 1){
        if (rowtotals[1] != 0){
          s <- c(rowtotals[1], rowtotals[j]-rowtotals[j-1])
          C_tmp <- Cgen(s, P_real)
          C_tmp[2,] <-  C_tmp[2,] + Cs[[1]][[j-1]][2,]
          Cs[[1]][[j]] <- C_tmp
        }
        else {
          C_tmp <- rbind(rep(0,nrow(P_real)), t(rmultinom(1, rowtotals[j]-rowtotals[j-1],P_real[2,]))) 
          C_tmp[2,] <- C_tmp[2,] + Cs[[1]][[j-1]][2,]
          Cs[[1]][[j]] <- C_tmp
        }
      }
      else{
        if (rowtotals[1] == 0){
          Cs[[1]][[1]] <- diag(0, nrow = nrow(P_real), ncol = nrow(P_real))
        }
        else {
          s <- c(rowtotals[1], rowtotals[1])
          Cs[[1]][[1]] <- Cgen(s, P_real)
        }
      }
    }
  }
}



# Main run - ignore R messages about "No id variables ..."
for (i in 1:ntest){
  mu_draws_list[[i]] <- list()
  mu_draws_hex_list[[i]] <- list()
  P_draws_hex_list[[i]] <- list()
  P_draws_list[[i]] <- list()
  outputs[[i]] <- list()
  for (j in 1:ntest){
    res.fit <- list(rej_sample(Wobs2,c(Cs[[i]][[j]][1,1],Cs[[i]][[j]][2,2]),ns=c(rowtotals[i],rowtotals[j]),N_iter=1e5))
    outputs[[i]][[j]] <- res.fit
    dfmu.fit <- data.frame(res.fit[[1]][,1:2])
    dfp.fit <- data.frame(res.fit[[1]][,3:4])
    names(dfmu.fit) <- coef_names_mu
    names(dfp.fit) <- coef_names_P
    dfmu.plt <- melt(dfmu.fit) # ignore R message about "No id variables ..."
    dfp.plt <- melt(dfp.fit)
    dfmu.hex.plt <- melt(dfmu.fit, id.vars = c("mu1","mu2"))
    dfp.hex.plt <- melt(dfp.fit, id.vars = c("p1","p2"))
    dfmu.plt$true <- unlist(lapply(1:2, function(x) {rep(mu_real[x], nrow(dfmu.plt)/2)}))
    dfmu.plt$cat1 <- rep(rowtotals[i], nrow(dfmu.plt))
    dfmu.plt$cat2 <- rep(rowtotals[j], nrow(dfmu.plt))
    dfp.plt$true <- unlist(lapply(1:2, function(x) {rep(P_real[x,x], nrow(dfp.plt)/2)}))
    dfp.plt$cat1 <- rep(rowtotals[i], nrow(dfp.plt))
    dfp.plt$cat2 <- rep(rowtotals[j], nrow(dfp.plt))
    dfmu.hex.plt$cat1 <- rep(rowtotals[i], nrow(dfmu.hex.plt))
    dfmu.hex.plt$cat2 <- rep(rowtotals[j], nrow(dfmu.hex.plt))
    dfp.hex.plt$cat1 <- rep(rowtotals[i], nrow(dfp.hex.plt))
    dfp.hex.plt$cat2 <- rep(rowtotals[j], nrow(dfp.hex.plt))
    mu_draws_list[[i]][[j]] <- dfmu.plt
    P_draws_list[[i]][[j]] <- dfp.plt
    mu_draws_hex_list[[i]][[j]] <- dfmu.hex.plt
    P_draws_hex_list[[i]][[j]] <- dfp.hex.plt
  }
}

###########Plots

#Contour plots (5x5)

dfmu.hex.plt <- do.call(rbind, lapply(mu_draws_hex_list, function(x) {do.call(rbind, x)}))
mu_hexplt <- ggplot(dfmu.hex.plt) + geom_density_2d_filled(aes(mu1,mu2), bins = 11, alpha = 0.8, contour_var = "ndensity", adjust = 2, show.legend = FALSE) + scale_fill_brewer(type = "seq", palette = "RdBu", direction = -1) + geom_point(x=20,y=80,shape = 21, fill = "yellow", color = "black", size = 1.5) + xlab(expression(mu[1])) + ylab(expression(mu[2])) + scale_x_continuous(sec.axis = sec_axis(~., name = expression(s[2]), breaks = NULL, labels = NULL)) + scale_y_continuous(sec.axis = sec_axis(~., name = expression(s[1]), breaks = NULL, labels = NULL)) + facet_grid(rows = vars(cat1), cols = vars(cat2))
#jpeg("mu_contour_withpoints_yellowpoints_1e5.jpg", height=25, width=25, res=600, pointsize=1, units="cm")
mu_hexplt
#dev.off()

dfp.hex.plt <- do.call(rbind, lapply(P_draws_hex_list, function(x) {do.call(rbind, x)}))
P_hexplt <- ggplot(dfp.hex.plt) + geom_density_2d_filled(aes(x=p1,y=p2), bins = 11, alpha = 0.8, contour_var = "ndensity", adjust = 2, show.legend = FALSE) + scale_fill_brewer(type = "seq", palette = "RdBu", direction = -1) + geom_point(x=0.9,y=0.8,shape = 21, fill = "yellow", color="black", size = 1.5) + xlab(expression(P[list(1,1)])) + ylab(expression(P[list(2,2)])) + scale_x_continuous(sec.axis = sec_axis(~., name = expression(s[2]), breaks = NULL, labels = NULL)) + scale_y_continuous(sec.axis = sec_axis(~., name = expression(s[1]), breaks = NULL, labels = NULL)) + facet_grid(rows = vars(cat1), cols = vars(cat2))
#jpeg("P_contour_withpoints_yellowpoints_1e5.jpg", height=25, width=25, res=600, pointsize=1, units="cm")
P_hexplt
#dev.off()

#Reduced contour plots (3x3)

# Figure 4 in paper
mu_draws_hex_list_reduced <- lapply(mu_draws_hex_list, function(x) {x[c(1,3,5)]})
mu_draws_hex_list_reduced <- mu_draws_hex_list_reduced[c(1,3,5)]
dfmu.hex.plt2 <- do.call(rbind, lapply(mu_draws_hex_list_reduced, function(x) {do.call(rbind, x)}))
mu_hexplt2 <- ggplot(dfmu.hex.plt2) + geom_density_2d_filled(aes(mu1,mu2), bins = 11, alpha = 0.8, contour_var = "ndensity", adjust = 2, show.legend = FALSE) + scale_fill_brewer(type = "seq", palette = "RdBu", direction = -1) + geom_point(x=20,y=80, shape = 21, fill = "yellow", color = "black", size = 1.5) + xlab(expression(mu[1])) + ylab(expression(mu[2])) + scale_x_continuous(sec.axis = sec_axis(~., name = expression(s[2]), breaks = NULL, labels = NULL)) + scale_y_continuous(sec.axis = sec_axis(~., name = expression(s[1]), breaks = NULL, labels = NULL)) + facet_grid(rows = vars(cat1), cols = vars(cat2))
#jpeg("mu_contour_withpoints_reduced_yellowpoints_1e5.jpg", height=25, width=25, res=600, pointsize=5, units="cm")
mu_hexplt2
#dev.off()

#Figure 5 in paper
P_draws_hex_list_reduced <- lapply(P_draws_hex_list, function(x) {x[c(1,3,5)]})
P_draws_hex_list_reduced <- P_draws_hex_list_reduced[c(1,3,5)]
dfp.hex.plt2 <- do.call(rbind, lapply(P_draws_hex_list_reduced, function(x) {do.call(rbind, x)}))
P_hexplt2 <- ggplot(dfp.hex.plt2) + geom_density_2d_filled(aes(x=p1,y=p2), bins = 11, alpha = 0.8, contour_var = "ndensity", adjust = 2, show.legend = FALSE) + scale_fill_brewer(type = "seq", palette = "RdBu", direction = -1) + geom_point(x=0.9,y=0.8, shape = 21, fill = "yellow", color="black", size = 1.5) + xlab(expression(P[list(1,1)])) + ylab(expression(P[list(2,2)])) + scale_x_continuous(sec.axis = sec_axis(~., name = expression(s[2]), breaks = NULL, labels = NULL)) + scale_y_continuous(sec.axis = sec_axis(~., name = expression(s[1]), breaks = NULL, labels = NULL)) + facet_grid(rows = vars(cat1), cols = vars(cat2))
#jpeg("P_contour_withpoints_reduced_yellowpoints_1e5.jpg", height=25, width=25, res=600, pointsize=5, units="cm")
P_hexplt2
#dev.off()

