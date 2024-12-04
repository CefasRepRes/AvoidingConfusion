
fixedC.gen <- function(alpha, K) {
  ####################################################################################################
  # Generates a confusion matrix such that each row is a draw from a symmetric Dirichlet distribution#
  ####################################################################################################
  #K:             Dimension of generated confusion matrix                                            #
  #alpha:        Vector of length K containing symmetric Dirichlet parameter for each row            #
  ####################################################################################################
  Conf <- matrix(0, K, K)
  for (i in 1:K){
    Conf[i,] <- rdirichlet(1, rep(alpha[i], K))
  }
  Conf
}


observedW.gen = function(mu, pmat, distribution="Poisson",
                         k, q) {
  #####################################################################
  # Generates a single realisation of observed data W.
  # Note that mu, k and distribution need to be vectors.
  # Note that the "VGAM" package must be loaded to use the zero-inflated Poisson distribution
  #
  # pmat is the correct probability matrix P
  # mu is the mean of the Poisson or Neg Bin distribution
  # k is the size parameter of the Neg Bin distribution
  #q is the zero-inflation parameter of the zero-inflated Poisson distribution
  # mu, k and distribution need to be vectors of the same length
  # distribution defines whether each generating process is Poisson,
  # negative binomial, or zero-inflated Poisson. It needs to be of the same length as mu
  # and k.
  #####################################################################
  
 # mu = mus
  #pmat=pmat
  #distribution = distribution
  #q=q
  
  nmu = length(mu)
  Yvec = rep(0, nmu)
  
  nclassified = matrix(0, nmu, nmu)
  
  nnbin <- 0
  for (i in 1:nmu) {
    if (distribution[i] == "Poisson") {
      Yvec[i] <- rpois(1, mu[i])
    }
    if (distribution[i] == "Negbin") {
      nnbin <- nnbin + 1
      Yvec[i] <- rnbinom(1, size = k[nnbin], mu = mu[i])
    }
    if (distribution[i] =="ZIPoisson") {
      Yvec[i] <- rzipois(1, mu[i], pstr0 = q[i])
    }
  }
  #    }
  
  for (i in 1:nmu) {
    nclassified[i, ] <- rmultinom(1, Yvec[i], pmat[i,])
  }
  
  Wvec = colSums(nclassified)
  
  list(Y=Yvec, W=Wvec, X=nclassified)
}



observedW.gen.multi = function(n, mu, pmat, distribution, k) {
##################################################################
# Generates more than 1 realisation of observed counts W
##################################################################
  
    Wmat = matrix(0, ncol=dim(pmat)[1], nrow=n)
  
  for (j in 1:n) {
    Wmat[j,] = observedW.gen(mu=mu, pmat=pmat,
                             distribution = distribution, k=k)$W
  }
  Wmat
}

observedW.gen.multi.covar = function(n, pmat, distribution, k, q, covs, intercepts, betas) {
##################################################################
# Generates multiple realisations of observed counts w in the case where covariates are included as a simple linear model with a single coefficient.
#
# covs is a vector of length n (i.e. a one-dimensional covariate observation for each location)
#intercepts is a vector of intercepts for each category
#betas is a vector of coefficients for each category
##################################################################
  
 # n = N
 # pmat = P_real
 # distribution = distribution
 # q = q_real
 # covs = covs
 # intercepts = intercepts_real
 # betas = betas_real
  
  Xmat = matrix(0, ncol=dim(pmat)[1], nrow=n)
  for (j in 1:n) {
    mus = sapply(1:dim(pmat)[1], function(k) {exp(intercepts[k] + betas[k]*covs[j])})
    Xmat[j,] = observedW.gen(mus, pmat=pmat,
                             distribution = distribution, q=q)$W
  }
  Xmat
}


Cgen <- function(N, P) {
  M <- length(N)
  C <- matrix(0, nrow = M, ncol = M)
  for (i in 1:M){
    C[i,] <- rmultinom(1, N[i], P[i,])
  }
  return(C)
}



geomprobs.gen <- function(Xdata, delta = 0.1, Conf, Unif = F, sig = 2) {
#############################################################################
##  Generates the probabilities used in a geometric generation to choose the
## 'jump' for updating the X classification matrix. The geometric
##  parameters are chosen such that the probability of the jump in a given
##  column exceeding the column sum divided by the number of classes is no
##  more than delta. We used a bound similar to Cantelli's inequality here (note
##  that this is only one choice, other options are possible).
##  The argument delta allows the user to tune the average jump size - 
##  larger delta will lead to smaller jumps on average.
#############################################################################
##  Xdata: Dataset to supply - columns must correspond to categories
##  delta: Upper bound (see above)
##  Conf: True confusion matrix to supply (optional)
##  Unif: If true, choose geometric parameter uniformly between delta and 0.5
##  sig: Number of decimal places to return for each probability
#############################################################################
  
  M <- dim(Xdata)[2]
  nreps <- dim(Xdata)[1]
  gprobs <- matrix(0, nrow = nreps, ncol = M)
  
  if (Unif == T) {
    
    if (missing(Conf)) {
      
      maxjump <- round(Xdata/M)
      
      for (i in 1:nreps) {
        for (j in 1:M) {
          gprobs[i,j] <- round(runif(1,min(((2*delta - 1) + sqrt((2*delta - 1)**2 + 4*delta*(maxjump[i,j] + 1)*(1-delta)))/(2*delta*(maxjump[i,j] + 1)),0.5),0.5),sig)
        }
      }
    }
    
    else {
      
      Cmax <- apply(Conf, 1, max)
      maxjump <- t(apply(Xdata, 1, function(x) {Cmax*x}))
      
      for (i in 1:nreps) {
        for (j in 1:M) {
          gprobs[i,j] <- round(runif(1,min(((2*delta - 1) + sqrt((2*delta - 1)**2 + 4*delta*(maxjump[i,j] + 1)*(1-delta)))/(2*delta*(maxjump[i,j] + 1)),0.5),0.5),sig)
        }
      }
      
    }
  }
  
  else {
    
    if (missing(Conf)) {
      
      maxjump <- round(Xdata/M)
      
      for (i in 1:nreps) {
        for (j in 1:M) {
          gprobs[i,j] <- round(min(((2*delta - 1) + sqrt((2*delta - 1)**2 + 4*delta*(maxjump[i,j] + 1)*(1-delta)))/(2*delta*(maxjump[i,j] + 1)),0.5),sig)
        }
      }
    }
    
    else {
      
      Cmax <- apply(Conf, 1, max)
      maxjump <- t(apply(Xdata, 1, function(x) {Cmax*x}))
      
      for (i in 1:nreps) {
        for (j in 1:M) {
          gprobs[i,j] <- round(min(((2*delta - 1) + sqrt((2*delta - 1)**2 + 4*delta*(maxjump[i,j] + 1)*(1-delta)))/(2*delta*(maxjump[i,j] + 1)),0.5),sig)
        }
      }
      
    }
  }
  
  return(gprobs)
}


Xinitial = function(Xsin, Pin) {
####################################################################
## Creates an initial X matrix from observed W values. This is used
## in the MCMC algorithm
####################################################################
  M = length(Xsin)
  X = matrix(0, ncol=M, nrow=M)
  for (j in 1:M) {
    X[,j] = rmultinom(1, Xsin[j], Pin[,j])
  }
  X
}


knewprior <- function(x,sigma,p){
  ###############################################################
  # Calculates log-likelihood for prior on the ratio mu/k.
  # Used in the mixed Poisson - negative binomial example
  ################################################################
  if(x < 0){
    return(log(p) + log(2) + dnorm(x,0,sigma,log=T))
  }else{
    cauchy_sigma <- (dnorm(0)/dcauchy(0)) * 1/sigma * (p/(1-p))
    return(log(1- p) + log(2) + dcauchy(x,scale=1/cauchy_sigma, log=T))
  }
}


dlike_log_NB_Po <- function(Xlik, mulik, klik, Wslik, Plik, NB_inds) {
  ##################################################################
  # Log-likelihood for a single X matrix given mixed Poisson - 
  # negative binomial component distributions
  ##################################################################
  
  if (any(colSums(Xlik)!=Wslik)){
    return(-Inf)
  }
  M = length(Wslik)
  Mk = length(NB_inds)
  ret <- 0
  for (i in 1:M){
    ret <- ret + dmultinom(Xlik[i,], sum(Xlik[i,]), Plik[i,], log = T)
  }
  if (Mk == 1 & M >= 3) {
    ret <- ret + sum(dpois(rowSums(Xlik[-NB_inds,]),lambda=mulik[-NB_inds],log=T)) + dnbinom(sum(Xlik[NB_inds,]),size=klik, mu=mulik[NB_inds],log=T)
  }
  if (Mk >= 2 & M == (Mk+1)) {
    ret <- ret + dpois(sum(Xlik[-NB_inds,]),lambda=mulik[-NB_inds],log=T) + sum(dnbinom(rowSums(Xlik[NB_inds,]),size=klik, mu=mulik[NB_inds],log=T))
  }
  if (Mk >=2 & M >= (Mk+2)) {
    ret <- ret + sum(dpois(rowSums(Xlik[-NB_inds,]),lambda=mulik[-NB_inds],log=T)) + sum(dnbinom(rowSums(Xlik[NB_inds,]),size=klik, mu=mulik[NB_inds],log=T))
  }
  else {
    ret <- ret + dpois(sum(Xlik[-NB_inds,]),lambda=mulik[-NB_inds],log=T) + dnbinom(sum(Xlik[NB_inds,]),size=klik, mu=mulik[NB_inds],log=T)
  }
  return(ret)
  
}


dlike_log_NB_Po_rep <- function(Xmatslik, mulik, klik, Wobslik, Plik, NB_inds){
  #######################################################################
  ## Mixed Poisson and negative binomial log-likelihood for multiple reps
  #######################################################################
  
  R = dim(Xmatslik)[3]
  loglik = 0
  for (j in 1:R) {
    loglik = loglik + dlike_log_NB_Po(Xmatslik[,,j], mulik, klik, Wobslik[j,], Plik, NB_inds)
  }
  return(loglik)
}

mcmcML_NB_Po2 <- function(Wdata, C_matrix, distributions = c("Poisson", "Poisson", "Negbin"), Nits=50000, thin=5, geom.probs,
                          cand_sd_mu, cand_sd_k, log_prop_mu, log_prop_k, chooseX=1:nrow(Wdata), fit_init = F, C_sam=T) {
  ###########################################################################
  # MCMC for mixed Poisson and negative binomial count distributions
  ###########################################################################
  # Wdata : N by M array of observed data
  # C_matrix : Confusion matrix
  # distributions: Character or character vector of distribution names (each 
  #                entry either "Poisson" or "Negbin") to indicate which rows
  #                of C_matrix correspond to which family of distributions.
  #                If a vector, length must match number of rows of C_matrix.
  # Nits : Number of MCMC iterations
  # thin : save X, mu and k every "thinth" iteration. Must divide exactly
  #        into Nits.
  # geom.probs : Probability for each variable (columns) and each set of
  #              observations (rows) for generating geometric proposal random
  #              variable for jittering X matrices
  # cand_sd_mu : Standard deviation of Normal proposal distribution around
  #              current value of mu
  # cand_sd_k  : Standard deviation of Normal proposal distribution around
  #              current value of k
  # log_prop_mu : Whether mu should be generated on log scale (T) or orginal
  #               scale (F)
  # chooseX    : The numbers in this vector denote which X variable matrices
  #              should be saved for investigation after the fit 
  # fit_init   : Whether the intial MCMC parameters should be estimated (T) by
  #              method of moments or chosen randomly (F)
  # C_sam     : Whether the C matrix should be sampled or fixed throughout
  #             the MCMC run
  ########################################################################
  
  nreps <- dim(Wdata)[1]               # Number of sets of observations
  M <- dim(Wdata)[2]                   # Number of variables
  NB_inds <- which(distributions == "Negbin") # Indices of negative binomial categories
  Mk <- length(NB_inds)
  
  ## Create the initial P matrices and X matrices (one per rep)
  
  if(C_sam == T){
    alpha_C <- C_matrix + 1
    P <- rdirichlet(M,alpha_C)
    
  }else{
    P <- C_matrix / rowSums(C_matrix)
  }
  
  Xmats <- array(0,dim=c(M, M, nreps))
  
  for (j in 1:nreps) {
    Xmats[,,j] = Xinitial(Xsin=as.numeric(Wdata[j,]), Pin=P)
  }
  Xmatsup = Xmats
  
  ## Initialise the counters for updating X and mu
  
  countmuk = rep(0, Mk)
  countX = rep(0, nreps)
  
  ### Initialise mu and k
  
  if (fit_init == F) {
    
    munow <- exp(runif(M, -2, 2))
    know <- exp(runif(Mk, -2, 2))
    
  }
  
  else {
    
    munow <- rep(0, M)
    know <- rep(0, Mk)
    
    for (j in 1:Mk) {
      tempj = fitdistr(Wdata[,NB_inds[j]], "Negative Binomial")
      munow[NB_inds[j]] = tempj$estimate[2]
      know[j] = tempj$estimate[1]
    }
    
    if (Mk <= (M-2)) {
      munow[-NB_inds] <- colMeans(Wdata[,-NB_inds])
    }
    if (Mk == (M-1)) {
      munow[-NB_inds] <- mean(Wdata[,-NB_inds])
    }
    
  }
  
  ## Make sure thin divides exactly into N
  
  Nthin = Nits / thin
  
  ## Arrays for saving generated posterior values
  
  save_mu <- array(0, dim=c(Nthin, M))
  save_k <- array(0, dim = c(Nthin, Mk))
  
  nsaveX = length(chooseX)
  save_X = array(0, dim=c(M,M,Nthin,nsaveX))
  
  save_P = array(0,dim=c(M,M,Nthin))
  
  ## Start of MCMC loop
  
  for (i in 1:Nits){
    
    ## update each X matrix
    
    for (j in 1:nreps) {
      
      temp = updateX_NB_Po(Xup = Xmats[,,j], Wsup=as.numeric(Wdata[j,]),
                           munowup = munow, knowup = know, Pup=P,
                           counter = countX[j], gprobs=geom.probs[j,], NB_inds = NB_inds)
      Xmatsup[,,j] = temp$X
      countX[j] = temp$count
    }
    Xmats = Xmatsup
    
    # Update mu and k together
    temp = updatemuk_NB_Po(Xmatsup=Xmats, knowup = know, munowup=munow, Wobs=Wdata, P=P,
                           log_prop_mu=log_prop_mu,log_prop_k=log_prop_k,
                           cand_sd_k=cand_sd_k,
                           cand_sd_mu=cand_sd_mu,
                           countmuk = countmuk, NB_inds = NB_inds, nreps)
    munow = temp$mu
    know = temp$k
    countmuk = temp$count
    
    # update C
    
    if (C_sam == T){
      tmpC <- alpha_C + apply(Xmats,c(1,2),sum)
      P <- rdirichlet(M, tmpC)
    }
    
    if (i %% thin == 0) {
      ii = as.integer(i/thin)
      save_mu[ii,] = munow
      save_k[ii,] = know
      save_P[,,ii] = P
      for (s in 1:nsaveX) {
        save_X[,,ii,s] <- Xmats[,,chooseX[s]]
      }
    }
  } # end of (i in 1:Nits) loop
  
  Xprob = 100 * countX / (M*Nits)
  mukprob = 100 * countmuk/ (Mk*Nits)
  
  output = list(mu=save_mu, k= save_k, X=save_X, P=save_P,
                Xprob=Xprob, mukprob = mukprob, chooseX = chooseX)
  output
}


updatemuk_NB_Po <- function(Xmatsup, munowup, knowup, Wobs, P, log_prop_mu = log_prop_mu,
                            log_prop_k=log_prop_k,
                            cand_sd_mu = cand_sd_mu, cand_sd_k=cand_sd_k, countmuk = countmuk, NB_inds, nreps) {
  ################################################################
  ## Update mu and k together
  ################################################################
  
  M = length(munowup)
  Mk = length(NB_inds)
  
  if (M > Mk) {
    if (M >= (Mk + 2)) {
      if (nreps > 1) {
        Yvec = apply(Xmatsup[-NB_inds,,], c(3,1), sum)
        Ymeans = colMeans(Yvec)
      }
      else {
        Yvec = apply(Xmatsup[-NB_inds,,], 1, sum)
        Ymeans = Yvec
      }
    }
    else {
      if (nreps > 1) {
        Yvec = apply(Xmatsup[-NB_inds,,], 2, sum)
      }
      else {
        Yvec = Xmatsup[-NB_inds,,]
      }
      Ymeans = mean(Yvec)
    }
    
    munowup[-NB_inds] = rgamma((M-Mk), shape=1 + nreps*Ymeans, rate=0.0001 + nreps)
    
  }
  
  munowup_ = munowup
  knowup_ = knowup
  
  for (j in 1:Mk) {
    
    if (log_prop_mu[j]==FALSE) {
      muj_ <- rnorm(1, mean=munowup[NB_inds[j]] , sd=cand_sd_mu[j])
    }
    else {
      muj_ <- exp(rnorm(1, mean=log(munowup[NB_inds[j]]) , sd=cand_sd_mu[j]))
    }
    if (log_prop_k[j]==FALSE) {
      kj_ <- rnorm(1, mean=knowup[j] , sd=cand_sd_k[j])
    }
    else {
      kj_ <- exp(rnorm(1, mean=log(knowup[j]) , sd=cand_sd_k[j]))
    }
    
    if (muj_ > 0 & kj_ > 0) {
      
      munowup_[NB_inds[j]] = muj_
      knowup_[j] = kj_
      log_like = dlike_log_NB_Po_rep(Xmatsup, munowup, knowup, Wobs, P, NB_inds)
      log_like_ = dlike_log_NB_Po_rep(Xmatsup, munowup_, knowup_, Wobs, P, NB_inds)
      muprior = dgamma(munowup[NB_inds[j]], shape=1, rate=0.0001, log=T)
      muprior_ = dgamma(muj_, shape=1, rate=0.0001, log=T)
      kprior = knewprior(log(munowup[j]/knowup[j]), 0.5, 0.1)
      kprior_ = knewprior(log(muj_/knowup_[j]), 0.5, 0.1)
      
      alpha <- log_like_ - log_like + muprior_ - muprior + kprior_ -kprior
      
      if (log(runif(1)) < alpha){
        munowup = munowup_
        knowup = knowup_
        countmuk[j] = countmuk[j] + 1
      }
      else {
        
        munowup_ = munowup
        knowup_ = knowup
        
      }
      
    }
  }
  list(mu = munowup, k=knowup, count=countmuk)
}

updateX_NB_Po = function(Xup, Wsup, munowup, knowup, Pup = P, counter,
                         gprobs, NB_inds) {
#############################################################################
## Update the X matrix.
#############################################################################
  M = dim(Xup)[1]
  
  for(j in 1:M){
    X_ <- Xup
    limit = rgeom(1, gprobs[j]) + 1
    n0 <- which(Xup[,j] >= limit)
    if (length(n0) >= 1) {
      
      if (length(n0) > 1){
        minus_1 <- sample(n0,size=1)
      }else{
        minus_1 <- n0
      }
      
      X_[minus_1,j] <- X_[minus_1,j] -  limit
      plus_1 <- sample((1:M)[-minus_1],size=1)
      X_[plus_1,j] <- X_[plus_1,j] +  limit
      
      log_like <- dlike_log_NB_Po(Xup, mulik = munowup, klik = knowup, Wsup, Pup, NB_inds)
      log_like_ <- dlike_log_NB_Po(X_, mulik = munowup, klik= knowup, Wsup, Pup, NB_inds)
      hast_ratio <- log(length(n0)) - log(sum(X_[,j]>=limit))
      
      alpha <- log_like_ - log_like + hast_ratio 
      
      if (log(runif(1)) < alpha){
        Xup <- X_
        counter = counter + 1
      }
    }       # end of "if n0>1 loop"
  }         # end of "i=1 to N loop"
  list(X=Xup, count=counter)
}




