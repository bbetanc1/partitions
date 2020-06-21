#' @useDynLib partitions
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rbinom rgamma runif

NULL

## Wrapper function ##
# Inputs: N - the sample size based on the dataset           #
#         Prior - Specify Prior Type: ESCNB, DP,PY           #
#         param - Parameter for prior type specified         #
# Output: Parameters for probability of existing cluster (A) #
#         and new cluster (B) for reseating algorithms       #
SetParam <- function(Prior, param, N){
    if (Prior == "DP"){
    	theta <- param[1]
    	A <- seq(N)
		  B <- rep(theta,N)
    }
    if (Prior == "PY"){
    	theta <- param[1]
    	delta <- param[2]
    	A <- seq(N) - delta
    	B <- rep(theta, N) + seq(N) * delta
    }
    if (Prior == "ESCNB"){
        rnb <- param[1]
        pnb <- param[2]
        beta <- ((1-pnb)^rnb)/(1-((1-pnb)^rnb))
        A <- (seq(N) + rnb)
        B <- (seq(N) + 1) * beta * rnb
    }
    return(list(A=A, B=B))
}

SetParamESCD <- function(mus, N, M){
    A <- (seq(N) + 1) * (c((mus[-1])/mus[-(M+1)], rep(0, N-M)))
    B <- (seq(N) + 1) * mus[1]
    return(list(A=A, B=B))
}


## Wrapper function ##
# Inputs: Prior - Specify Prior Type: DP, PY, ESCNB, ESCD      #
# Output: Lower and Upper bound on support of the distribution #
#         of Prior parameters                                  #

lowupSlice <- function(Prior){
    if (Prior == "DP"){
    	lo <- 0
    	up <- Inf
    }
    if (Prior == "PY"){
    	lo <- c(0, 0)         
    	up <- c(Inf, 1) 
    }
    if (Prior == "ESCNB"){
        lo <- c(0, 0)
        up <- c(Inf, 1)
    }
    if (Prior == "ESCD"){
        lo <- c(0, 0, 0)
        up <- c(Inf, Inf, 1)
    }
   return(list(lo=lo, up=up))
}

## Wrapper function ##
# Inputs: Prior - Specify Prior Type: DP, PY, ESCNB, ESCD       #
#         N - number of record in data                          #
# Output: Initial values for hyperpameters                      #

SetInit <- function(Prior, N){
	if (Prior == "DP"){
    	#concentration parameter
        theta <- N/2
        x1 <- as.vector(theta)
      }  
    if (Prior == "PY"){
    	#concentration parameter
        theta <- N/2 
        # discount parameter
        beta <- 0.5 
        x1 <- c(theta, beta)
      }
     if (Prior == "ESCNB"){
        rnb <- 1
        pnb <- 0.5
        x1 <- c(rnb, pnb)
       }
     if (Prior == "ESCD"){
        alpha <- 1
        rnb <- 1
        pnb <- 0.5
        x1 <- c(alpha, rnb, pnb)
       }
       return(x1)
}

## Wrapper function ##
# Inputs: Prior - Specify Prior Type: DP, PY, , ESCNB, ESCD #
# Output: Indicator for hyperparameter sampling             #

SamInd <- function(Prior){
    if (Prior == "DP"){
        #concentration parameter
        samind <- c(1)
    }
    if (Prior == "PY"){
        # concentration parameter
        # discount parameter
        samind <- c(1, 1)
    }
    if (Prior == "ESCNB"){
        # NB parameters r and p
        samind <- c(1, 1)
    }
    if (Prior == "ESCD"){
        # NB parameters r and p, alpha fixed
        samind <- c(0, 1, 1)
    }
    return(samind)
}


## Wrapper function ##
# Inputs: Prior - Specify Prior Type: DP, PY, ESCNB, ESCD       #
#         N - number of record in data                          #
# Output: Vector of hyperpameter values for prior on parameters #
#         of specified partition prior                          #

hpriorpar <- function(Prior, N){
	 if (Prior == "DP"){
    	# Gamma(ag1, bg1) prior over concentration parameter
        ag <- 1
        bg <- 1/(N/2)
        hpriorpar <- c(ag, bg)
      }  
    if (Prior == "PY"){
    	# Gamma(ag, bg) prior over concentration parameter
        bg <- 1/(N/2)
        ag <- 1
      # Beta(ab,bb) prior over discount parameter
        ab <- 1
        bb <- 1
        hpriorpar <- c(ag, bg, ab, bb)
      }
     if (Prior == "ESCNB" || Prior == "ESCD"){
        # Gamma(ag, bg) prior over r
        bg <- 1
        ag <- 1
        # Beta(ab,bb) prior over  p
        ab <- 2
        bb <- 2
        hpriorpar <- c(ag, bg, ab, bb)
       }
    return(hpriorpar)
}

## Wrapper function ##
# Inputs: Prior - Specify Prior Type: DP, PY, NBNB, NBD      #
# Output: Parameter names used to save posterior samples     #

NameParam <- function(Prior, nfields){
    if (Prior=="DP"){
        cnames <- c("theta")
        for(i in 1:nfields){
            cnames <- c(cnames, sprintf("beta_%d",i))
        }
    }
    if (Prior=="PY"){
        cnames <- c("theta", "delta")
        for(i in 1:nfields){
            cnames <- c(cnames, sprintf("beta_%d",i))
        }
    }
    if (Prior=="ESCNB"){
        cnames <- c("r", "p")
        for(i in 1:nfields){
            cnames <- c(cnames, sprintf("beta_%d",i))
        }
    }
    if (Prior=="ESCD"){
        cnames <- c("alpha", "r", "p")
        for(i in 1:nfields){
            cnames <- c(cnames, sprintf("beta_%d",i))
        }
    }
    return(cnames)
}

# Computes proportion of each values in a field
# Used to compute empirical distribution of data
calcProp <- function(i,data) {
    return(table(data[,i])/sum(table(data[,i])))
}

# For each field, remap values to a consecutive 
# list of integers starting with 1
remap <- function(i, data) {
    return(as.numeric(factor(data[,i], labels = seq(1,length(unique(data[,i]))))))
}

# Sample from a Dirichlet
rdirichletf <- function(params) {
    p <- sapply(params, function(x) rgamma(1,x))
    return(p/sum(p))
}

# Find Beta distribution parameters for distortion
# probailities, given the mean and sd
abeta <- function(bmean, bsd){
    return (bmean*(1-bmean)*(bmean/bsd^2) - bmean)
}


# selects a subset of points (called "block") that agree on some
# randomly choosen features. Chaperones are chosen within the block 
# for a few consecutive iterations.
select_block<-function(x, size_block){
  N<-dim(x)[1] #  number of datapoints
  L<-dim(x)[2] #  number of fields
  n_fields <- rbinom(1, L, prob=0.75) # select how many fields the chaperones will agree on
  block<-c(1:N)
  if(n_fields>0){
    fields_to_agree<-sample(c(1:L),size = n_fields,replace = FALSE)# select which fields the chaperones will agree on
    x0<-x[sample(c(1:N),1),] #select a point at random (to be used as "pivot" point)
    for(l in 1:n_fields){
      if(length(block)>size_block){
        block<-block[which(x[block,fields_to_agree[l]]==x0[fields_to_agree[l]])]
      }
    }
  }
  return(block)
}

#' Remap data to a list of consecutive integers per field
#'
#' @param data Data frame containing only categorical variables 
#' @return Data frame of remapped values to consecutive 
#'         list of integers
#' @export
DataRemap <- function(data) {
	nfields <- ncol(data)
	x <- sapply(1:nfields,remap,data=data)
	return(x)
}

## Function Details ##
# Main sampler function: chaperones algorithm for cluster    #
# assignments and Slice sampler for priors' parameters and   #
# distortion probablities of the fields.                     #
#                                                            #
# Inputs: Data - e.g "Italy"                                 #
#         Prior - specify Prior Type: DP,PY,ESCNB, ESCD      #
#         burn - MCMC burn-in period                         #
#         nsamples - MCMC iterations                         #
#         truncds - truncation value for field distortions   #
#         samds - indicator to sample field distortions      #
#         spacing - chaperones restricted Gibbs moves        #
#         w - step size for slice sampler                    #
#         m - number of steps                                #
#         hpriorpar - vector of prior hyperparameters        #
#         block_flag - TRUE for non-uniform chaperones       #
# Output: Posterior samples for prior parameters and         #
#         cluster assingments                                #
# Summary: This is a wrapper function that requires the user #
#          specify the data, prior and MCMC parameters       #
#          and returns  posterior samples of cluster         #
#          assignments and hyperparameters.                  #

#' Posterior samples of cluster assignments and Prior parameters
#'
#' @param data Data frame containing only categorical variables 
#' @param Prior Specify partition prior: "DP", "PY", "ESCNB"
#' @param burn MCMC burn-in period
#' @param nsamples MCMC iterations after burn-in
#' @param spacing Thinning for chaperones algorithm (default 1000)
#' @param block_flag TRUE for non-uniform chaperones (default)
#' @return List with posterior samples for cluster assignments (Z),
#'         Prior parameters and distortion probabilities (Params)
#' @export
SampleCluster <- function(data, Prior, burn, nsamples, spacing=1000, block_flag=TRUE){
    x <- DataRemap(data)
    N <- nrow(x) 
    nfields <- ncol(x)
    proportions <- lapply(1:nfields,calcProp,data=x)
    # Initial clustering
    z0 <- sample(1:(N), N, replace=TRUE)
    # No gaps
    z <- as.numeric(factor(z0 , labels = seq(1,length(unique(z0)))))
    # Initial Prior parameter values
    x1 <- SetInit(Prior, N)
    samind <- SamInd(Prior)
    # Beta prior for distortion parameters
    bmean <- 0.005
    bsd <- 0.01
    cb <- abeta(bmean, bsd)
    db <- cb/bmean - cb
    hpriords <- c(cb, db)
    pdist <- cb/(cb+db)
    # Initial value for distortion parameters
    betas <- rep(pdist,nfields)
    truncds <- 0.1
    
    # hyperparameter values for priors on Prior parameters
    hpriorpar <- hpriorpar(Prior, N)
    
    # lo - lower limit distribution support (slice)
    # up - upper limit distribution support (slice)
    loup <- lowupSlice(Prior)
    lo <- loup$lo
    up <- loup$up
    
    s <- matrix(0, nsamples, length(x1) + nfields)
    zM <- matrix(0,  nsamples, N)
    
    lx <- loglikxSP(betas, x, z, proportions)
    ll <- 1
    # number of iterations for each random block of chaperons
    itbl <- 50
    
    if(Prior=="DP"|Prior=="PY"|Prior=="ESCNB"){
      for (i in 1:(burn + nsamples)){
        if ((i %% 10)==0){
          print(paste("iter=",i))
        }
        Khat <- length(unique(z))
        Nk <- as.vector(table(z))
        # univariate slice sampler for all parameters
        w <- 1 # step size for slice sampler
        m <- 10 # limit on steps for slice sampler
        x1 <- unislicem(x1, N, Khat,lx, Nk, hpriorpar, w, m, lo, up, Prior, samind)
        lods <- 0
        upds <- truncds
        betas <- unislicespb1(betas, x, z, proportions, hpriords, w, m, lods, upds,
                              x1, N, Khat, Nk, hpriorpar, Prior)
        
        # parameters of reallocation algortihm #
        AB <- SetParam(Prior, x1, N)
        A <- AB$A
        B <- AB$B
        if(block_flag){
          for(ss in 1:(spacing/itbl)){    
            bl <- select_block(x, 200)        
            z1 <- Web_SamplerSP_fbl(x, z, bl, A, B, betas, proportions, 1, itbl)
            # no gaps in cluster assigment #
            z <- as.numeric(factor(z1[1,], labels = seq(1,length(unique(z1[1,])))))
          }
        }else{
          z1 <- Web_SamplerSP(x, z, A, B, betas, proportions, 1, spacing)
          # no gaps in cluster assigment #
          z <- as.numeric(factor(z1[1,], labels = seq(1,length(unique(z1[1,])))))
        }
        lx <- loglikxSP(betas, x, z, proportions)
        #loglik[i] <- lx
        if((i > burn) ){
          zM[ll, ] <- z
          s[ll,] <- c(x1,betas)
          ll=ll+1
        }
      }
    }      
    
    if(Prior=="ESCD"){
      M <- 100 # max cluster size
      Khat <- length(unique(z))
      tz <- table(z)
      Nk <- as.vector(tz)
      if(max(tz) >= M){
        M <- max(tz)
      }
      Lm <- table(factor(tz, levels=1:M))
      alpha0 = x1[1];
      r0 = x1[2];
      p0 = x1[3];
      lmu0 = lgamma(seq(M) + r0) + (seq(M))*log(p0) + r0*log(1-p0) - log(1-(1-p0)^r0)-
        lgamma(r0) - lgamma(seq(M) + 1)
      mu0 = exp(lmu0)
      for (i in 1:(burn + nsamples)) {
        if ((i %% 10)==0){
          print(paste("iter=",i))
        }
        # univariate slice sampler for all parameters
        w <- 1 # step size for slice sampler
        m <- 10 # limit on steps for slice sampler
        x1 <- unislicemESCD(x1, lx, Lm, mu0, hpriorpar, w, m, lo, up, samind)
        lods <- 0
        upds <- truncds
        betas <- unislicespb1(betas, x, z, proportions, hpriords, w, m, lods, upds,
                                x1, N, Khat, Nk, hpriorpar, Prior)
  
        # parameters of reseating algortihm 
        alpha0 <- x1[1]
        r0 <- x1[2]
        p0 <- x1[3]
        
        # Sampling mu from dirichlet of size max cluster size
        lmu0 <-lgamma(seq(M) + r0) + (seq(M))*log(p0) + r0*log(1-p0) - log(1-(1-p0)^r0)-
          lgamma(r0) - lgamma(seq(M) + 1)
        mu0 <- exp(lmu0)
        alphasDM <- c(alpha0 * mu0 + Lm, abs(alpha0 * (1 - sum(mu0))))
        mus <- rdirichletf(alphasDM)
        AB <- SetParamESCD(mus, N, M)
        A <- AB$A
        B <- AB$B
        if(block_flag){
          for(ss in 1:(spacing/itbl)){    
            bl <- select_block(x, 200)        
            z1 <- Web_SamplerSP_fbl(x, z, bl, A, B, betas, proportions, 1, itbl)
            # no gaps in cluster assigment #
            z <- as.numeric(factor(z1[1,], labels = seq(1,length(unique(z1[1,])))))
          }
        }else{
          z1 <- Web_SamplerSP(x, z, A, B, betas, proportions, 1, spacing)
          # no gaps in cluster assigment #
          z <- as.numeric(factor(z1[1,], labels = seq(1,length(unique(z1[1,])))))
        }
        Khat <- length(unique(z))
        tz <- table(z)
        Nk <- as.vector(tz)
        if (max(tz) >= M) {
          M <- max(tz)
        }
        Lm <- table(factor(tz, levels = 1:M))
        lx <- loglikxSP(betas, x, z, proportions)
        if ((i > burn)) {
          zM[ll,] <- z # cluster samples
          s[ll, ] <- c(x1, betas) # parameters samples
          ll = ll + 1
        }
      }
    }  
    colnames(s) <- NameParam(Prior, nfields)
    return(list(Z=zM,Params=s))
}