#Steorts, R. 2015-03-03 (From Steorts (2014))

##### FUNCTIONS TO ANALYZE GIBBS SAMPLER OUTPUT #####

#' Calculates FNR and FDR when ground truth is available
#'
#' @param z Matrix with posterior samples
#' @param id Data identifiers (ground truth)
#' @return Two component vector with FNR and FDR
#' @export
calcError <- function(z, id){
    est.links <- links(z)
    est.links.pair <- pairwise(est.links)
    true.links <- links(matrix(id,nrow=1))
    true.links.pair <- pairwise(true.links)
  # Correct, incorrect, and missing links
    comparison <- links.compare(est.links.pair,true.links.pair,counts.only=TRUE)
    missing.links <- comparison$missing
    true.links<-comparison$correct
    false.links <- comparison$incorrect
    fnr <- missing.links/(true.links+missing.links)
    if(true.links+false.links==0){
        fdr <- 0
    }else{
        fdr <- false.links/(true.links+false.links)
    }
    output <- c(fnr,fdr)
    names(output) <- c("FNR","FDR")
    return(output)
}

#' @export
estLinks <- function(z){
est.links <- links(z)
return(est.links)
}

#' @export
StatsCluster <- function(z){
    K <- apply(z, 1, function(x) {length(unique(x))})
    if(length(unique(K))==1){
        Mnk <- apply(z, 1, table)
        Nks <- split(Mnk, col(Mnk))
    }else{
        Nks <- apply(z, 1, table)
    }
    maxNk <- unlist(lapply(Nks ,max)) # Maximum cluster size #
    meanNk <- unlist(lapply(Nks ,mean)) # Mean cluster size #
    singles <- unlist(lapply(Nks,function(nk)sum(nk==1))) # Number of singleton clusters #
    # Quantiles 90% and 95%
    q90 <- unlist(lapply(Nks,quantile,prob=0.9))
    q95 <- unlist(lapply(Nks,quantile,prob=0.95))
    stats <- cbind(K, maxNk, meanNk, singles, q90, q95)
    return(stats)
}

# Function to compute a record's MMS based on a single linkage structure
mms <- function(lambda,record){
	lambda. <- drop(lambda)
	return(which(lambda.==lambda.[record]))
}

# Function to compute a record's MPMMS based on a Gibbs sample
# Note: It returns a list of the MPMMS ($mpmms) and its probability ($prob)
mpmms <- function(lam.gs,record){
	G <- dim(lam.gs)[1]
	mms.out <- NULL
	num.mms <- 0
	for(g in 1:G){
		lambda.g <- drop(lam.gs[g,])
		mms.g <- which(lambda.g==lambda.g[record])
		mms.g.str <- paste(mms.g,collapse="-")
		if(is.null(mms.out[as.character(mms.g.str)]) || 
			is.na(mms.out[as.character(mms.g.str)])){
			num.mms <- num.mms+1
			mms.out <- append(mms.out,1)
			names(mms.out)[num.mms] <- mms.g.str
		}else{
			mms.out[as.character(mms.g.str)] <-
				mms.out[as.character(mms.g.str)]+1
		}
	}
	mms.out <- mms.out/G
	mpmms.info <- mms.out[which.max(mms.out)]
	mpmms.out <- as.numeric(strsplit(names(mpmms.info),
					split="-",fixed=TRUE)[[1]])
	prob.out <- as.numeric(mpmms.info)
	return(list(mpmms=mpmms.out,prob=prob.out))
}




# Function that returns the shared MPMMS (except with an easier condition to code than JASA paper)
# Function to make a list of vectors of estimated links by "P(MPMMS)>0.5" method
# Note: The default settings return only MPMMSs with multiple members
links <- function(lam.gs,include.singles=FALSE,show.as.multiple=FALSE){
	N <- dim(lam.gs)[2]
	links.out <- list()
	for(r in 1:N){
		mpmms.list.r <- mpmms(lam.gs,r)
		mpmms.r <- mpmms.list.r$mpmms
		if(length(mpmms.r)>1 || include.singles){
			prob.r <- mpmms.list.r$prob
			if(prob.r>0.5){
				if(r==min(mpmms.r) || show.as.multiple){
					links.out[[length(links.out)+1]] <- mpmms.r
				}
			}else{
				if(include.singles){links.out[[length(links.out)+1]] <- r}
			}
		}
	}
	return(links.out)
}


# Function to take links list that may contain 3-way, 4-way, etc.
# and reduce it to pairwise only (e.g., a 3-way link 12-45-78 is
# changed to 2-way links: 12-45, 12-78, 45-78
pairwise <- function(.links){
	links.out <- list()
	if(length(.links)>0){
		# more efficient to use combn
		for(l in 1:length(.links)){
			link.l <- .links[[l]]
			for(r1 in 1:(length(link.l)-1)){
				for(r2 in (r1+1):length(link.l)){
					links.out[[length(links.out)+1]] <- c(link.l[r1],link.l[r2])
				}
			}
		}
	}
	return(links.out)
}

# This function takes a set of pairwise links and identifies correct, incorrect,
# and missing links
# (correct = estimated and true, incorrect = estimated but not true,
# missing = true but not estimated)
   # TODO: consider using %in% to remove for loops?
links.compare <- function(est.links.pair,true.links.pair,counts.only=TRUE){
	correct.out <- list()
	incorrect.out <- list()
	missing.out <- list()
	if(length(est.links.pair)>0 && length(true.links.pair)>0){
		for(l in 1:length(est.links.pair)){
			link.l <- est.links.pair[[l]]
			found <- FALSE
			for(t in 1:length(true.links.pair)){
				if(all(link.l==true.links.pair[[t]])){
					correct.out[[length(correct.out)+1]] <- link.l
					found <- TRUE
				}
			}
			if(!found){
				incorrect.out[[length(incorrect.out)+1]] <- link.l
			}
		}
		for(m in 1:length(true.links.pair)){
			link.m <- true.links.pair[[m]]
			missing <- TRUE
			for(e in 1:length(est.links.pair)){
				if(all(link.m==est.links.pair[[e]])){
					missing <- FALSE
				}
			}
			if(missing){
				missing.out[[length(missing.out)+1]] <- link.m
			}
		}
	}
	if(length(est.links.pair)==0 && length(true.links.pair)>0){
		missing.out <- true.links.pair
	}
	if(length(est.links.pair)>0 && length(true.links.pair)==0){
		incorrect.out <- est.links.pair
	}
	if(counts.only){
		return(list(correct=length(correct.out),incorrect=length(incorrect.out),
			missing=length(missing.out)))
	}else{
		return(list(correct=correct.out,incorrect=incorrect.out,
			missing=missing.out))
	}
}

# Check whether 2 records which are estimated to be linked have the same IDs
check_IDs <- function(recpair, identity_vector) {
	rec1 <- recpair[1]
	rec2 <- recpair[2]
	return(identity_vector[rec1] == identity_vector[rec2])
}


