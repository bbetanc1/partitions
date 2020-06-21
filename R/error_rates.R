## functions to compute FNR and FDR ##

#' Calculates FDR when ground truth is available
#'
#' @param z Vector of cluster assignments
#' @param id Vector of true cluster assignments (ground truth)
#' @return FDR
#' @export
fdr_fun <- function(z,id){
  if(n_matches_fun(z)==0){
    return(0)
  }else{
    return(
      sum(vapply(X = c(1:length(z)), FUN = i_false_det_fun, FUN.VALUE = 1, z.=z,id.=id))/(2*n_matches_fun(z))
    )
  }
}

#' Calculates FNR when ground truth is available
#'
#' @param z Vector of cluster assignments
#' @param id Vector of true cluster assignments (ground truth)
#' @return FNR
#' @export
fnr_fun <- function(z,id){
  if(n_matches_fun(id)==0){
    return(0)
  }else{
    return(
      sum(vapply(X = c(1:length(z)),FUN = i_false_neg_fun,FUN.VALUE = 1,z.=z,id.=id))/(2*n_matches_fun(id))
    )
  }
}

#' Calculates average FDR when ground truth is available
#'
#' @param zm Matrix with posterior samples of cluster assignments
#' @param id Vector of true cluster assignments (ground truth)
#' @return Average FDR over posterior samples
#' @export
mean_fdr <- function(zm,id){
  fdr_vec <- apply(X = zm, MARGIN = 1, FUN = fdr_fun, id=id)
  return(mean(fdr_vec))
}

#' Calculates average FNR when ground truth is available
#'
#' @param zm Matrix with posterior samples of cluster assignments
#' @param id Vector of true cluster assignments (ground truth)
#' @return Average FNR over posterior samples
#' @export
mean_fnr <- function(zm,id){
  fnr_vec <- apply(X = zm, MARGIN = 1, FUN = fnr_fun, id=id)
  return(mean(fnr_vec))
}

# auxiliary functions to compute fnr and fdr

i_false_det_fun <- function(i,z.,id.){
  i_match <- which(z.[-i]==z.[i])
  i_match_true <- which(id.[-i]==id.[i])
  i_common <- length(intersect(i_match,i_match_true))
  i_false_det <- length(i_match)-i_common
  return(i_false_det)
}

i_false_neg_fun <- function(i,z.,id.){
  i_match <- which(z.[-i]==z.[i])
  i_match_true <- which(id.[-i]==id.[i])
  i_common <- length(intersect(i_match,i_match_true))
  i_false_neg <- length(i_match_true)-i_common
  return(i_false_neg)
}

n_matches_fun <- function(z){
  sizes <- table(z)
  n_match <- sum(vapply(X = sizes,FUN = function(s){return(choose(s,2))},FUN.VALUE = 1))
  return(sum(n_match))
}