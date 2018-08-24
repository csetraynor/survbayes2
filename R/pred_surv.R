#' Posterior prediction of survival outcome
#'
#' This function predicts the survival from the draws of the posterior.
#'
#' @param post posterior draws in matrix format \cr
#' @param strata strata
#' @param obs original dataframe defaults to ic2surv
#' @return d a plot dataset
#' @export
#' @importFrom magrittr %>%
pred_surv <- function(x = NULL, long_x = NULL,  mod, unix = FALSE, patient_id = "patient_id"){
  if(is.null(long_x) ) {
    long_x <- gen_stan_dat(x)
  }
  post <- suppressWarnings( rstanarm::posterior_linpred(object = mod, newdata = long_x) )
  
  if(unix){
    ncores <- parallel::detectCores()
  } else {
    ncores <-  1
  }
  
  if(patient_id %in% colnames(long_x)){
    plot.list <- parallel::mclapply(1:nrow(post), function(i){
      long_x$loghaz <- post[i, ]
      long_x %>%
        dplyr::group_by(patient_id) %>%
        dplyr::mutate(surv = exp( - cumsum(exp(loghaz))) ) %>%
        dplyr::ungroup() %>%
        dplyr::select(surv)
    }, mc.cores = ncores )
  } else {
    plot.list <- parallel::mclapply(1:nrow(post), function(i){
      long_x$loghaz <- post[i, ]
      long_x %>%
        dplyr::mutate(surv = exp( - cumsum(exp(loghaz))) ) %>%
        dplyr::select(surv)
    },  mc.cores = ncores)
  }
  post.matrix <- do.call(cbind, plot.list)
  as.matrix(post.matrix)
}

get_survProb <- function(newdat, time = "time", surv = "surv", patient_id = "patient_id"){
  test <- newdat[, grepl(paste(c(time, surv, patient_id), collapse = "|"), colnames(newdat))] 
                 
  probs <- split(test, as.factor(test[[patient_id]]) ) %>%
    Reduce(function(dtf1,dtf2) full_join(dtf1, dtf2, by= time), .) %>%
    select(contains(surv)) %>%
    as.matrix() 
  probs <- t(probs)
  probs <- cbind( rep(1, nrow(probs) ), probs) 
  colnames(probs) <- 1:ncol(probs)
  probs
  
}
