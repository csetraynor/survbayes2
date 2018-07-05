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
pred_surv <- function(post, longdata){
  plot.list <- lapply(1:nrow(post), function(i){
    longdata$loghaz <- post[i, ]
    longdata %>%
      dplyr::group_by(sample_id) %>%
      dplyr::mutate(surv = exp( - cumsum(exp(loghaz))) ) %>%
      dplyr::ungroup() %>%
      dplyr::select(surv)
  })
  post.matrix <- do.call(cbind, plot.list)
  as.matrix(post.matrix)
}

get_survProb <- function(newdat, time = "time", surv = "surv", patient_id = "sample_id"){
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
