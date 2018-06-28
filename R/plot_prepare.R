#' Prepare plot frame
#'
#' This function prepares plot frame
#'
#' @param post posterior draws in matrix format \cr
#' @param strata strata
#' @param obs original dataframe defaults to ic2surv
#' @return d a plot dataset
#' @export
#' @importFrom magrittr %>%
get_plot.frame <- function(post, strata, obs = ic2surv){
  plot.frame <- tibble::tibble(
    postmean = apply(plot.matrix, 1, mean),
    lower = apply(plot.matrix, 1, quantile, probs = 0.055),
    upper = apply(plot.matrix, 1, quantile, probs = 0.945),
    time = obs %>% arrange(time) %>% select(time) %>% 
      unlist %>% as.double
  ) %>%
    cbind ( obs %>% 
              dplyr::arrange(time) %>% 
              dplyr::select(strata) )
  
  zeros.frame <- tibble::tibble(time=0, postmean=1, 
                            lower = 1, upper = 1) %>%
    cbind(plot.frame %>% select(strata) %>% unique)
  
  plot.frame <- rbind(plot.frame, zeros.frame)
  plot.frame
}

get_km.frame <- function(obs, strata ){
  form <- as.formula(paste0("Surv( time , status) ~ ", paste(strata, collapse = "+")))
  mle.surv <- survfit(form, data = obs  )
  obs.mortality <- data.frame(time = mle.surv$time,
                              surv = mle.surv$surv,
                              strata = summary(mle.surv, censored=T)$strata)
  
  zeros <- data.frame(time=0, surv=1, strata=unique(obs.mortality$strata))
  obs.mortality <- rbind(obs.mortality, zeros)
  
  strata <- str_strata(obs.mortality)
  obs.mortality$strata <- NULL
  obs.mortality <- cbind( obs.mortality, strata)
  obs.mortality 
}

#' Prepare plot with new dataframe
#'
#' This function prepares plot frame
#'
#' @param post posterior draws in matrix format \cr
#' @param strata strata
#' @param obs original dataframe defaults to ic2surv
#' @return d a plot dataset
#' @export
#' @importFrom magrittr %>%
get_plot_new_frame <- function(mod,newdat, timepoints, treatment){
  post <- suppressWarnings(rstanarm::posterior_linpred(m1.map2stan, newdata = newdat))
  
  survdata <- apply(post, MARGIN = 1, get_survival_function)
  
  mean.surv <- apply(survdata, MARGIN = 1, mean)
  mean.upper <- apply(survdata, 1, quantile, 0.945)
  mean.lower <- apply(survdata, 1 , quantile, 0.055)
  
  bc.plot <- data.frame(
    survmean = c(1, mean.surv),
    survlower = c(1, mean.lower),
    survupper = c(1, mean.upper),
    time = c(0, newdat$timepoints)
  ) %>%
    dplyr::mutate(treatment = treatment)
  bc.plot
}

