#' Kaplan and Meier estimates
#'
#' This function plots the TKaplan and Meier estimates
#'
#' @param
#' d a dataset \cr
#' x customise plot title \cr
#' event_type codification of event type , default DECEASED
#' @return d a clean dataset
#' @export
#' @importFrom rlang !!
#' @import ggfortify
#' @importFrom magrittr %>%
plot_km <- function(dat, time = time, status = status, event_type = 1, strata = "1" ){
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)
  
  form = as.formula(paste0("Surv( time , status) ~ ", strata))
  mle.surv <- survival::survfit(Surv( time , status) ~ 1,
                                data = dat  )
  obs.mortality <- data.frame(time = mle.surv$time,
                              surv = mle.surv$surv)
  ggplot2::ggplot(obs.mortality, aes(time, surv)) +
    geom_step()
}
