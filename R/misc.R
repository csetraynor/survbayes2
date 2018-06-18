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
gen_stan_dat <- function(dat, status = "status", time = "time", timepoints = T) {
  # prepare for longdat formating
  dat$sample_id <- 1:nrow(dat)  #create sample id
  # get unique times: only event times equivalent to Cox model
  if(length(timepoints) > 1 ){
    times <- timepoints
  } else{
    times <- dat[dat[[status]], ]
    times <- times[order(unique(unlist(times[, time]))), time]
  }

  form <- as.formula(paste0("Surv(", time, " ,", status, " )", "~."))
  longdat <- survival::survSplit(form, data = dat, cut = times)
  # create time point id
  longdat <- longdat %>% dplyr::group_by(sample_id) %>% dplyr::mutate(t_id = seq(n()))
  # calculate log duration for off-set variable
  longdat$dtime <- longdat[time] - longdat[["tstart"]]
  longdat$log_dtime <- as.double( unlist( log(longdat$dtime) ) )

  longdat

  # stan_data <- list(N = nrow(longdata), S =
  # dplyr::n_distinct(longdata$sample), T = length(times), s =
  # array(as.numeric(longdata$s)), t = array(longdata$t), event =
  # array(longdata$deceased), t_obs = array(longdata$os_months), t_dur =
  # array(longdata$t_dur)) stan_data
}
