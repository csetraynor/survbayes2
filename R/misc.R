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

#' Generate data for stan
#'
#' This function generates a long format dataframe for
#' repeated tte, time dependent covariates.
#'
#' @param d a dataset \cr
#' @return d a long dataset
#' @export
#' @importFrom rlang !!
#' @importFrom magrittr %>%
gen_stan_surv <- function(dat, status = "status", time = "time", timepoints = T) {
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
}
#' String method for strata in surv analysis
#'
#' This function is a string method to facilitate ploting.
#'
#' @param c a character frame \cr
#' @return d a clean dataset
#' @export
#' @importFrom magrittr %>%
str_strata <- function(c){
  test <- strsplit(as.character( c$strata ), "\\,")
  test <-as_data_frame( do.call(rbind, test) )
  colnames(test) <- gsub("=.*| ","",  test[1, ] ) 
  strata_replace <- function(x){
    x <- gsub(".*=| ", "", x) 
    x
  }
  test <- test %>%
    dplyr::mutate_all(funs(strata_replace))
  test
}

#' Expand dataset accordin to calendar time
#'
#' @param dat a data frame \cr
#' @return d a long dataset
#' @export
#' @importFrom magrittr %>%
gen_new.frame <- function(dat, time = "time", id = "patient_id", timepoints){

  
  data_cal <- lapply(seq_along( dat[[id]] ), function(x) {
    id_sq <- rep(x, length(timepoints) - 1)
    dat.list <- data.frame(t_start = timepoints[-length(timepoints)],
                           time = timepoints[-1],
                      sample_id = id_sq) 
    dat.list$dtime <- dat.list$time - dat.list$t_start;
    dat.list$log_dtime <- log(dat.list$dtime)
      
    test_data <- dat[ ,!grepl(time, colnames(dat))]
    merge(dat.list, test_data[x, ], sort = FALSE)
  }
  )
  do.call(rbind, data_cal)
}

get_survival_function  <- function(loghaz) {
  exp( - cumsum(exp(loghaz))) 
}

my_replace <- function(x){
  
  x <- gsub("1-Sep", "Sep_1", x)
  x <- gsub("\\-", "_", x)
  gsub("`", "", x)
}
