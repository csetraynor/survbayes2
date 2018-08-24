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
#' @export gen_stan_dat
#' @importFrom rlang !!
#' @importFrom magrittr %>%
gen_stan_dat <- function(dat, status = "status", time = "time", timepoints = T) {
  # prepare for longdat formating
  dat$numeric_id <- 1:nrow(dat)  #create sample id
  # get unique times: only event times equivalent to Cox model
  if(length(timepoints) > 1 ){
    nodes <- timepoints
  } else{
    timepoints <- dat[dat[[status]], ]
    times <- sort( unique( timepoints[[time]] ) )
     nodes <- spc::quadrature.nodes.weights(10, type="GL", x1=min(times), x2=max(times) )$nodes
    # print(nodes)
    #  nodes <- ceiling(seq(1, length(times), length.out = 10))
    # # 
    # timepoints <- times[nodes]
    # nodes <- geomSeries(2, max(times))
  }
  form <- as.formula(paste0("Surv(", time, " ,", status, " )", "~."))
  longdat <- survival::survSplit(form, data = dat, cut = nodes)
  # create time point id
  longdat <- longdat %>% dplyr::group_by(numeric_id) %>% dplyr::mutate(t_id = seq(n()))
  # calculate log duration for off-set variable
  longdat$dtime <- longdat[[time]] - longdat[["tstart"]]
  longdat$log_dtime <- as.double( unlist( log(longdat$dtime) ) )
  longdat
}

geomSeries <- function(base, max) {
  base^(0:floor(log(max, base)))
}



#' Generate long data format for stan
#'
#' This function generates a long format dataframe for
#' repeated tte, time dependent covariates.
#'
#' @param d a dataset \cr
#' @return d a long dataset
#' @export gen_long_dat
#' @importFrom rlang !!
#' @importFrom magrittr %>%
gen_long_dat <- function(dat, status = "status", time = "time") {
  # prepare for longdat formating
  dat$numeric_id <- 1:nrow(dat)  #create sample id
  
  long_time <- lapply(seq_along(dat$numeric_id), function(i){
    obs_time = dat[i, time]
    obs_status = dat[i, status]
    dummy_time = spc::quadrature.nodes.weights(10, type="GL", x1= 0 , x2= obs_time )$nodes
    dummy_status = rep(0, 10)
    data.frame(numeric_id = i,
               time = c(dummy_time, obs_time),
               status = c(dummy_status, obs_status),
               dtime = diff(c(0,dummy_time, obs_time )))
  } )
  
  long_dat <- do.call(rbind.data.frame, long_time)
  long_dat$log_dtime <- log(long_dat$dtime)
  
  dat[time] <- NULL
  dat[status] <- NULL
  
  merge(dat, long_dat, by = "numeric_id")
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
gen_new_frame <- function(dat, time = "time", id = "patient_id", timepoints){

  
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

merge.formula <- function(form1, form2, ...){
  
  # get character strings of the names for the responses 
  # (i.e. left hand sides, lhs)
  lhs1 <- deparse(form1[[2]])
  lhs2 <- deparse(form2[[2]])
  if(lhs1 != lhs2) stop('both formulas must have the same response')
  
  # get character strings of the right hand sides
  rhs1 <- strsplit(deparse(form1[[3]]), " \\+ ")[[1]]
  rhs2 <- strsplit(deparse(form2[[3]]), " \\+ ")[[1]]
  
  # create the merged rhs and lhs in character string form
  rhs <- c(rhs1, rhs2)
  lhs <- lhs1
  
  # put the two sides together with the amazing 
  # reformulate function
  out <- reformulate(rhs, lhs)
  
  # set the environment of the formula (i.e. where should
  # R look for variables when data aren't specified?)
  environment(out) <- parent.frame()
  
  return(out)
}


# this is how you get the addition operator working for formulas
Ops.formula <- function(e1, e2){
  FUN <- .Generic
  if(FUN == '+'){
    out <- merge(e1, e2)
    environment(out) <- parent.frame()
    return(out)
  }
  else stop('can not yet subtract formula objects')
}

resample <- function(x, ...) x[sample.int(length(x), ...)]

resample_stratified <- function(x, sizeto = NULL, strata = NULL, replace_logic = F){
  if(is.null(strata) | is.null(sizeto)){
    stop(print( "Error strata and sizeto should be specified") )
  }else{
    splits <- lapply(seq_along(strata), function(z){
      sizeto <- floor( sizeto / length(unique(x[[strata[z] ]])) / length(strata) ) 
      splits_k <- lapply(unlist( unique(x[[strata[z] ]]) ), function(i){
        probabilities <-  rep(0, nrow(x))
        probabilities[ grep(i, x[[strata[z] ]])] <- 1 
        x[resample(x = 1:nrow(x), size = sizeto, replace = replace_logic,
                   prob = probabilities), ]
      })
      do.call(rbind, splits_k)
    })
    do.call(rbind, splits)
  }
}


convert_blank_to_na <- function(x) {
  if(!purrr::is_character(x)){
    return(x)
    print("Error not character")
  } else {
    ifelse(x == "" | x == "[Not Available]" | x == "--" | x == "not reported" | x == "[Not Applicable]", NA, x)
  }
}