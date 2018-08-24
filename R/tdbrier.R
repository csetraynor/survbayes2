#'  Calculate tdBrier
#'
#' These functions calculate the survival analysis metric measured of a
#' system compared to a hold-out test set. The measurement and the "truth"
#' have a survival time and a censoring indicator 0/1 indicating if the event
#' result or the event.
#'
#'
#' The Brier score is defined as the squared distance between the
#' expected survival probability and the observed survival.
#' Therefore, it measures the discrepancy between observation
#' and model-based prediction.
#'
#' The integrated Brier Score summarises the Brier Score over the range
#' of observed events.Similar to the original Brier score [40] the iBrier:
#' ranges from 0 to 1; the model with an out-of-training sample value closer
#' to 0 outperforms the rest.
#' @aliases tdbrier tdbrier.model.list tdbrier.int.matrix
#' stdbrier.int.reference
#' @param data For the default functions, a datframe containing survival
#' (time), and status (0:censored/1:event), and the explanatory variables.
#' @param mod Coxph model object fitted with coxph (survival).
#' @return A tdBrier object
#' @seealso [iBrier]
#' @keywords brier
#' @examples
#' require(survival)
#' require(dplyr)
#' data(lung)
#' lung <- lung %>%
#' mutate(status = (status == 2))
#'
#' mod <- coxph(Surv(time, status)~ age, data = lung)
#'
#' tdbrier <- get_tdbrier(lung, mod)
#' integrate_tdbrier(tdroc)
#'
#' @export tdbrier
#' @author Carlos S Traynor
#' @references
#'
#'  Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012).
#' Evaluating Random Forests for Survival Analysis Using Prediction Error
#' Curves. Journal of Statistical Software, 50(11), 1-23.
#' URL http://www.jstatsoft.org/v50/i11/.
#' @export tdbrier
tdbrier <- function(data, mod,...)
  UseMethod("tdbrier")

#' @export
#' @rdname tdbrier

get_tdbrier <-
  function(data, mod,  ...) {
    x <- rsample::analysis(data)

    features <- names(mod$coefficients)


    pred_dat <- rsample::assessment(data)


    #Create grid of equidistant time points for testing
    timepoints <-  seq(0, max(x$time),
                       length.out = 100L)
    probs <- pec::predictSurvProb(mod,
                                  newdata = pred_dat,
                                  times = timepoints)
    #Calculate brier score
    suppressMessages(brier <- pec::pec(probs, Surv(time, status) ~ 1,
                                       data = pred_dat,
                                       maxtime = max(timepoints),
                                       exact = FALSE,
                                       exactness = 99L))

    return(brier)
  }

#' @rdname tdbrier
#' @export
integrate_tdbrier <-
  function(x, ...) {
    stop <- max(x$time[!is.na(x$AppErr$matrix)])
    ibrier <- pec::crps(x, models = "matrix", times = stop)[1]
    ibrier <- unlist(ibrier)
    return(ibrier)
  }
#' @export
#' @rdname tdbrier
integrate_tdbrier_reference <-
  function(x, ...) {
    stop <- max(x$time[!is.na(x$AppErr$Reference)])
    ibrier <- pec::crps(x, models = "Reference", times = stop)[1]
    ibrier <- unlist(ibrier)
    return(ibrier)
  }

#'  Calculate survival Brier score
#'
#' These functions calculate the survival analysis metric measured of a
#' system compared to a hold-out test set. The measurement and the "truth"
#' have a survival time and a censoring indicator 0/1 indicating if the event
#' result or the event.
#'
#'
#' The Brier score is defined as the squared distance between the
#' expected survival probability and the observed survival.
#' Therefore, it measures the discrepancy between observation
#' and model-based prediction.
#'
#' The integrated Brier Score summarises the Brier Score over the range
#' of observed events.Similar to the original Brier score [40] the iBrier:
#' ranges from 0 to 1; the model with an out-of-training sample value closer
#' to 0 outperforms the rest.
#' @aliases tdbrier tdbrier.model.list tdbrier.int.matrix
#' stdbrier.int.reference
#' @param data For the default functions, a datframe containing survival
#' (time), and status (0:censored/1:event), and the explanatory variables.
#' @param mod Coxph model object fitted with coxph (survival).
#' @return A tdBrier object
#' @seealso [iBrier]
#' @keywords brier
#' @examples
#' require(survival)
#' require(dplyr)
#' data(lung)
#' lung <- lung %>%
#' mutate(status = (status == 2))
#'
#' mod <- coxph(Surv(time, status)~ age, data = lung)
#'
#' tdbrier <- get_tdbrier(lung, mod)
#' integrate_tdbrier(tdroc)
#'
#' @export get_survbrier
#' @author Carlos S Traynor
#' @references
#'
#'  Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012).
#' Evaluating Random Forests for Survival Analysis Using Prediction Error
#' Curves. Journal of Statistical Software, 50(11), 1-23.
#' URL http://www.jstatsoft.org/v50/i11/.
get_survbrier <- function(test, mod, time = "time", status = "status", Unix = FALSE){
  
  timepoints <-   spc::quadrature.nodes.weights(100, type="GL", x1= 0 , x2= max(test[[time]][test[[status]] ] ) )$nodes
  
  if(any( class(mod)  == "coxph") ){
    probs <- pec::predictSurvProb(mod, newdata = test, times = timepoints)
    
   cox.brier <- suppressMessages( pec::pec(probs, Surv(time, status) ~ 1,
                                        data = test,
                                        maxtime = max(timepoints),
                                        exact = FALSE,
                                        exactness = 99L) )
   apperror <- cox.brier$AppErr$matrix
   time <- timepoints
   ref <- cox.brier$AppErr$Reference
   ibrier <- integrate_tdbrier(cox.brier)
   ibrier_ref <- integrate_tdbrier_reference(cox.brier)
   rsquaredbs <- 1-(ibrier/ibrier_ref)
   out <- list(rsquaredbs, ibrier, apperror, ref, time)
   names(out) <- c( "rsquaredbs","ibrier", "apperror", "reference","time")
   out
   
  } else {
    if(any( class(mod) == "stanreg" ) ){
      
      newdat <-  gen_new_frame(dat = test,
                               timepoints = timepoints)
      pred_frame <- pred_surv(long_x = newdat, 
                              mod = mod,
                              unix = Unix)
      
      models.brier <- lapply(seq_along(colnames(pred_frame)), function(i){
        newdat$surv <- as.vector( pred_frame[ ,i] ) 
        
        probs <- get_survProb(newdat = newdat)
        
        brier <- suppressMessages( pec::pec(probs, Surv(time, status) ~ 1,
                                            data = test,
                                            maxtime = max(timepoints),
                                            exact = FALSE,
                                            exactness = 99L) )
        brier
      })
      models.brier.error <- lapply(models.brier, function(b){
        b$AppErr$matrix
      })
      models.ibrier <- sapply(models.brier, function(b){
        integrate_tdbrier(b)
      })
      models.brier.error <- do.call(cbind, models.brier.error )
      apperror <- apply(models.brier.error, 1, mean)
      time <- timepoints
      ibrier <- mean(models.ibrier)
      ref <- models.brier[[1]]$AppErr$Reference
      ibrier_ref <- integrate_tdbrier_reference(models.brier[[1]])
      rsquaredbs <- 1-(ibrier/ibrier_ref)
      out <- list(rsquaredbs, ibrier, apperror, ref, time)
      names(out) <- c( "rsquaredbs","ibrier", "apperror", "reference","time")
      
      out
      
    } else {
      "Print object must be either coxph or stanreg"
    }
  }
}

#' Get brier score
#'
#'
#' Get brier score from a survbayes model
#'
#' @param x data
#' @param mod Coxph model object fitted with coxph (survival).
#' @return mod
#' @seealso [coxph]
#' @keywords coxph
#'
#' @author Carlos S Traynor
#' @references
#'
#'  Terry M. Therneau and Patricia M. Grambsch (2000).
#'   Modeling Survival Data: Extending the Cox Model.
#'   Springer, New York. ISBN 0-387-98784-3.
#'@importFrom prodlim Hist
#'@export  get_bs2

get_bs2 <- function(pred_frame, x, surv_form , preProcValues = NULL) {
  
  train_data <- rsample::analysis(x)
  
  timepoints <-  seq(0, max(train_data$time),
                     length.out = 100L)
  
  test_data <- rsample::assessment(x)
  #Standardise test data 
  if(is.null(preProcValues)){
    preProcValues <- caret::preProcess(train_data[c("age_at_diagnosis", "npi")],
                                       method = c("center", "scale") )
  }
  
  testTransformed <- predict(preProcValues, test_data[c("age_at_diagnosis", "npi")])
  colnames(testTransformed) = c("age_std", "npi_std")
  
  test_data <- cbind(test_data, testTransformed)
  
  newdat <-  gen_new.frame(dat = test_data, timepoints = timepoints)
  newdat <- newdat[ ,match(c(surv_form[!grepl(":", surv_form)], "log_dtime", "time", "sample_id"), colnames(newdat))]
  
  models.brier <- lapply(seq_along(colnames(pred_frame)), function(i){
    newdat$surv <- as.vector( pred_frame[ ,i] ) 
    
    probs <- get_survProb(newdat = newdat)
    
    brier <- suppressMessages( pec::pec(probs, Surv(time, status) ~ 1,
                                        data = test_data,
                                        maxtime = max(timepoints),
                                        exact = FALSE,
                                        exactness = 99L) )
    brier
  })
  models.brier.error <- lapply(models.brier, function(b){
    b$AppErr$matrix
  })
  models.ibrier <- sapply(models.brier, function(b){
    integrate_tdbrier(b)
  })
  models.brier.error <- do.call(cbind, models.brier.error )
  apperror <- apply(x = models.brier.error, MARGIN = 1, mean)
  time <- timepoints
  ibrier <- mean(models.ibrier)
  ref <- models.brier[[1]]$AppErr$Reference
  out <- list(ibrier, apperror, ref, time)
  names(out) <- c( "ibrier", "apperror", "reference","time")
  
}
