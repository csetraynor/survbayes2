#' Get c-index
#'
#'
#' The efficacy of the survival model can be measured by
#'  the concordance statistic
#'
#' @param data For the default functions, a datframe containing survival
#' (time), and status (0:censored/1:event), and the explanatory variables.
#' @param mod Coxph model object fitted with coxph (survival).
#' @return A cindex object
#' @seealso [coxph]
#' @keywords cindex
#' @examples
#' require(survival)
#' require(dplyr)
#' data(lung)
#' lung <- lung %>%
#' mutate(status = (status == 2))
#'
#' mod <- coxph(Surv(time, status)~ age, data = lung)
#'
#' get_cindex(lung, mod)
#'
#' @export get_cindex
#' @author Carlos S Traynor
#' @references
#'
#'  Terry M. Therneau and Patricia M. Grambsch (2000).
#'   _Modeling Survival Data: Extending the Cox Model_.
#'   Springer, New York. ISBN 0-387-98784-3.
#'   @export get_cindex
#'
get_cindex <- function(data, mod,...)
  UseMethod("get_cindex")

#' @export
#' @rdname get_cindex
get_cindex <-
  function(test, mod, time = "time", status = "status", Unix = FALSE){
    
    timepoints <-   spc::quadrature.nodes.weights(100, type="GL", x1= 0 , x2= max(test[[time]][test[[status]] ]) )$nodes
    
    if(any( class(mod)  == "coxph") ){
      probs <- pec::predictSurvProb(mod, newdata = test, times = timepoints)
      
      statistic <- suppressMessages( pec::cindex(probs, 
                                    Surv(time, status) ~ 1,
                                    data = test,
                                    pred.times = timepoints,
                                    eval.times = timepoints ) )
      
      apperror <- unlist(statistic$AppCindex)
      time <- timepoints
      concordance <-  tail(apperror, n = 1)
      ref <- 0.5
      out <- list(concordance, apperror, ref, time)
      names(out) <- c( "concordance", "apperror", "reference","time")
      
      out
      
    } else {
      if(any( class(mod) == "stanreg" ) ){
        
        newdat <-  gen_new_frame(dat = test,
                                 timepoints = timepoints)
        pred_frame <- pred_surv(long_x = newdat, 
                                mod = mod,
                                unix = Unix)
        
        models.cindex <- lapply(seq_along(colnames(pred_frame)), function(i){
          newdat$surv <- as.vector( pred_frame[ ,i] ) 
          
          probs <- get_survProb(newdat = newdat)
          
          c.index <-       suppressMessages( pec::cindex(probs, 
                                                       Surv(time, status) ~ 1,
                                                       data = test,
                                                       pred.times = timepoints,
                                                       eval.times = timepoints ) )
          c.index
        })
        models.cindex.app <- lapply(models.brier, function(b){
          b$AppCindex
        })
        models.concordance <- sapply(models.cindex.app, function(b){
          tail(unlist(b), n = 1)
        })
        models.cindex.app <- do.call(cbind, models.cindex.app )
        apperror <- apply(models.cindex.app, 1, mean)
        time <- timepoints
        concordance <- mean(models.concordance)
        ref <- 0.5
        out <- list(concordance, apperror, ref, time)
        names(out) <- c( "concordance", "apperror", "reference","time")
        
        out
        
      } else {
        "Print object must be either coxph or stanreg"
      }
    }
  }
