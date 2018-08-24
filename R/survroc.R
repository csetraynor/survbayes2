#' Calculate tdROC
#'
#'
#' The time dependent ROC analysis (tdROC) allows
#' the user to visualise the sensitivity (true positive rate) and
#' specificity (false positive rate) for all possible cut-offs for the
#' predicted survival.
#' #'
#' @param data For the default functions, a datframe containing survival
#' (time), and status (0:censored/1:event), and the explanatory variables.
#' @param mod Coxph model object fitted with coxph (survival).
#' @return A tdROC object
#' @seealso [iROC]
#' @keywords tdroc
#' @examples
#' require(survival)
#' require(dplyr)
#' data(lung)
#' lung <- lung %>%
#' mutate(status = (status == 2))
#'
#' mod <- coxph(Surv(time, status)~ age, data = lung)
#'
#' tdroc <- get_tdroc(lung, mod)
#' integrate_tdroc(tdroc)
#'
#' @export tdroc
#' @author Carlos S Traynor
#' @references
#'  Liang Li, Cai Wu Department of Biostatistics and
#'  The University of Texas MD Anderson Cancer Center (2016).
#'  tdROC: Nonparametric Estimation of Time-Dependent ROC Curve
#'   from Right Censored Survival Data. R package version 1.0.
#'  https://CRAN.R-project.org/package=tdROC
#' @export tdroc
tdroc <- function(data, mod,...)
  UseMethod("tdroc")

#' @export
#' @rdname tdroc
get_tdroc <-
  function(data, mod, ...) {
    pred_dat <- rsample::assessment(data)
    probs <- predict(mod, newdata = pred_dat, type = "lp")
    
    roc <- tdROC::tdROC(X = probs[!is.na(probs)],
                        Y = pred_dat$time[!is.na(probs)],
                        delta = pred_dat$status[!is.na(probs)],
                        tau = max(pred_dat$time),
                        n.grid = 1000)
    return(roc)
  }
#' @export
#' @rdname tdroc
integrate_tdroc  <-
  function( mod, ...) {
    mod$AUC[1] %>% unlist
  }


#'  Calculate survival ROC
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
#' @export get_survroc
#' @author Carlos S Traynor
#' @references
#'
#'  Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012).
#' Evaluating Random Forests for Survival Analysis Using Prediction Error
#' Curves. Journal of Statistical Software, 50(11), 1-23.
#' URL http://www.jstatsoft.org/v50/i11/.
get_survroc <- function(test, mod, time = "time", status = "status", Unix = FALSE){
  
  timepoints <-   spc::quadrature.nodes.weights(100, type="GL", x1= 0 , x2= max(test[[time]][test[[status]] ] ) )$nodes
  
  if(any( class(mod)  == "coxph") ){
    probs <- predict(mod, newdata = test, type = "lp")
    
    statistic <- tdROC::tdROC(X = probs[!is.na(probs)],
                        Y = test[[time]][!is.na(probs)],
                        delta = test[[status]][!is.na(probs)],
                        tau = max(test[[time]][test[[status]] ]),
                        n.grid = 1000)
    
    iroc <- integrate_tdroc(statistic)

    sens <-  statistic$ROC$sens
    
    spec <-  statistic$ROC$spec
    
    grid <- statistic[[1]]$ROC$grid
    
    out <- list(iroc, sens, spec, grid)
    names(out) <- c( "iroc", "sens", "spec","grid")
    
    out
    
    
    
  } else {
    if(any( class(mod) == "stanreg" ) ){
      
      mod.samples <- as.data.frame(mod)
      vars <- attr( mod$terms, "term.labels")
      vars <- vars[-match("time", vars)]
      
      x <- test_luscs[vars]
      
      X <- data.frame(id = seq_along(1:nrow(x)))
      
      for(i in seq_along(colnames(x)) ){
        
        if(is.numeric(x[[i]])){
          matrixna <- x[i]
          X <- cbind(X, matrixna)
        } else{
          x[[i]] <- as.factor(x[[i]])
        coluna <- x[i]
        matrixna <- model.matrix(~ ., coluna)[ ,-1]
        matrixna <- as.data.frame(matrixna)
        if('matrixna' %in% colnames(matrixna) ){
          colnames(matrixna) <- paste0(colnames(coluna), "1")
        }
        X <- cbind(X, matrixna)
        }
      }
      X <- X[-match("id", colnames(X) )]
      mod.samples <- mod.samples[colnames(X)]
      
      
      roc.out <- lapply(seq_along(1:nrow(mod.samples)), function(samp){
        probs <- as.matrix(X) %*% as.vector(unlist( mod.samples[samp, ] ))
        tdROC::tdROC(X = probs[!is.na(probs)],
                     Y = test[[time]][!is.na(probs)],
                     delta = test[[status]][!is.na(probs)],
                     tau = max(test[[time]][test[[status]] ]),
                     n.grid = 1000)
      } )
      
      iroc <- sapply(roc.out, function(statistic){
        integrate_tdroc(statistic)
      })
      
      iroc <- mean(unlist( iroc))
      
      sens <- lapply(roc.out, function(sensibility){
        sensibility$ROC$sens
      })
      
      sens <- do.call(cbind, sens)
      sens <- apply(sens, 1, mean)
      
      spec <- lapply(roc.out, function(specificity){
        sensibility$ROC$spec
      })
      
      spec <- do.call(cbind, spec)
      spec <- apply(spec, 1, mean)
      
      grid <- roc.out[[1]]$ROC$grid
      
      out <- list(iroc, sens, spec, grid)
      names(out) <- c( "iroc", "sens", "spec","grid")
      
      out
      
    } else {
      "Print object must be either coxph or stanreg"
    }
  }
}
