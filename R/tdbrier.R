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
