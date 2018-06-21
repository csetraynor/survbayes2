#' Model fitting
#'
#'
#' Fits coxph model.
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
#'@export mod_coxfit

mod_coxfit <- function(x, surv_form, iter = 20, inits = NA_character_, ...) {
  x <- rsample::analysis(x)

  if(is.character(inits)){
    form <- as.formula( paste0(c( "Surv(time, status)", surv_form), collapse = ""))
    coxph(form , data = x)
  } else {
      coef <- paste(surv_form, collapse = "+")
      form <- as.formula(paste(c( "Surv(time, status)~", coef), collapse = "") )
      coxph(form , data = x, init = inits, control = coxph.control(iter.max = iter) ) 
    }
}
