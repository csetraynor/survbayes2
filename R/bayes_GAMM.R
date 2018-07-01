#' Posterior Hazard
#'
#'
#' Bayesian semiparametric survival model with GAMM
#'
#' @param x data
#' @param surv_mod vector with covariates
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
#'@export post_surv

post_surv <- function(x, surv_form, prior = NULL, iter = 12000, ...) {
  
  form <- as.formula( paste0(c( "status~1+offset(log_dtime)+s(time)+", paste(surv_form, collapse = "+")), collapse = ""))
  
  #prepare long format dataset
  long_x <- gen_stan_dat(x)
  
  if(is.list(prior)){
    m1.stan_gam <- brms::brm(brms::bf(form),
                             data = long_x, family = poisson(),  seed = 17,
                             iter = iter, warmup = 1000, thin = 10, 
                             control = list(adapt_delta = 0.99),
                             prior = prior)
  } else {
    m1.stan_gam <- rstanarm::stan_gamm4(form, data = long_x , family= poisson(),
                                        seed = 17,
                                        iter = iter, warmup = 1000, thin = 10, 
                                        control = list(adapt_delta = 0.99))
  }
  return(m1.stan_gam)
}



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
#'@export get_ibrier_bm

get_ibrier_bm <- function(x, surv_form, prior = NULL, ...) {
  x <- rsample::analysis(x)


  m1.stan_gam <- post_surv(x = x , surv_form = surv_form, prior = prior)

  timepoints <-  seq(0, max(x$time),
                     length.out = 100L)
  
  pred_dat <- rsample::assessment(x)
  newdat <-  gen_new.frame(dat = pred_dat, timepoints = timepoints)
  newdat <- newdat[ ,match(c(surv_form, "log_dtime", "time", "sample_id"), colnames(newdat))]
  post <- posterior_linpred(m1.stan_gam, newdata = newdat)
  pred_frame <- pred_surv(post = post, longdata = newdat)
  newdat$surv_mean <- as.vector( apply(pred_frame, 1, mean) ) 
  
  probs <- get_survProb(newdat = newdat)

  brier <- pec::pec(probs, Surv(time, status) ~ 1,
                    data = pred_dat,
                    maxtime = max(timepoints),
                    exact = FALSE,
                    exactness = 99L)
  brier
}
  