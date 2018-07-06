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

post_surv <- function(x, surv_form, prior = NULL, ...) {
  
  form <- as.formula( paste0(c( "status~1+offset(log_dtime)+s(time)+", paste(surv_form, collapse = "+")), collapse = ""))
  
  #prepare long format dataset
  long_x <- gen_stan_dat(x)
  
  if(is.list(prior)){
    m1.stan_gam <- brms::brm(brms::bf(form),
                             data = long_x, family = poisson(), ...)
  } else {
    m1.stan_gam <- rstanarm::stan_gamm4(form, data = long_x ,
                                        family= poisson(),
                                        ... )
  }
  return(m1.stan_gam)
}



#' Prediction of survival from model
#'
#'
#' Fits bayesian semi-parametric survival models
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
#'@export pred_sbm

pred_sbm <- function(x, surv_form, prior = NULL, ...) {
  train_data <- rsample::analysis(x)
  
  model.map2stan <- post_surv(x = train_data , surv_form = surv_form, prior = prior, iter = 12000, thin = 10, adapt_delta = 0.99)

  timepoints <-  seq(0, max(train_data$time),
                     length.out = 100L)
  
  test_data <- rsample::assessment(x)
  newdat <-  gen_new.frame(dat = test_data, timepoints = timepoints)
  newdat <- newdat[ ,match(c(surv_form[!grepl(":", surv_form)], "log_dtime", "time", "sample_id"), colnames(newdat))]
  post <- suppressWarnings(posterior_linpred(model.map2stan, newdata = newdat) )
  pred_frame <- pred_surv(post = post, longdata = newdat)

  pred_frame
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
#'@export  get_bs2

get_bs2 <- function(pred_frame, x ) {

  train_data <- rsample::analysis(x)
  
  timepoints <-  seq(0, max(train_data$time),
                     length.out = 100L)
  
  test_data <- rsample::assessment(x)
  newdat <-  gen_new.frame(dat = test_data, timepoints = timepoints)
  newdat <- newdat[ ,match(c(surv_form[!grepl(":", surv_form)], "log_dtime", "time", "sample_id"), colnames(newdat))]
  
  ibrier <- sapply(seq_along(colnames(pred_frame)), function(i){
    newdat$surv <- as.vector( pred_frame[ ,i] ) 
    
    probs <- get_survProb(newdat = newdat)
    
    brier <- suppressMessages( pec::pec(probs, Surv(time, status) ~ 1,
                      data = test_data,
                      maxtime = max(timepoints),
                      exact = FALSE,
                      exactness = 99L) )
    integrate_tdbrier(brier)
    
  })
  mean(ibrier)
}
 

#' Get concordance index
#'
#'
#' Get concordance index from a survbayes model
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
#'@export  get_ci2

get_ci2 <- function(pred_frame, x ) {
  
  train_data <- rsample::analysis(x)
  
  timepoints <-  seq(0, max(train_data$time),
                     length.out = 100L)
  
  test_data <- rsample::assessment(x)
  newdat <-  gen_new.frame(dat = test_data, timepoints = timepoints)
  newdat <- newdat[ ,match(c(surv_form[!grepl(":", surv_form)], "log_dtime", "time", "sample_id"), colnames(newdat))]
  
  concordance <- sapply(seq_along(colnames(pred_frame)), function(i){
    newdat$surv <- as.vector( pred_frame[ ,i] ) 
    
    probs <- get_survProb(newdat = newdat)
    
    
    tail(suppressWarnings(suppressMessages( (pec::cindex(probs, 
                                           Surv(time, status) ~ 1,
                                        data = test_data,
                                  pred.times = timepoints,
                                  eval.times = timepoints))))$AppCindex$matrix , n=1)
    
  })
} 