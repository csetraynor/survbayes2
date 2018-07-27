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

#' Prediction of Survival probability from a Bayesian model
#'
#'
#' Predicts the survival probability for each drawn of the posterior distribution.
#'
#' @param x data
#' @param bgam Bayesian Generalized Additive Models
#' @param preProcess Pre-processing information from training set to Test set.
#' @param ... Ignored
#' @return mod
#' @seealso [gamair]
#' @keywords gamair
#'
#' @author Carlos S Traynor
#' @references
#'
#'  Wood, S.N. (2006) Generalized Additive Models: An Introduction with R 
#'  Chapman & Hall/CRC, Boca Raton, Florida. ISBN 1-58488-474-6
#'@export surv_pred_bgam

surv_pred_bgam <- function(x, bgam, preProcValues = NULL, ...) {
  train_data <- rsample::analysis(x)
  
  #get timepoints
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
  
  #get new dataset
  newdat <-  gen_new.frame(dat = test_data, timepoints = timepoints)
  
  # get surv formula
  form <-  names(bgam$coefficients)
  toMatch <- c(":", "Intercept", "s\\(time")
  toMatch <- form[ !grepl(paste(toMatch,collapse="|"), 
                               form) ]
  surv_form <- gsub("TRUE", "", toMatch)
  
  newdat <- newdat[ ,match(c(surv_form, "log_dtime", "time", "sample_id"), colnames(newdat))]
  #predict linear predictor
  post <- suppressWarnings(posterior_linpred(bgam, newdata = newdat) )
  pred_frame <- pred_surv(post = post, longdata = newdat)
  
  pred_frame
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