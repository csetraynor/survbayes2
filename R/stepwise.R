#' Cox model step-wise fitting 
#'
#' Fit model by step-wise variable selection.
#'
#' @param
#' train : train data \cr
#' time : time variable \cr
#' status : status variable \cr
#' method: forward, backward or both \cr
#' fit : fit model: "Stepwise", "Random forest"
#' @return a coxph forward Stepwise fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import partykit
#' @import grid
#' @import libcoin
#' @import mvtnorm
#' @import rpart
#' @import doMC

cox_step <- function(x, time = "time", status = "status", surv_form = NULL, penalty = "AIC" ){
  
  #dummy fit
  fitnull <- survival::coxph(Surv(time, status) ~ 1,
                                           data = x) 
  
  
  ##### create formula

  testBeta <- as.formula(paste("~", paste(surv_form, collapse = "+"), sep = " " ) )
  if(penalty == "AIC"){
    k = 2
  }else{
    if(penalty == "BIC"){
      k = log(nrow(trainX))
    }
  }
  #forward Stepwise selection
  mod <- suppressWarnings(MASS::stepAIC(fitnull,
                                        scope = list(upper = testBeta, lower = ~1) ,
                                        trace = TRUE,
                                        direction = "forward",
                                        k = k,
                                        steps =  nrow(x)-1 ))
}


#' Bayes model step-wise fitting 
#'
#' Fit model by step-wise variable selection.
#'
#' @param
#' train : train data \cr
#' time : time variable \cr
#' status : status variable \cr
#' method: forward, backward or both \cr
#' fit : fit model: "Stepwise", "Random forest"
#' @return stan_step fit model
#' @export stan_surv_step
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import grid
#' @import libcoin
#' @import mvtnorm
#' @import rpart
#' @import doMC

stan_surv_step <- function(x, time = "time", status = "status", scope = NULL, 
                      criterion = "loo", verbose = TRUE, fit = NULL, 
                      stan_function = NULL){
  
  var <- NULL
  scope.full <- scope
  if(is.null(fit)){
    #Fit null model
    current.model <- post_surv(x = x) 
  } else {
    current.model <- fit; rm(fit);
  }
  #calculate loo for null model, next.dev.estimate is initial model or null model
  current.loo <- rstanarm::loo(current.model, cores = parallel::detectCores() )
  next.dev.estimate <- -2*current.loo$estimates["elpd_loo","Estimate"] # working in deviance scale
  print("Null model"); print(current.loo ); print(paste0(c("Deviance: ", next.dev.estimate )) )
  current.dev.estimate <- Inf
  i <- 0
  while(current.dev.estimate > next.dev.estimate){
    #remove var chosen in previous step
    if(!is.null(var)){
      current.model <- fit.models[[match(max.loo.estimate, loo.estimates)]]
      current.loo <- loo.models[[match(max.loo.estimate, loo.estimates)]]
      print(paste0(c("Step ", i, " : ", var ) ) ); print(current.loo ); print(paste0(c("Deviance: ", next.dev.estimate )) );
      #update current model
      if(isTRUE(all.equal(var, scope.full) ) ){
        break
      }
      scope <- scope.full[-match(x = var, table = scope.full)]  
    } else {
      scope <- scope.full
    }; 
    
    #Model fitting with rStan
    fit.models <- lapply(scope, function(p){
      if(!is.null(var)){
        par <- c(p, var)
      } else {
        #add var
        par <- p
      };  
      #train model
      post_surv(x = x, surv_form = par)
    })
    #extract loo for each model
    loo.models <- lapply(fit.models, function(m){
      rstanarm::loo(m, cores = detectCores()) 
    })
    loo.estimates <- sapply(loo.models, function(l){
      l$estimates["elpd_loo","Estimate"]
    });max.loo.estimate <- max(unlist(loo.estimates))
    
    #update var
    var <- c(var, scope[match(max.loo.estimate, loo.estimates)] ); 
    #update current dev
    inter.dev.estimate <- next.dev.estimate;
    next.dev.estimate <- -2*max.loo.estimate; 
    current.dev.estimate <- inter.dev.estimate; rm(inter.dev.estimate);
    i <- i + 1;
  }
  out <- list(current.model, current.loo)
  names(out) <- c("model", "loo")
  return(out)
}



for.setattr <- function(static_name, myList) {
  for (i in seq_along(myList)) {
    var_name <- paste(static_name, i, sep = "_")
    setattr(myList[[i]], name = var_name, value = static_name)
  }
  
  ### Just a BACKUP
 
}
