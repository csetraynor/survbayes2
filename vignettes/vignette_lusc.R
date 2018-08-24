# !diagnostics off
#load libraries
library(readr)
library(dplyr)
library(tidyprojectAZ)
library(ggplot2)
library(survival)
theme_set(theme_bw())
devtools::document()

#load data
data("luscc")

luscc <- luscc[complete.cases(luscc), ] %>%
  dplyr::mutate(time = os_months,
                status = os_status) 
luscc$time[luscc$time == 0] <- 1e-3

glimpse(luscc); summary(luscc)

set.seed(10)
lusc_data <- resample_dataset(luscc, strata = "status", sizeto = 100)

cox.model <- coxph(Surv(time = time , event = status) ~ ipi , data = lusc_data)

gam.model <- post_surv( lusc_data,
                        surv_form = c("ipi"))

#Assesing the predictive performance : the survival ROC curve

timepoints <-  seq(0, max(lusc_data$time),
                   length.out = 100L)

newdat <-  gen_new.frame(dat = lusc_data, timepoints = timepoints)
surv_form <- c("ipi")
newdat <- newdat[ ,match(c(surv_form[!grepl(":", surv_form)], "log_dtime", "time", "sample_id"), colnames(newdat))]

post <- suppressWarnings(posterior_linpred(gam.model, newdata = newdat) )
plot.list <- lapply(1:nrow(post), function(i){
  newdat$loghaz <- post[i, ]
  newdat %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(surv = exp( - cumsum(exp(loghaz))) ) %>%
    dplyr::ungroup() %>%
    dplyr::select(surv)
})
post.matrix <- do.call(cbind, plot.list)
pred_frame <- as.matrix(post.matrix)

#Assesing the predictive performance: survival Brier score

ibrier <- sapply(seq_along(colnames(pred_frame)), function(i){
  newdat$surv <- as.vector( pred_frame[ ,i] ) 
  
  probs <- get_survProb(newdat = newdat)
  
  brier <- suppressMessages( pec::pec(probs, Surv(time, status) ~ 1,
                                      data = lusc_data,
                                      maxtime = max(timepoints),
                                      exact = FALSE,
                                      exactness = 99L) )
  integrate_tdbrier(brier)
})
mean(ibrier)

probs <- predictSurvProb(cox.model, newdata = lusc_data, times = timepoints)

brier <- suppressMessages( pec::pec(probs, Surv(time, status) ~ 1,
                                    data = lusc_data,
                                    maxtime = max(timepoints),
                                    exact = FALSE,
                                    exactness = 99L) )
integrate_tdbrier(brier)



# Test split validation

library(rsample)

set.seed(10)
vfold <- rsample::mc_cv(lusc_data, times = 1, strata = "status")$splits$`1`
training_data <- rsample::analysis(vfold)
test_data <- rsample::assessment(vfold)

glimpse(training_data)

cox.model <- coxph(Surv(time , event = status) ~ ipi , data = training_data)

gam.model <- post_surv( training_data,
                        surv_form = c("ipi") )

timepoints <-  seq(0, max(training_data$time),
                   length.out = 100L)

newdat <-  gen_new.frame(dat = test_data, timepoints = timepoints)
surv_form <- c("ipi")
newdat <- newdat[ ,match(c(surv_form[!grepl(":", surv_form)], "log_dtime", "time", "sample_id"), colnames(newdat))]

post <- suppressWarnings(posterior_linpred(gam.model, newdata = newdat) )
plot.list <- lapply(1:nrow(post), function(i){
  newdat$loghaz <- post[i, ]
  newdat %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(surv = exp( - cumsum(exp(loghaz))) ) %>%
    dplyr::ungroup() %>%
    dplyr::select(surv)
})
post.matrix <- do.call(cbind, plot.list)
pred_frame <- as.matrix(post.matrix)

#Assesing the predictive performance: survival Brier score

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

probs <- pec::predictSurvProb(cox.model, newdata = test_data, times = timepoints)

brier <- suppressMessages( pec::pec(probs, Surv(time, status) ~ 1,
                                    data = test_data,
                                    maxtime = max(timepoints),
                                    exact = FALSE,
                                    exactness = 99L) )
integrate_tdbrier(brier)
