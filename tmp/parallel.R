library(brms)
library(methods)
library(rstan)
library(survival)
library(magrittr)
library(plyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("gen_stan_data.R")


#Load data --------------
mc_samp <- readRDS("mc_samp.RDS")

# # Cox clinical model -------------
# surv_form <-  '~ age_std + npi_std + mastectomy + hormone_therapy + chemotherapy + radio_therapy + her2pos + erpos'
# cox_test <- lapply(mc_samp$splits, function(x){
#   x <- as.data.frame(x$data)[x$in_id,]
#   form <- as.formula( paste0(c( "Surv(time, status)", surv_form), collapse = ""))
#   coxph(form , data = x)
# })
# mc_samp$cox_clinical <- cox_test
# rm(cox_test)

#Bayesian clinical model ----------
surv_form <- c( "age_std", "npi_std", "mastectomy", "hormone_therapy", "chemotherapy", "radio_therapy", "her2pos", "erpos")


# create brms object -----
x <- mc_samp$splits$`1`
x <- as.data.frame(x$data)[x$in_id,]
long_x <- gen_stan_dat(x)

form <- paste0(c( "status~1+offset(log_dtime)+s(time)+", paste(surv_form, collapse = "+")), collapse = "")
bayes_model <-  brm(bf(form),
                         data = long_x, family = poisson() ) 
rm(x);rm(long_x);

# loop object through splits
bayes_test <- lapply(mc_samp$splits, function(x){
  x <- as.data.frame(x$data)[x$in_id,]
  long_x <- gen_stan_dat(x)
  brms::update(bayes_model, newdata = x)
})
mc_samp$bayes_clinical <- bayes_test
rm(bayes_test)

saveRDS(mc_samp, "mc_samp_clinical.RDS")

rm(list=ls())
