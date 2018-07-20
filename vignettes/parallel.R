#R script
#load libraries
library(methods)
library(rstan)
library(survival)
library(magrittr)
library(plyr)
library(rstanarm)
library(caret)

#set options for Stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Read arguments listed on command line
args = (commandArgs(TRUE))
i = as.integer(args[1])
print(i)

#source gen_stan_dat function
source("misc.R")

#Load data --------------
mc_samp <- readRDS("data-raw/mc_samp.RDS")

# get data -----

x <- mc_samp$splits[[i]] #choose splits
x <- as.data.frame(x$data)[x$in_id,]


#Apply standardisation to continuos variables

preProcValues <- caret::preProcess(x[c("age_at_diagnosis", "npi")], method = c("center", "scale") )
trainTransformed <- predict(preProcValues, x[c("age_at_diagnosis", "npi")])
x <- cbind(x, trainTransformed)

#Keep info about training data
saveRDS(preProcValues, paste0("trainTransformed_clin", i, ".RDS") )
rm(list = c("preProcValues", "trainTransformed") )


#Create lonf format
long_x <- gen_stan_dat(x)
#Bayesian clinical model ----------
surv_form <- c( "age_std", "npi_std", "her2pos", "erpos")

form <- as.formula( paste0(c( "status~1+offset(log_dtime)+s(time)+", paste(surv_form, collapse = "+")), collapse = ""))

# loop object through splits

bayes_test <- rstanarm::stan_gamm4(form, data = long_x ,
                                   #set likelihood (poisson)
                                   family= poisson(),
                                   #set priors (default)
                                   prior = normal(),
                                   prior_intercept = normal(),
                                   prior_smooth = exponential(autoscale = FALSE),
                                   prior_aux = exponential(),
                                   prior_covariance = decov(),
                                   prior_PD = FALSE,
                                   #arguments passed to sampling algorithm
                                   seed = 7,
                                   iter = 10000, thin = 5, refresh = 0, chains = 4,
                                   control = list(adapt_delta = 0.9)
                                   )

saveRDS(bayes_test, paste0("bayes_model_brms_clinical_full_", i, ".RDS") )

#clear envir
rm(list=ls())



# # Cox clinical model -------------
# surv_form <-  '~ age_std + npi_std + her2pos + erpos'
# cox_test <- lapply(mc_samp$splits, function(x){
#   x <- as.data.frame(x$data)[x$in_id,]
#   form <- as.formula( paste0(c( "Surv(time, status)", surv_form), collapse = ""))
#   coxph(form , data = x)
# })
# mc_samp$cox_clinical <- cox_test
# rm(cox_test)

