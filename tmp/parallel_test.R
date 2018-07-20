library(brms)
library(methods)
library(rstan)
library(survival)
library(magrittr)
library(plyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Read arguments listed on command line
args = (commandArgs(TRUE))
i = as.integer(args[1])

print(i)

#source gen_stan_dat function
source("R/misc.R")


#Load data --------------
mc_samp <- readRDS("data-raw/mc_samp.RDS")


#Bayesian clinical model ----------
surv_form <- c( "age_std", "npi_std", "her2pos", "erpos")

form <- as.formula( paste0(c( "status~1+offset(log_dtime)+s(time)+", paste(surv_form, collapse = "+")), collapse = ""))

# loop object through splits

x <- mc_samp$splits[[i]] #choose splits
x <- as.data.frame(x$data)[x$in_id,]
long_x <- gen_stan_dat(x)
bayes_test <- brm(bf(form),
                  data = long_x, family = poisson(),  cores = 4, seed = 17,
                  iter = 14000, thin = 10, refresh = 0,
                  control = list(adapt_delta = 0.99) ) 

saveRDS(bayes_test, paste0("bayes_model_brms_clinical_full_", i, ".RDS") )


rm(list=ls())

