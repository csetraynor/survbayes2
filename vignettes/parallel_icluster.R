#load libraries
library(methods)
library(rstan)
library(survival)
library(magrittr)
library(plyr)
library(rstanarm)
library(rsample)
library(utils)


#Read arguments listed on command line
args = (commandArgs(TRUE))
i = as.integer(args[1])
print(i)

c = as.integer(args[2])
print(c)

#set options for Stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#source gen_stan_dat function
source("misc.R")

#Load data --------------
data(bric_data_clinical)

#choose cluster
bric_data_clinical <- bric_data_clinical[bric_data_clinical$intclust_int == c, ]

#create sample
#create samples
set.seed(7)
mc_samp <- rsample::mc_cv(bric_data_clinical, strata =  "status", times = 100)
rm(bric_data_clinical)
# get data -----

x <- mc_samp$splits[[i]] #choose splits
x <- rsample::analysis(x)

rm(mc_samp)

#Apply standardisation to continuos variables

preProcValues <- caret::preProcess(x[c("age_at_diagnosis", "npi")], method = c("center", "scale") )
trainTransformed <- predict(preProcValues, x[c("age_at_diagnosis", "npi")])
colnames(trainTransformed) = c("age_std", "npi_std")

x <- cbind(x, trainTransformed)

#Keep info about training data
saveRDS(preProcValues, paste0("modelfiles_clust", c, "/preProcTrain_clin_iclust", c, "_sample_", i, ".RDS") )
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
                                 # iter = 10000,
				 # thin = 5,
				 # refresh = 0,
                                 # control = list(adapt_delta = 0.9),
				   chains = 4 
                                   )

saveRDS(bayes_test, paste0("modelfiles_clust", c, "/bayes_model_clinical_iclust_", c, "_sample_", i, ".RDS") )

#clear envir
#rm(list=ls())



# # Cox clinical model -------------
surv_form <-  '~ age_std + npi_std + her2pos + erpos'

form <- as.formula( paste0(c( "Surv(time, status)", surv_form), collapse = ""))
cox_model <- coxph(form , data = x)

saveRDS(cox_model, paste0("modelfiles_clust", c, "/cox_model_clinical_iclust_", c, "_sample_", i, ".RDS") )

rm(cox_model)

