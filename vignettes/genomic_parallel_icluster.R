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

#set options for Stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#source gen_stan_dat function
source("R/misc.R")

#Load data --------------
data("bric_data_clinical")
data("cna_matrix")
data(gene_matrix)
colnames(cna_matrix) <- my_replace(colnames(cna_matrix))
colnames(gene_matrix) <- my_replace(colnames(gene_matrix))
gene_names <- colnames(gene_matrix)[-match("patient_id", colnames(gene_matrix))]
cna_names <- colnames(cna_matrix)[-match("patient_id", colnames(cna_matrix))]

bric_data <- dplyr::left_join(bric_data_clinical, gene_matrix, by = "patient_id")
bric_data <- dplyr::left_join(bric_data, cna_matrix, by = "patient_id")

rm(list = c("bric_data_clinical", "cna_matrix", "gene_matrix") )

#create samples
set.seed(7)
mc_samp_genomic <- rsample::mc_cv(bric_data, strata =  "status", times = 100)
rm(bric_data)

devtools::use_data(mc_samp)
# get data -----

train <- mc_samp$splits[[i]] #choose splits
train <- rsample::analysis(train)

rm(mc_samp)

#Apply standardisation to continuos variables

preProcValues <- caret::preProcess(train[c("age_at_diagnosis", "npi")], method = c("center", "scale") )
trainTransformed <- predict(preProcValues, train[c("age_at_diagnosis", "npi")])
colnames(trainTransformed) = c("age_std", "npi_std")

train <- cbind(train, trainTransformed)

#Keep info about training data
saveRDS(preProcValues, paste0("modelfiles_gene_clust", c, "/preProcTrain_clin_iclust_", c, "_sample_", i, ".RDS") )
rm(list = c("preProcValues", "trainTransformed") )

#Apply standardisation to gene _matrix

preProcGeneValues <-  caret::preProcess(train[gene_names],method = c("center", "scale", "knnImpute") )
geneTransformed <- predict(preProcGeneValues, train[gene_names])
train <- train[ ,- match(colnames(geneTransformed), colnames(train))]
train <- cbind(train, geneTransformed)

saveRDS(preProcGeneValues, paste0("modelfiles_gene_clust", c, "/trainTransformed_gene_iclust_",c,"_sample_", i, ".RDS") ); 
rm(list = c("preProcGeneValues", "geneTransformed"))

#Create lonf format
long_train <- gen_stan_surv(train)

#Bayesian clinical model ----------

X <- long_train[ ,c( "age_std", "npi_std", "her2pos", "erpos", "time", "status", "log_dtime", gene_names, cna_names)]
# loop object through splits

#Rough guess of global_scale (gs)
nz = 100
z = 70000
N = nrow(long_train)
gs = (nz/z)/N

# rstanarm needs check interpret.gam function
# bayes_test <- rstanarm::stan_gamm4(status~1+offset(log_dtime)+s(time)+.-time-log_dtime, 
#                                    data = X ,
#                                    #set likelihood (poisson)
#                                    family= poisson(),
#                                    #set priors (default)
#                                     prior = hs(df = 1, global_df = 1,
#                                                global_scale = gs,
#                                                slab_df = 4,
#                                                slab_scale = 2.5),
#                                     prior_intercept = normal(),
#                                     prior_smooth = exponential(allow_autoscale = FALSE),
#                                     prior_aux = exponential(),
#                                     prior_covariance = decov(),
#                                     prior_PD = FALSE,
#                                    #arguments passed to sampling algorithm
#                                    seed = 7,
#                                    iter = 1,
#                                    # thin = 5,
#                                    # refresh = 0,
#                                    # control = list(adapt_delta = 0.9),
#                                    chains = 1 
# )

parameter_ratio <- (nz/z)
library(brms)
brms_model <- brm(bf(status~1+offset(log_dtime)+s(time)+.-time-log_dtime), 
                     data = X ,
                     #set likelihood
                     family = poisson(),
                     #prior
                     prior = set_prior(horseshoe(df = 3, par_ratio = parameter_ratio) ) ,
                     seed = 7,
                     iter = 4000, 
                     thin = 10,
                     warmup = 1000)
                     
                     


saveRDS(brms_model, paste0("modelfiles_gene_clust", c, "/bayes_brms_model_clinical_iclust_", c, "_sample_", i, ".RDS") )
rm(brms_model)
#clear envir
#rm(list=ls())



# # Cox clinical model -------------

model_fit <- function(dat){

# create predictor matrix
x <- dat[ ,c( "age_std", "npi_std", "her2pos", "erpos", gene_names, cna_names)]

###### Fit models
##################################################
#grouped enet
#p.fac = rep(1, ncol(x))
#p.fac[match(c("npi_std","age_std"), colnames(x))] = 0
#prepare
x <- as.matrix(x)
y <- as.matrix(dat[ ,c( "time", "status") ], ncol = 2)
require(doMC)
registerDoMC(cores=4)
### Apply elastic net with a= 0.8
mod <-  glmnet::cv.glmnet(x, y, family = "cox",
                          grouped = TRUE,
                          alpha = 0.8,
			# penalty.factor = p.fac,
                          parallel = TRUE)
return(mod)
}

cox_model <- model_fit(train)

saveRDS(cox_model, paste0("modelfiles_clust", c, "/cox_model_clinical_iclust_", c, "_sample_", i, ".RDS") )

rm(cox_model)




                            
