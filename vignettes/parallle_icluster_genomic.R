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
data(cna_matrix)
data(gene_matrix)
colnames(cna_matrix) <- my_replace(colnames(cna_matrix))
colnames(gene_matrix) <- my_replace(colnames(gene_matrix))
gene_names <- colnames(gene_matrix)[-match("patient_id", colnames(gene_matrix))]
gene_code <- paste0(rep("gene", length(gene_names)),  seq_along(gene_names))
cna_names <- colnames(cna_matrix)[-match("patient_id", colnames(cna_matrix))]
cna_code <- paste0(rep("cna", length(cna_names)),  seq_along(cna_names))

bric_data <- dplyr::left_join(bric_data_clinical, gene_matrix, by = "patient_id")
bric_data <- dplyr::left_join(bric_data, cna_matrix, by = "patient_id")

rm(list = c("bric_data_clinical", "cna_matrix", "gene_matrix") )
#choose cluster
bric_data <- bric_data[bric_data$intclust_int == c, ]

#create samples
set.seed(7)
mc_samp <- rsample::mc_cv(bric_data, strata =  "status", times = 100)
rm(bric_data)
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

#Apply standardisation to gene _matrix

preProcGeneValues <-  caret::preProcess(x[gene_names],method = c("center", "scale", "knnImpute") )
gene_matrix <- predict(preProcGeneValues, x[gene_names])
x <- x[ ,- match(colnames(gene_matrix), colnames(x))]
x <- cbind(x, gene_matrix)

saveRDS(preProcGeneValues, paste0("trainTransformed_gene_cluster_",c,"_sample_", i, ".RDS") ); 
rm(list = c("preProcGeneValues", "gene_matrix"))

#Create lonf format
long_x <- gen_stan_dat(x)
#Bayesian clinical model ----------

gene_form <- paste(gene_code, collapse = "+")
cna_form <- paste(cna_code, collapse = "+")
clin_form <- paste(c( "age_std", "npi_std", "her2pos", "erpos"), collapse = "+")
form1 <- as.formula("status~1+offset(log_dtime)+s(time)+.-time-log_dtime")
surv_form <- paste(c(gene_form, cna_form, clin_form, form1), collapse = "+" )

reformulate(surv_form, "status")

X <- long_x[ ,c( "age_std", "npi_std", "her2pos", "erpos", "time", "status", "log_dtime")]
# loop object through splits

#Rough guess of global_scale (gs)
nz = 20
z = 70000
N = nrow(long_x)
gs = (nz/z)/N

# library(brms)
bayes_test <- rstanarm::stan_gamm4(status~1+offset(log_dtime)+s(time)+.-time-log_dtime, 
                                   data = X ,
                                   #set likelihood (poisson)
                                   family= poisson(),
                                   #set priors (default)
                                    prior = hs(df = 1, global_df = 1,
                                               global_scale = gs,
                                               slab_df = 4,
                                               slab_scale = 2.5),
                                    prior_intercept = normal(),
                                    prior_smooth = exponential(allow_autoscale = FALSE),
                                    prior_aux = exponential(),
                                    prior_covariance = decov(),
                                    prior_PD = FALSE,
                                   #arguments passed to sampling algorithm
                                   seed = 7,
                                   iter = 1,
                                   # thin = 5,
                                   # refresh = 0,
                                   # control = list(adapt_delta = 0.9),
                                   chains = 1 
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



library(rstanarm)
library(dplyr)

library(rstanarm)
data("womensrole", package = "HSAUR3")
womensrole_gamm4_1 <- rstanarm::stan_gamm4(cbind(agree, disagree) ~ offset(gender) + .,
                                           data = womensrole,
                                           family = binomial(link = "logit"), 
                                           prior = student_t(df = 7), 
                                           prior_intercept = student_t(df = 7))



# from example(gamm4, package = "gamm4"), prefixing gamm4() call with stan_
dat <- mgcv::gamSim(1, n = 400, scale = 2) ## simulate 4 term additive truth
#> Gu & Wahba 4 term additive model
## Now add 20 level random effect `fac'...
dat$fac <- fac <- as.factor(sample(1:20, 400, replace = TRUE))
dat$y <- dat$y + model.matrix(~ fac - 1) %*% rnorm(20) * .5

br <- stan_gamm4(y ~ s(x0) + s(x2) + ., data = dat, random = ~ (1 | fac),
                 chains = 1, iter = 200) # for example speed
