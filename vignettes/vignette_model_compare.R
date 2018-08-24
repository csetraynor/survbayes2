start_time <- Sys.time()

library(readr)
library(plyr)
library(dplyr)
library(ggplot2)
library(survival)
library(parallel)
library(caret)
library(rsample)
library(VIM)
library(rstan)
library(rstanarm)
theme_set(theme_bw())
devtools::document()
options(mc.cores = parallel::detectCores() )
rstan_options(auto_write = TRUE)

data(lusc)

#convert variable names to lower case
colnames(lusc) <- tolower(colnames(lusc))
#convert all string NA to NA format
lusc <- lusc %>%
  mutate_all(funs(convert_blank_to_na))

lusc %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE,
            sortVars = TRUE, sortCombs = TRUE, plot = FALSE, only.miss = TRUE)

lusc <- lusc %>%
  dplyr::filter(tobacco_smoking_history_indicator != "1" )

lusc <- lusc %>% 
  dplyr::select( patient_id, os_status, os_months,  ajcc_nodes_pathologic_pn, ajcc_tumor_pathologic_pt, age,  icd_10, sex , ajcc_pathologic_tumor_stage, ajcc_metastasis_pathologic_pm, tobacco_smoking_history_indicator ) %>%
  dplyr::mutate(age = as.numeric(age),
                os_months = as.numeric(os_months),
                os_status = as.logical(os_status == "DECEASED"),
                age_d = ifelse(age < 65, 0, 1),
                stage = 0,
                nodes = 0,
                tumor = 0,
                metastasis = ifelse(ajcc_metastasis_pathologic_pm == "M0", 0, 1),
                smoke = 0) 

lusc$stage[grep("III|IV", lusc$ajcc_pathologic_tumor_stage)] <- 1

lusc$nodes[grep("2|3", lusc$ajcc_nodes_pathologic_pn)] <- 1

lusc$tumor[grep("3|4", lusc$ajcc_tumor_pathologic_pt)] <- 1

lusc$smoke[grep("2", lusc$tobacco_smoking_history_indicator)] <- 1

# An International Prognostic Index variable is created which has levels low, medium and high. 
lusc <- lusc %>%
  dplyr::mutate(ipi_n = stage + nodes + 
                  tumor + metastasis,
                ipi = ifelse( ipi_n < 1 , "low", 
                              ifelse( ipi_n < 2 , "medium",
                                      "hihg") ) )

lusc <- lusc %>%
  dplyr::mutate_at(c("stage", "nodes", "tumor", "smoke"), as.factor)

lusc <- lusc[complete.cases(lusc), ] 

preProcValues <- caret::preProcess(lusc[c("age")], method = c("center", "scale") )
trainTransformed <- as.data.frame( predict(preProcValues, lusc[c("age")]) )
colnames(trainTransformed) = c("age_std")
luscs <- cbind(lusc, trainTransformed)
rm(list = c("trainTransformed", "preProcValues")) # but we keep preProcValues

#Train null model

system.time( survbayes.mod1 <- post_surv(x = luscs %>%
                              dplyr::mutate(time = os_months,
                                            status = os_status)) )

#To create step stan function
deparse(survbayes.mod1$formula)
survbayes.mod1$stan_function
#Stepwise selection
system.time( survbayes.step.fit1 <- stan_surv_step(x = luscs %>%
                                                     dplyr::mutate(time = os_months,
                                                                   status = os_status), 
                                                   fit = survbayes.mod1,
                                                   scope = c("smoke", "age_std", "stage", "nodes", "tumor", "metastasis", "s(age_std)"), verbose = TRUE)
)

system.time( survbayes.step.fit2 <- stan_surv_step(x = luscs %>%
                                                     dplyr::mutate(time = os_months,
                                                                   status = os_status), 
                                                   fit = survbayes.mod1,
                                                   scope = c("smoke", "age_std", "ipi", "s(age_std)"), verbose = TRUE)
)

pem.fits <- list(survbayes.step.fit1, survbayes.step.fit2)
pem.fits.dev <- sapply(fits, function(m){
  -2*m[[2]]$estimates["elpd_loo","Estimate"]
});max(unlist(fits.dev))

best.pem.fit <- pem.fits[[ match( max(unlist(pem.fits.dev)), pem.fits.dev ) ]][[1]]
class(best.pem.fit)

## Train Cox model
cox.mod1 <- cox_step(x = luscs %>% 
                       mutate(time = os_months,
                              status = os_status), surv_form = c("smoke", "age_std", "stage", "nodes", "tumor", "metastasis") )

cox.mod2 <- cox_step(x = luscs %>% 
                       mutate(time = os_months,
                              status = os_status), surv_form = c("smoke", "age_std", "ipi", "s(age_std)" ) )

cox.fits <- list(cox.mod1, cox.mod2)
cox.aic <- sapply(cox.fits , function(c){
  -2*c$loglik[[2]] + length( c$coefficients )
})

best.cox.fit <- cox.fits[[ match( max(unlist(cox.aic)), cox.aic ) ]]
class(best.cox.fit)

# Train test validation
# set.seed(99)
# slusc <- resample_stratified(lusc, strata = c( "age_d", "smoke", "ipi"), sizeto = 100)

set.seed(10)
vfold <- rsample::mc_cv(lusc, times = 1, strata = "os_status")$splits$`1`
train_lusc <- rsample::analysis(vfold)
test_lusc <- rsample::assessment(vfold)

preProcValues <- caret::preProcess(train_lusc[c("age")], method = c("center", "scale") )
trainTransformed <- as.data.frame( predict(preProcValues, train_lusc[c("age")]) )
colnames(trainTransformed) = c("age_std")
train_luscs <- cbind(train_lusc, trainTransformed)
rm(trainTransformed) # but we keep preProcValues
testTransformed <- as.data.frame( predict(preProcValues, test_lusc[c("age")]) )
colnames(testTransformed) = c("age_std")
test_luscs <- cbind(test_lusc, testTransformed)
rm(list = c("preProcValues", "testTransformed")) 

pem.train <- post_surv(x = train_luscs %>%
                        dplyr::mutate(time = os_months,
                                      status = os_status),
                      form = best.pem.fit$formula) 

#Train cox model
cox.train <- coxph( best.cox.fit$formula, data = train_luscs %>%
                        dplyr::mutate(time = os_months,
                                      status = os_status) )


## Obtain Brier score

cox.brier <- get_survbrier(test = test_luscs %>%
                dplyr::mutate(time = os_months,
                              status = os_status),
              mod = cox.train
              )
cox.apperror <- cox.brier$AppErr$matrix

pem.brier <- get_survbrier(test = test_luscs %>%
                             dplyr::mutate(time = os_months,
                                           status = os_status),
                           mod = pem.train
)


## Obtain c-index

cox.cindex <- get_cindex(test = test_luscs %>%
                             dplyr::mutate(time = os_months,
                                           status = os_status),
                           mod = cox.train
)

brier.cindex <- get_cindex(test = test_luscs %>%
                           dplyr::mutate(time = os_months,
                                         status = os_status),
                         mod = pem.train
)


## Obtain ROC


cox.roc <- get_survroc(test = test_luscs %>%
                           dplyr::mutate(time = os_months,
                                         status = os_status),
                         mod = cox.train
)

pem.roc <- get_survroc(test = test_luscs %>%
                             dplyr::mutate(time = os_months,
                                           status = os_status),
                           mod = pem.train
)

end_time <- Sys.time()
end_time - start_time

# 
# 
# survbayes.mod1 <- post_surv(x = train_lusc %>%
#                               dplyr::mutate(time = os_months,
#                                             status = os_status))
# 
# loo.mod1 <- loo(survbayes.mod1, cores = detectCores())
# 
# survbayes.mod2 <- post_surv(x = train_lusc %>%
#                               dplyr::mutate(time = os_months,
#                                             status = os_status),
#                             surv_form = c("smoke"))
# 
# timepoints = seq(1 , max(train_lusc$os_months), length.out = 75)
# 
# dtime = diff(c(0, timepoints))
# 
# #create a dataset for simulation
# nosmoke_data <- data.frame(
#   time = timepoints,
#   log_dtime = log(dtime), 
#   smoke = 0
# )
# smoke_data <- data.frame(
#   time = timepoints,
#   log_dtime = log(dtime), 
#   smoke = 1
# )
# 
# sk.plot.frame <- get_plot_frame(mod = survbayes.mod2, long_x = smoke_data, unix = TRUE)  %>%
#   dplyr::mutate(treatment =  "Currently smoking")
# n.plot.frame <- get_plot_frame(mod = survbayes.mod2, long_x = nosmoke_data, 
#                                unix = TRUE) %>%
#   dplyr::mutate(treatment = "Reformed smoking")
# 
# data.plot <- rbind(sk.plot.frame, n.plot.frame)
# 
# # Kaplan and Meier ML estimate
# obs.mortality <- get_km.frame(train_lusc, strata = c("smoke"), time = "os_months", status = "os_status" ) %>%
#   dplyr::mutate_at(vars("smoke"), as.numeric) %>%
#   dplyr::mutate_at(vars("smoke"), as.logical)
# head(obs.mortality)
# 
# p1 <- ggplot2::ggplot(data.plot %>%
#                         dplyr::rename(History = treatment),
#                       aes(time, survmean,
#                           group = History, colour = History)) + 
#   geom_line(mapping=aes(group = History, colour = History),  alpha = 0.8) +
#   geom_step(data= obs.mortality %>%
#               mutate(History = ifelse(smoke, "Currently smoking", "Reformed smoking") ), aes(time,                  surv,  colour = History),  alpha = 0.9)+
#   labs(title = "Smoking LUSC patients",
#        subtitle = "Survival counterfactual plot",
#        caption = "4000 posterior sample size") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   ylab("Survival probability") +
#   xlab("Time (months)") +
#   theme_bw() + theme(legend.position=c(0.7, 0.72)) 
# p1
# 
# 
# #### Overlay plot
# 
# newdata <- gen_long_dat(dat = train_lusc %>% 
#                           mutate(time = os_months,
#                                  status = os_status))
# 
# post <- posterior_linpred(survbayes.mod2, newdata = newdata)
# plot.matrix <- link.surv(post = post, longdata = newdata %>%
#                            mutate(sample_id = numeric_id))
# 
# plot.frame <- get_plot.frame(post = plot.matrix, strata = c("smoke"), obs = rbind(train_lusc %>%                                                                               dplyr::mutate(
#   time = os_months,
#   status = os_status
# )))
# 
# plot.frame$smoke <- as.logical(as.numeric(as.character(plot.frame$smoke)))
# 
# p2 <- ggplot2::ggplot(plot.frame %>%
#                        mutate(History = ifelse(smoke, "Currently smoking", "Reformed smoking") ),
#                      aes(time, postmean, group = History, colour = History)) +
#   geom_line(mapping=aes(group = History, colour = History),  alpha = 0.8)  +
#   # geom_ribbon(aes( ymin = lower,
#   #                  ymax = upper,
#   #                  fill = History), alpha = 0.5, colour=NA) +
#   geom_step(data=obs.mortality  %>%
#               mutate(History = ifelse(smoke, "Currently smoking", "Reformed smoking") ) , aes(time,
#                                      surv,
#                                      colour = History),  alpha = 0.9)+
#   labs(title = "Smoking LUSC patients",
#        subtitle = "Survival counterfactual plot",
#        caption = "4000 posterior sample size") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   ylab("Survival probability") +
#   xlab("Time (months)") +
#   theme_bw() + theme(legend.position=c(0.7, 0.72)) 
# p2
# 
# survbayes.mod3 <- post_surv(x = train_lusc %>%
#                               dplyr::mutate(time = os_months,
#                                             status = os_status),
#                             surv_form = c("smoke", "age_std"))
# survbayes.mod4 <- post_surv(x = train_lusc %>%
#                               dplyr::mutate(time = os_months,
#                                             status = os_status),
#                             surv_form = c("smoke", "age_std", "s(age_std)"))
# 
# 
# loo.mod2 <- loo(survbayes.mod2, cores = detectCores())
# 
# survbayes.mod5 <- post_surv(x = train_lusc %>%
#                               dplyr::mutate(time = os_months,
#                                             status = os_status),
#                             surv_form = c("smoke", "age_std", "stage", "nodes", "tumor", "metastasis"))
# 
# survbayes.mod6 <- post_surv(x = train_lusc %>%
#                               dplyr::mutate(time = os_months,
#                                             status = os_status),
#                             surv_form = c("smoke", "age_std", "ipi"))
# 
# survbayes.mod7 <- post_surv(x = train_lusc %>%
#                               dplyr::mutate(time = os_months,
#                                             status = os_status),
#                             surv_form = c("smoke", "age_std", "stage", "nodes", "tumor", "metastasis", "s(age_std)"),
#                             prior =  hs())
# 
# survbayes.mod8 <- post_surv(x = train_lusc %>%
#                               dplyr::mutate(time = os_months,
#                                             status = os_status),
#                             surv_form = c("smoke", "age_std", "ipi", "s(age_std)"))
# 
# loo.mod1 <- loo(survbayes.mod1, cores = detectCores())
# loo.mod2 <- loo(survbayes.mod2, cores = detectCores())
# loo.mod3 <- loo(survbayes.mod3, cores = detectCores())
# loo.mod4 <- loo(survbayes.mod4, cores = detectCores())
# loo.mod5 <- loo(survbayes.mod5, cores = detectCores())
# loo.mod6 <- loo(survbayes.mod6, cores = detectCores())
# loo.mod7 <- loo(survbayes.mod7, cores = detectCores())
# loo.mod8 <- loo(survbayes.mod8, cores = detectCores())
# 
# print( loo::compare(loo.mod1, loo.mod2, loo.mod3 , loo.mod4, loo.mod5, loo.mod6, loo.mod7, loo.mod8
# ) , digits = 3)
# 
# 
# system.time( survbayes.step.fit1 <- stan_surv_step(x = train_lusc %>%
#             dplyr::mutate(time = os_months,
#                           status = os_status), 
#             fit = survbayes.mod1,
#             scope = c("smoke", "age_std", "stage", "nodes", "tumor", "metastasis", "s(age_std)"), verbose = TRUE)
# )
# 
# system.time( survbayes.step.fit2 <- stan_surv_step(x = train_lusc %>%
#                                                dplyr::mutate(time = os_months,
#                                                              status = os_status), 
#                                              fit = survbayes.mod1,
#                                              surv_form = c("smoke", "age_std", "npi", "s(age_std)"), verbose = TRUE)
# )
#             
# #Cox model stepwise
# 
# cox.mod1 <- cox_step(x = train_lusc %>% 
#                        mutate(time = os_months,
#                               status = os_status), surv_form = c("smoke", "age_std", "stage", "nodes", "tumor", "metastasis") )
# 
# cox.mod2 <- cox_step(x = train_lusc %>% 
#                        dplyr::mutate(time = os_months,
#                               status = os_status), 
#                      surv_form = c("smoke", "age_std", "ipi") )
# 
# 
# ##Create test train splits
# 
# set.seed(100)
# vfold <- rsample::mc_cv(lusc, times = 1, strata = "os_status")$splits$`1`
# train_lusc <- rsample::analysis(vfold)
# test_lusc <- rsample::assessment(vfold)
# 
# ##tdbrier
# 
# timepoint <-  seq(0, max(train_lusc$os_months[train_lusc$os_status]),
#                    length.out = 100L)
# timepoint <-   sort (unique ( test_lusc$os_months[test_lusc$os_status] ) )
# 
# 
# newdat <-  gen_new_frame(dat = test_lusc %>% 
#                           dplyr::mutate(time = os_months,
#                                         status = os_status),
#                          timepoints = timepoint)
# pred_frame <- pred_surv(long_x = newdat, 
#                        mod = survbayes.mod2,
#                        unix = TRUE)
# newdat$surv_mean <- as.vector( apply(pred_frame, 1, mean) ) 
# 
# probs <- get_survProb(newdat = newdat, surv = "surv_mean", patient_id = "patient_id")
# 
# brier.gam <- pec::pec(probs, Surv(time, status) ~ 1,
#                   data = test_lusc %>% 
#                     dplyr::mutate(time = os_months,
#                                   status = os_status),
#                   maxtime = max(timepoints),
#                   exact = FALSE,
#                   exactness = 49L)
# 
# probs <- pec::predictSurvProb(cox.model1, newdata = test_luscs, times = timepoints)
# 
# brier.cox <- suppressMessages( pec::pec(probs, Surv(os_months, os_status) ~ 1,
#                                     data = test_luscs,
#                                     maxtime = max(timepoints),
#                                     exact = FALSE,
#                                     exactness = 99L) )
