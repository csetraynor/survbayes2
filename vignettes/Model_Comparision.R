if(!require("devtools")) install.packages("devtools")
if(!require("iclust2prog")) devtools::install_github("iclust2prog")

library(glmnet)
library(purrr)
library(dplyr)
library(tidyr)
library(rsample)
# library(tidyposterior)
library(gridExtra)
library(bindrcpp)
library(iclust2prog)
theme_set(theme_bw())
library(iclust2prog)
library(predsurv)
library(magrittr)
library(dplyr)
library(forcats)
library(tidyr)
library(purrr)
library(modelr)
library(tidybayes)
library(ggplot2)
library(ggstance)
library(ggridges)
library(rstan)
library(rstanarm)
library(survbayes2)
import::from(LaplacesDemon, invlogit)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
devtools::document()

###Load data
data("ic2_lasso_gwt")
data("iclust2_glmnet")
#Plot glmnet fits
glmnet::plot.cv.glmnet(ic2_lasso_gwt)

#Extract features, see my_replace in utils
iclust2_features <- extract_features(iclust2_glmnet)
iclust2_features$feature <-my_replace(iclust2_features$feature)

iclust2_features_wt <- extract_features(ic2_lasso_gwt)
iclust2_features_wt$feature <- iclust2prog::my_replace(iclust2_features_wt$feature)

############ Survival analysis vignette

data("ic2surv")
ic2surv <- ic2surv %>%
  dplyr::mutate(mastectomy = breast_surgery == "MASTECTOMY",
                crtherapy = (chemotherapy == "YES" | radio_therapy == "YES") )
set.seed(9666)
mc_samp <- mc_cv(ic2surv, strata =  "status", times = 200)

cens_rate <- function(x) mean(analysis(x)$status == 1)
summary(map_dbl(mc_samp$splits, cens_rate))

############### Create models
mc_samp$mod_clin_cox <- pmap(list(mc_samp$splits),
                            function(data){
                              mod_coxfit(x = data,
                                    surv_form = '~ age_std + npi'
                                    )
                            })
mc_samp$mod_clin_gen_cox <- pmap(list(mc_samp$splits),
                         function(data){
                           mod_coxfit(x = data,
                                   surv_form = iclust2_features$feature,
                                   inits = iclust2_features$coef
                           )
                         })
mc_samp$mod_clin_gen_treat_cox <- pmap(list(mc_samp$splits),
                             function(data){
                               mod_coxfit(x = data,
                                          surv_form = c( iclust2_features$feature, "mastectomy", "hormone_therapy"),
                                          inits = c(iclust2_features$coef,rep(0,2) )
                               )
                             })



############### Get Brier -------------
mc_samp$brier_clin_cox  <- pmap(list(mc_samp$splits, mc_samp$mod_clin_cox),
                                    function(data, model){
                                      get_tdbrier(data = data,
                                                  mod = model
                                      )
                                    })
mc_samp$brier_clin_gen_cox  <- pmap(list(mc_samp$splits, mc_samp$mod_clin_gen_cox),
                            function(data, model){
                              get_tdbrier(data = data,
                                          mod = model
                              )
                            })
mc_samp$brier_gen_treat_cox <- pmap(list(mc_samp$splits, mc_samp$mod_clin_gen_treat_cox),
                              function(data, model){
                                get_tdbrier(data = data,
                                            mod = model)
                              })


###integrate Brier
mc_samp$'C' <- map_dbl(mc_samp$brier_clin_cox, integrate_tdbrier)
mc_samp$'CG' <- map_dbl(mc_samp$brier_clin_gen_cox, integrate_tdbrier)
mc_samp$'CGwT' <- map_dbl(mc_samp$brier_gen_treat_cox, integrate_tdbrier)
mc_samp$Null <- map_dbl(mc_samp$brier_clin_cox, integrate_tdbrier_reference)
mc_samp$'Bayes_clin' <- map_dbl(mc_samp$mod_clin_bayes, integrate_tdbrier)
mc_samp$'Bayes_gene' <- map_dbl(mc_samp$mod_clin_gen_bayes, integrate_tdbrier)


int_brier <- mc_samp %>%
  dplyr::select(-matches("^mod"), -starts_with("brier"))


resamples <- gather(int_brier) %>%
  dplyr::mutate(statistic = transform$func(statistic))
make_stancode(statistic ~  model + (1 | id),
              data = resamples)


mod <- stan_glmer(statistic ~  model + (model + 0 | id),
                  data = resamples, ...)
int_brier %>%
  dplyr::select(-splits) %>%
  gather() %>%
  ggplot(aes(x = statistic, col = model, fill = model)) +
  geom_line(stat = "density") +
  theme_bw() +
  theme(legend.position = "top")

int_brier <- perf_mod(int_brier, seed = 6507, iter = 5000, transform = logit_trans,
                      hetero_var = TRUE)

pdf <- ggplot(tidy(int_brier)%>%
                filter(model != "Reference")) +
  theme_bw() +
  ylab("Posterior probability of BS")
pdf
ibrier_tab <- tidy(int_brier) %>%
  filter(model != "Reference") %>%
  group_by(model) %>%
  summarise(lower = quantile(posterior, 0.05),
            mean = mean(posterior),
            upper = quantile(posterior, 0.95))
as.data.frame(ibrier_tab) %>% mutate_all(my_round)

comparisons <- contrast_models(
  int_brier,
  list_1 = rep("CG", 3),
  list_2 = c( "Null", "CGwT", "C"),
  seed = 20
)

compare <- ggplot(comparisons, size =  0.01) +
  theme_bw()
compare <- compare + ylab("") + xlab("DeltaBS")

diff_tab <- summary(comparisons, size = 0.01) %>%
  dplyr::select(contrast, starts_with("pract"))
diff_tab

ibrier_Tab <- post_tab(diff_tab, ibrier_tab)
ibrier_Tab <- ibrier_Tab %>% mutate_all(my_round)
pdf("MC_Results.pdf", 7 ,5)
grid.arrange(pdf, compare, nrow = 1)
dev.off()


# Bayesian Model -------------------

genomic_features <- iclust2_features$feature[-match(c("age_std", "npi"), iclust2_features$feature)]

genomic_prior <- c(set_prior(horseshoe(df = 3, par_ratio = 0.1)))
mc_samp$mod_clin_bayes <- pmap(list(mc_samp$splits),
                               function(data){
                                 get_brier_bGAMM(x = data,
                                                 surv_form = c("age_std", "npi")
                                                 
                                 )
                               })

mc_samp$mod_clin_gen_bayes <- pmap(list(mc_samp$splits),
                                   function(data){
                                     get_brier_bGAMM(x = data,
                                                     surv_form = iclust2_features$feature,
                                                     prior = genomic_prior
                                     )
                                   })
mc_samp$mod_clin_gen_treat_bayes <- pmap(list(mc_samp$splits),
                                         function(data){
                                           get_brier_bGAMM(x = data,
                                                           surv_form = c( iclust2_features_wt$feature, "chemotherapy*radio_therapy", "mastectomy", "hormone_therapy"),
                                                           prior = genomic_prior
                                           )
                                         })
