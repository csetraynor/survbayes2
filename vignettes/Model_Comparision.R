## Survival analysis vignette -------------------
# if(!require("devtools")) install.packages("devtools")
# if(!require("survbayes2")) devtools::install_github("survbayes2")
# if(!require("iclust2prog")) devtools::install_github("iclust2prog")
devtools::document()
library(iclust2prog)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###Load data
bric_data <- readRDS("/media/mtr/SeagateExpansionDrive/brca_data/bric_data.RDS")

#Preprocess
colnames(bric_data) <- my_replace(colnames(bric_data))

bric_data$mastectomy <-  bric_data$breast_surgery == "MASTECTOMY";
bric_data$chemotherapy <-  bric_data$chemotherapy == "YES";
bric_data$radio_therapy <- bric_data$radio_therapy == "YES";
bric_data$hormone_therapy <- bric_data$hormone_therapy == "YES";
bric_data$her2pos <- bric_data$her2_status == "+";
bric_data$erpos <- bric_data$er_status == "+";
#prepare time and status for pipeline
bric_data$status <- bric_data$os_status == "DECEASED";
bric_data$time <- bric_data$os_months
bric_data$time[bric_data$time == 0] <- 1e-3
#center variables
bric_data$age_std <- bric_data$age_at_diagnosis - mean(bric_data$age_at_diagnosis) / sd(bric_data$age_at_diagnosis) 

bric_data_short <- bric_data[, c("age_std", "npi", "mastectomy",
                         "hormone_therapy", "chemotherapy", "radio_therapy", "her2pos", "erpos", "status", "time")]

devtools::use_data(bric_data_short)
#create samples
set.seed(9666)
mc_samp <- mc_cv(bric_data, strata =  "status", times = 200)

memory.size(50000)
# Fit Cox models --------------------
mc_samp$mod_clin_cox <- pmap(list(mc_samp$splits),
                             function(data){
                               mod_coxfit(x = data,
                                          surv_form = '~ age_std + npi + mastectomy +
                                          hormone_therapy + chemotherapy + radio_therapy +
                                          her2pos + erpos'
                               )
                             })
mc_samp$mod_gen_cox <- pmap(list(mc_samp$splits),
                                 function(data){
                                   mod_coxfit(x = data,
                                              surv_form = iclust2_features$feature,
                                              inits = iclust2_features$coef
                                   )
                                 })

# Fit Bayes models --------------------
#Create fit object 
stan_model_object <- 

mc_samp$mod_clin_bym <- pmap(list(mc_samp$splits),
                             function(splits){
                               pred_sbm(x = splits,
                                        surv_form = c( "age_std", "npi", "mastectomy",
                                                       "hormone_therapy", "chemotherapy", "radio_therapy", "her2pos", "erpos")
                               )
                             })

############### Get Cox Brier -------------
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
library(tidyposterior)
int_brier <- perf_mod(int_brier, seed = 6507, iter = 5000, transform = logit_trans)

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
