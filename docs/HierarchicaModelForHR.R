devtools::install_github("csetraynor/condma")
library(condma)
library(cgdsr)
library(dplyr)
library(survival)
library(predsurv)
library(metafor)
library(lme4)
library(brms)
library(brmstools)
library(tidyverse)
library(ggplot2)
theme_set(theme_bw())
devtools::document()
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

the_study_list = getCancerStudies(mycgds)[25,1]
case_list = getCaseLists(mycgds, the_study_list)[2,1]
clinical_data <-  getClinicalData(mycgds, case_list)

colnames(clinical_data) <- tolower(colnames(clinical_data))
clinical_data <- tibble::rownames_to_column(clinical_data, var = "patient_id")
dplyr::glimpse(clinical_data)

clinical_data$patient_id <- gsub("\\.", "-", clinical_data$patient_id )

clinical_data <- clinical_data[complete.cases(clinical_data), ]
clinical_data <- std_dat(clinical_data)

treatment <- clinical_data %>%
  dplyr::mutate(time = os_months, status = os_status == "DECEASED") %>%
  dplyr::select(patient_id, time, status, chemotherapy, intclust, npi, radio_therapy, age_at_diagnosis )


treat.splits <- split(treatment, as.factor(treatment$intclust) )
coxph.splits <- lapply(seq_along(treat.splits), function(x){
  cox.fit <- coxph(Surv(time, status) ~ chemotherapy + npi + age_at_diagnosis , data = treat.splits[[x]]) 
  data.frame( hr = cox.fit$coefficients[1],
              se = summary(cox.fit)$coefficients[1,3],
              iclust = names(treat.splits)[[x]]
  )
})
coxph.frame <- do.call(rbind, coxph.splits)
coxph.frame <- coxph.frame[c(1,3:10,2), ]


#akin to random effect
bm <- condma::ic_meta(coxph.frame %>%
                          mutate(yi = hr,
                                 sei = se,
                                 study = iclust), hetero_var = TRUE, iter = 5000, warmup = 2000, cores = 4)

bm.out <- bm$stan
pdf("HierarchicalModel.pdf")
brmstools::forest(bm.out,
                  show_data = TRUE,
                  av_name = "Effect size",
                  sort = TRUE,
                  theme_forest = F) +
  xlab("log hazard ratio") +
  ylab("Integrative cluster") +
  xlim(c(-2,2)) +
  geom_vline(xintercept = 0, alpha = 0.5)+
  labs(title="Bayesian Hierarchical Modeling of survival outcome",
       subtitle="Comparison of treatment effects ")
dev.off()
