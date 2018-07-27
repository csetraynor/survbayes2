


#create samples
bric_data_clinical <- bric_data_clinical[bric_data_clinical$intclust_int == 1, ]

set.seed(7)
mc_samp <- rsample::mc_cv(bric_data_clinical, strata =  "status", times = 100)

library(purrr)
devtools::document()

mc_samp$bayes_model_iclust1 <- lapply(seq_along(1:100), function(i){
  readRDS(paste0("tmp/modelfiles_clust1/bayes_model_clinical_iclust_1_sample_", i , ".RDS") )
})

mc_samp$preProcTrain_iclust1 <- lapply(seq_along(1:100), function(i){
  readRDS(paste0("tmp/modelfiles_clust1/preProcTrain_clin", i , ".RDS") )
})

#Get survival probability
mc_samp$iclust1_pred <- pmap(list(mc_samp$splits, mc_samp$bayes_model_iclust1, mc_samp$preProcTrain_iclust1),
                             function(splits, bmodel, preProc ){
                               surv_pred_bgam(
                                 x = splits,
                                 bgam = bmodel,
                                 preProcValues = preProc
                                 
                               )
                             })

#Calculate Brier score
mc_samp$iclust1_bs2 <- pmap(list(mc_samp$iclust1_pred, mc_samp$splits ),
                             function(pred, splits ){
                               get_bs2(
                                 pred_frame = pred,
                                 x = splits,
                                 surv_form = c("age_std", "npi_std", "her2pos", "erpos")
                               )
                             })
