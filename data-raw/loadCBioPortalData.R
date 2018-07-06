#Load libraries ---
library(tidyverse)
library(predsurv)

#Download from CBioPortal -----
brca_gene <- readRDS("E:/brca_data/brca_data.RDS")
brca_cna <- readRDS("E:/brca_data/cna_expression.RDS")
gene_matrix <- cbind(brca_gene, brca_cna)
str(gene_matrix)
rm(list = c("brca_cna", "brca_gene"))
#Clinical data -------
library(cgdsr)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

the_study_list = getCancerStudies(mycgds)[25,1]
case_list = getCaseLists(mycgds, the_study_list)[2,1]
clinical_data <-  getClinicalData(mycgds, case_list)

#Preprocess data ------
colnames(clinical_data) <- tolower(colnames(clinical_data))
clinical_data <- tibble::rownames_to_column(clinical_data, var = "patient_id")
clinical_data$patient_id <-gsub("\\.", "-", clinical_data$patient_id)

#Select clinical vars -------
dplyr::glimpse(clinical_data)
brca_data <- clinical_data %>%
  select(patient_id, chemotherapy, radio_therapy, hormone_therapy, breast_surgery,
         npi, os_months, os_status, age_at_diagnosis, er_status, her2_status, threegene,
         tumor_size, tumor_stage, grade, intclust)

bric_data <- left_join(gene_matrix, brca_data, by = "patient_id")
predsurv::plot_na(brca_data[brca_data$patient_id %in% bric_data$patient_id, ])

#Save and back-up
saveRDS(bric_data, "E:/brca_data/bric_data.RDS")


