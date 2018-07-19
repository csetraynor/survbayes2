cna=readr::read_tsv("brca_metabric/data_CNA.txt")

clidata <- readRDS("brca_data.RDS")

#get cna expression matrix obtain id
cna_expression <- cnadata %>% dplyr::select(- Hugo_Symbol, -Entrez_cna_Id)
Hugo_Symbol <- cnadata %>% dplyr::select(Hugo_Symbol) %>% unlist
id_cna_expression <- colnames(cna_expression)
#from clinical data keep observations with cna expression
#obtain sample with cna expression measurements and match
Y <- clidata[clidata$patient_id %in% id_cna_expression,]
Y <- arrange(Y,  match(Y$patient_id, id_cna_expression))
assertthat::assert_that(all(id_cna_expression == Y$patient_id))

##Transpose cna expression matrix
cna_matrix <- as.matrix(cna_expression, ncol = ncol(cna_expression))
t_cna_matrix <- t(cna_matrix)

rm(cna_matrix)

preProcgm <-  caret::preProcess(t_cna_matrix, method = "knnImpute")
t_cna_matrix <- predict(preProcgm, t_cna_matrix)
rm(preProcgm)

rm(preProcgm)

t_cna_matrix <- model.matrix( ~ ., t_cna_matrix)

##Get dataframe back
t_cna_expression <- as.data.frame(t_cna_matrix)

rm(t_cna_matrix)

colnames(t_cna_expression) <- paste(Hugo_Symbol, "_cna", sep="")

### cbind clinical and genomic data
brca <- cbind(Y, t_cna_expression)

rm(Y)
rm(t_cna_matrix)


#Get brca iclust2
iclust2 <- brca[brca$intclust == 2,]
iclust2$intclust <- NULL
brca$intclust <- NULL
