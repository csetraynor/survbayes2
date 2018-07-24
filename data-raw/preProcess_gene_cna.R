library(predsurv)
library(caret)
library(dplyr)
#Load data
dna_microarray <- readr::read_tsv("E:/brca_data/brca_metabric/data_expression.txt")

# Keep gene names
gene_names <- dna_microarray$Hugo_Symbol

# Delete all namesc
gene_matrix <- dna_microarray[ ,!grepl("Hugo_Symbol|Entrez_Gene_Id", colnames(dna_microarray))]

# Check for missing data
sum(is.na(gene_matrix))

# Keep subject names
patient_id <- colnames(gene_matrix)

# Transpose gene matrix
colnames(gene_matrix) <- NULL
gene_matrix <- t(gene_matrix)
str(gene_matrix)

# Normalise gene matrix
colnames(gene_matrix) <- gene_names
preProcGeneValues <- caret::preProcess(gene_matrix, method = c("center", "scale", "knnImpute") )
gene_matrix <- predict(preProcGeneValues, gene_matrix)
saveRDS(preProcGeneValues, paste0("trainTransformed_clin", i, ".RDS") ); rm(preProcGeneValues)

# Impute Gene matrix
gene_matrix <- impute::impute.knn(gene_matrix, k = 10)$data

# add patient id
gene_matrix <- as.data.frame(gene_matrix)
gene_matrix$patient_id <- patient_id
str(gene_matrix)



##### CNA matrix

cna <- readr::read_tsv("E:/brca_data/brca_metabric/data_CNA.txt")

#get cna names
cna_names <- cna$Hugo_Symbol
cna_names <- paste(cna_names, "_cna", sep="")
# cna matrix
cna_matrix <- cna[ ,!grepl("Hugo_Symbol|Entrez_Gene_Id", colnames(cna))]

# Subject id 
cna_matrix <- cna_matrix[patient_id]
patient_id_cna <- colnames(cna_matrix)
assertthat::assert_that(all.equal(patient_id, patient_id_cna))


##Transpose cna expression matrix
colnames(cna_matrix) <- NULL
cna_matrix <- t(cna_matrix)
colnames(cna_matrix) <- cna_names
##Check zero var

test <- caret::nearZeroVar(cna_matrix)
if(length(test) > 0){
  cna_matrix <- cna_matrix[ ,-test]
}

#Convert to factor
cna_matrix <- as.data.frame(cna_matrix)
cna_matrix <- lapply(cna_matrix, as.factor)
cna_matrix <- do.call(cbind.data.frame, cna_matrix)
#add NA as factor
cna_matrix <- lapply(cna_matrix, addNA)
cna_matrix <- do.call(cbind.data.frame, cna_matrix)

glimpse(cna_matrix)

#Create model matrix
X <- model.matrix( ~ ., cna_matrix[,1:10])[ ,-1]

X <- as.data.frame(X)

for(i in 11:ncol(cna_matrix)){
  
  #In case of more than one level
  if(nlevels(cna_matrix[,i]) > 1){
    coluna <- as.data.frame(cna_matrix[,i])
    colnames(coluna) <- colnames(cna_matrix[i])
    matrixna <- model.matrix(~ ., coluna)[ ,-1]
    matrixna <- as.data.frame(matrixna)
    X <- cbind(X, matrixna)
  }
  #Transform to numeric
  else{
    coluna <- as.numeric(as.factor(cna_matrix[,i]))
    coluna <- as.data.frame(coluna)
    print(coluna)
    colnames(coluna) <- colnames(cna_matrix[i])
    X <-cbind(X, coluna)
  }
}


cna_matrix <- as.data.frame(X)


#add patient id
cna_matrix$patient_id <- patient_id_cna

test <- caret::nearZeroVar(cna_matrix)
if(length(test) > 0){
  cna_matrix <- cna_matrix[ ,-test]
}

saveRDS(cna_matrix, "cna_matrix.RDS")

###Load data
bric_data <- readRDS("E:/brca_data/bric_data.RDS")


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
#bric_data$age_std <- bric_data$age_at_diagnosis - mean(bric_data$age_at_diagnosis) / sd(bric_data$age_at_diagnosis) 
#bric_data$npi_std <- bric_data$npi - mean(bric_data$npi) / sd(bric_data$npi) 
if(assertthat::assert_that(all.equal(bric_data$intclust.y,bric_data$intclust.x )) ){
  bric_data$intclust <- bric_data$intclust.y
} else {
  print("Error not equal!")
}

bric_data_clinical <- bric_data[, c("age_at_diagnosis", "npi", "her2pos", "erpos", "status", "time","intclust" )]

#Transform intclust
unique(bric_data_clinical$intclust)
bric_data_clinical$intclust_fact <-  factor(bric_data_clinical$intclust, levels = c("1","2", "3", "4ER-", "4ER+", "5", "6", "7", "8", "9", "10") )
assertthat::assert_that( all.equal(as.character(bric_data_clinical$intclust_fact), bric_data_clinical$intclust) )
bric_data_clinical$intclust_int <- as.numeric(bric_data_clinical$intclust_fact) 

devtools::use_data(bric_data_clinical, overwrite = TRUE)
saveRDS(bric_data_clinical, "E:/brca_data/bric_data_clinical.RDS")
