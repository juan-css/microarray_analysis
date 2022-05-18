#### Microarray expression analysis
#### Affymetrix

# Author: Juan Carlo Santos e Silva
# Computational Systems BIolopgy Laboratories
# Advisor: Helder I Nakaya
# Date: 18/05/2022

# -- Get data
# -- Quality control
# -- Annotation
# -- Collapse/Summarization
# -- Differential Expression Analysis



########## 1. Get data  --------------

# Baixando dados do estudo GSE54992
# BiocManager::install("GEOquery")
# BiocManager::install("affy")
# BiocManager::install("arrayQualityMetrics")
library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(dplyr)
library(tidyr)
library(readr)
library(biomaRt)
library(data.table)
options(stringsAsFactors = FALSE)

# Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir.create("data")
setwd("./data")

gse_ID <- "GSE54992"
gpl_ID <- "GPL570"
path_gse = paste0('data/', gse_ID, '/')

#baixando os dados do GEO
gse <- getGEO(gse_ID)
show(gse_ID)

# baixando dados brutos do estudo (.CEL)
getGEOSuppFiles(gse_ID)

# alvando os dados de expressão em um arquivo
table_expression <- exprs(gse[[1]])
write.table(table_expression, "GSE54992/table_expression_array.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

# visualizar o phenodata
table_phenodata <- pData(phenoData(gse[[1]]))
write.table(table_phenodata, "GSE54992/table_phenodata.tsv", sep = "\t")

# obtendo dados da plataforma
probe_table <- gse[[1]]@featureData@data
write.table(probe_table, "GSE54992/annotation_platform.tsv", sep = "\t")



########## 2. Quality Control  --------------

# Go back to the parent directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



#### ----- a. Run AQM pre normalization -----####

# Find studies with CEL files (the file needs)
untar(tarfile = paste0(path_gse,'GSE54992_RAW.tar'), 
      exdir = paste0(path_gse,'GSE54992_RAW'))

cel_files <- list.files(path = paste0(path_gse,'GSE54992_RAW'),
                        pattern = ".cel", 
                        ignore.case = TRUE, 
                        recursive = TRUE, 
                        full.names = TRUE)

# Reading the CEL files (create Affybatch file) 
rawdata <- ReadAffy(filenames = cel_files)

# Create a directory to save AQM pre normalization (non-normalized)
dir_aqm_non_norm <- paste0("intermediate/aqm_teste/", gse_ID, "_", gpl_ID, "_AQM_non_norm")

# Running AQM
arrayQualityMetrics(expressionset = rawdata, 
                    outdir = dir_aqm_non_norm, 
                    force = TRUE, do.logtransform = TRUE)

raw_expr <- rawdata@assayData[["exprs"]]

# Check index.html file



### ----- b. Normalization -----####

# https://www.biostars.org/p/69570/
# https://stackoverflow.com/questions/25581769/creating-eset-object-from-preprocessed-expression-matrix

# Using RMA to normalize the data 
expr_norm <- rma(rawdata)



#### ----- c. Run AQM pos normalization -----####
dir_aqm_norm <- paste0("intermediate/aqm_teste/", gse_ID, "_", gpl_ID, "_AQM_norm")

# The expressionset needs to be the file that comes out ReadAffy (it's not just a simple dataframe, it's multiple lists together)
arrayQualityMetrics(expressionset = expr_norm, 
                    outdir = dir_aqm_norm, 
                    force = TRUE)

# Save the normalized expression file, if case:
expr_norm_final <- expr_norm@assayData[["exprs"]]
write.table(expr_norm_final, "intermediate/expr_final_norm.csv", sep = "\t", row.names = TRUE, col.names = TRUE)

############ Attention! ############
# Only remove outliers if the arrayQualityMetrics detects outlier in two or more tests.
# In case you detect this put the name of the samples in the vector outlier_samples
# Else, go to step 3



# ----- d. Remove outliers, renorm and control quality ----- ####
dir_aqm_renorm <- paste0("intermediate/aqm_teste/", gse_ID, "_", gpl_ID, "_AQM_renorm")

# Creating a vector with the outlier samples (see the table in the index)
outlier_samples <- c("GSE54992_RAW/XXXXXX_OUTLIER_SAMPLE_NAME.CEL.gz", 
                     "GSE54992_RAW/XXXXXX_OUTLIER_SAMPLE_NAME.CEL.gz", 
                     "GSE54992_RAW/XXXXXX_OUTLIER_SAMPLE_NAME.CEL.gz")

# Removing the outlier samples from cel_files (see the index and check the sample names)
cel_files2 <- cel_files[!cel_files %in% outlier_samples]

# Get the AffyBatch from cel_files without outliers (pode pular)
rawdata2 <- ReadAffy(filenames = cel_files2)

# Renorm without outliers
expr_renorm <- rma(rawdata2)

# Check the AQM again
arrayQualityMetrics(expressionset = expr_renorm, 
                    outdir = dir_aqm_norm, 
                    force = TRUE)

# Save the renormalized expression data  
expr_renorm_final <- expr_renorm@assayData[["exprs"]]
write.table(expr_renorm_final, "intermediate/expr_final_norm_noOutlier.csv", sep = "\t", row.names = TRUE, col.names = TRUE)

# Clean useless data
rm()



########## 3. Annotation  --------------
# Get the saved probe_table again
probe_table = read.delim("data/GSE54992/annotation_platform.tsv")


# Select Normalized by authors or normalized by you:
# expr = fread("data/GSE54992/table_expression_array.tsv")
expr = read.delim("intermediate/expr_final_norm.csv")



# ----- a. Gene names from authors ----
# Turn probe_names into a column
expr <- cbind(probeName = rownames(expr), expr)
rownames(expr) = NULL

# Filter probe_table to get the gene_ID
colnames(probe_table)
probe_table <- probe_table[,c("ID", "Gene.Symbol")]

# Use left_join to get the names of the genes from probe_table
expr <- expr %>% 
  left_join(probe_table, by = c("probeName" = "ID"))

# Separating using /// of "Gene.Symbol" e colocando em outra coluna ao lado (coluna "Delete")
# deletando a coluna "Delete" com a função select(-Delete)
expr <- expr %>% 
  separate("Gene.Symbol", c("geneName_gse", "Delete"), " /// ") %>% 
  dplyr::select(-Delete)

expr <- expr[,c(1, ncol(expr), 2:(ncol(expr)-1))]


#### ----- 4. Collapse   --------------

source("source.R")

# Collapsing
expr <- collapse.rows(expr = expr, 
                      probe.col = 'probeName', 
                      gene.col = 'geneName_gse', 
                      method = 'maxMean')

expr = expr[expr$geneName_gse != "", ] 
length(unique(expr$geneName_gse))
expr = expr[,-2]

write.table(expr, 
            "intermediate/expr_table_collapsed.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE)



