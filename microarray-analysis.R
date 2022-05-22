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

# Data analysis
# BiocManager::install("GEOquery")
library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(dplyr)
library(tidyr)
library(readr)
library(biomaRt)
library(data.table)
library(limma)
options(stringsAsFactors = FALSE)

# Plots
library(ggplot2)
library(ggrepel)

# Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir.create("data")
setwd("./data")

# Baixando dados do estudo GSE54992
gse_ID <- "GSE54992"
gpl_ID <- "GPL570"
path_gse = paste0('data/', gse_ID, '/')

#baixando os dados do GEO
gse <- getGEO(gse_ID)

# baixando dados brutos do estudo (.CEL)
getGEOSuppFiles(gse_ID)

# Salvando os dados de expressão em um arquivo
table_expression <- exprs(gse[[1]])
write.table(table_expression, 
            "GSE54992/table_expression_array.tsv", 
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE)

# visualizar o phenodata
metadata <- pData(phenoData(gse[[1]]))
write.table(metadata, 
            "GSE54992/metadata.tsv", 
            sep = "\t")

# obtendo dados da plataforma
probe_table <- gse[[1]]@featureData@data
write.table(probe_table, 
            "GSE54992/annotation_platform.tsv",
            sep = "\t")



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
dir_aqm_non_norm <- paste0("intermediate/aqm_teste/", 
                           gse_ID, 
                           "_", 
                           gpl_ID, 
                           "_AQM_non_norm")

# Running AQM
arrayQualityMetrics(expressionset = rawdata, 
                    outdir = dir_aqm_non_norm, 
                    force = TRUE, 
                    do.logtransform = TRUE)

raw_expr <- rawdata@assayData[["exprs"]]

# Check index.html file



### ----- b. Normalization -----####

# https://www.biostars.org/p/69570/
# https://stackoverflow.com/questions/25581769/creating-eset-object-from-preprocessed-expression-matrix

# Using RMA to normalize the data 
expr_norm <- rma(rawdata)



#### ----- c. Run AQM pos normalization -----####
dir_aqm_norm <- paste0("intermediate/aqm_teste/",
                       gse_ID, 
                       "_",
                       gpl_ID, 
                       "_AQM_norm")

# The expressionset needs to be the file that comes out ReadAffy (it's not just a simple dataframe, it's multiple lists together)
arrayQualityMetrics(expressionset = expr_norm, 
                    outdir = dir_aqm_norm, 
                    force = TRUE)

# Save the normalized expression file, if case:
expr_norm_final <- expr_norm@assayData[["exprs"]]
write.table(expr_norm_final,
            "intermediate/expr_final_norm.csv",
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE)

############ Attention! ############
# Only remove outliers if the arrayQualityMetrics detects outlier in two or more tests.
# In case you detect this put the name of the samples in the vector outlier_samples
# Else, go to step 3



# ----- d. Remove outliers, renorm and control quality ----- ####
dir_aqm_renorm <- paste0("intermediate/aqm_teste/", 
                         gse_ID, 
                         "_",
                         gpl_ID, 
                         "_AQM_renorm")

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
write.table(expr_renorm_final, 
            "intermediate/expr_final_norm_noOutlier.csv",
            sep = "\t",
            row.names = TRUE, 
            col.names = TRUE)

# Clean useless data
rm(list=ls())



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

# Clean environment
rm(list=ls())



#### ----- 5. DEGs   --------------

# ---- a. Set-up expression and metadata ----

# Filter the samples to the samples you want
metadata = read.delim("data/GSE54992/metadata.tsv")

# Select the class you want to compare
table(metadata$disease.state.ch1)
classes_to_deg = c("healthy donor", "tuberculosis")
metadata = metadata[metadata$disease.state.ch1 %in% classes_to_deg,]

# Set expression table names
expr = read.delim("intermediate/expr_table_collapsed.tsv", row.names = 1)
colnames(expr) = gsub("\\_.*","", colnames(expr))

# Filter expression table with the metadata
# VERY IMPORTANT!! ALWAYS DO THIS!!
print(paste("Total number of samples:",
            nrow(metadata[metadata$geo_accession %in% colnames(expr),])))

expr = expr[,colnames(expr) %in% metadata$geo_accession]

# reorder expression matrix according to the order of the samplesinfo object
# VERY IMPORTANT!! ALWAYS DO THIS!!
expr <- expr[, metadata$geo_accession]

# Check if all are in the same order
all(metadata$Sample %in% colnames(expr))

# Set the classes you need to compare
metadata$Class = metadata$source_name_ch1

# Remove useless information
metadata$Class = gsub("PBMC from ", "", metadata$Class)

# Classes
table(metadata$Class)

# ---- b. Limma ----

# Design matrix experiment
samples = metadata$Class
samples = factor(samples)
samples

design.mat = model.matrix(~0+samples)
colnames(design.mat) = levels(samples)
design.mat

# Contrast Matrix
contrast.mat = makeContrasts(
  test1 = TB - HC, # deg for treatment1
  # test2 = TB - LTB, # In case of compare a third class
  levels = design.mat
)

contrast.mat

# Fit Limmma
# Fit linear model to estimate T, N for each gene
fit = lmFit(expr, design.mat)

# Fit linear model to estimate a set of contrast, e.g. T-N 
fit = contrasts.fit(fit, contrast.mat)

# Given a microarray linear model fit, compute moderate t-statistics, 
# moderate F-statistic and log-odds of differential expression by empirical Bayes 
# moderation of the standard errors towards a common value.
fit = eBayes(fit)


# Comparing classes
degs.test1 = topTable(fit, coef = "test1", 
                      number = nrow(expr),
                      adjust.method = 'fdr', 
                      p.value=0.05, 
                      lfc=log2(1))
# deg.test2 = topTable(fit, coef = "test2", number = nrow(expr),
#                      adjust.method = 'fdr', p.value=0.05, lfc=log2(1))


dim(degs.test1)[1]
dir.create('results/degs')
write.table(degs.test1, 
            'results/degs/degs_TB_vs_HC.tsv',
            sep = '\t', 
            row.names = T)




#### ----- 5. PLOTs   --------------

# ---- a. Volcano Plot ----

#### Define DEGs
degs.test1['deg_status'] <- 'no'
degs.test1$deg_status[degs.test1$logFC > 1 & degs.test1$adj.P.Val < 0.05] <- "UP"
degs.test1$deg_status[degs.test1$logFC < -1 & degs.test1$adj.P.Val < 0.05] <- "DOWN"

# Set colors
mycolors <- c("royalblue3", "red2", "lightgrey")
names(mycolors) <- c("DOWN", "UP", "no")

# Theme
cleanup = 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    # legend.key = element_rect(fill = "white"),
    axis.title =
      ggplot2::element_text(size = 12),
    axis.text = 
      ggplot2::element_text(size = 12),
    legend.text =
      ggplot2::element_text(size = 12)
  )

# Set the labels
degs.test1$label <- row.names(degs.test1)
degs.test1$label[degs.test1$deg_status == "no"] <- NA

# Plot
pdf("results/degs/degs_volcano.pdf")
ggplot(degs.test1, aes(
  x = logFC,
  y = -log10(adj.P.Val),
  col = deg_status,
  label = label
)) +
  geom_point() +
  geom_text_repel(
    max.overlaps = 30,
    box.padding = 1,
    segment.color = "lightgrey"
  ) +
  labs(x = "log2(Fold-Change)", y = "-log10(P.adjusted)") +
  # geom_vline(xintercept=c(-1, 1), col="grey", linetype="dashed") +
  # geom_hline(yintercept=-log10(0.05), col="grey", linetype="dashed") +
  scale_colour_manual(values = mycolors) +
  theme_minimal() +
  cleanup
dev.off()


