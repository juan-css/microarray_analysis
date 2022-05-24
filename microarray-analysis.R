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
library(matrixStats)
options(stringsAsFactors = FALSE)

# Plots
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(gridtext)

# Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir.create("data")
setwd("./data")

# Baixando dados do estudo GSE54992
gse_ID = "GSE54992"
gpl_ID = "GPL570"
path_gse = paste0('data/', gse_ID, '/')

#baixando os dados do GEO
gse = getGEO(gse_ID)

# baixando dados brutos do estudo (.CEL)
getGEOSuppFiles(gse_ID)

# Salvando os dados de expressão em um arquivo
table_expression = exprs(gse[[1]])
write.table(table_expression, 
            paste0(gse_ID, "/table_expression_array.tsv"), 
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE)

# visualizar o phenodata
metadata <- pData(phenoData(gse[[1]]))
write.table(metadata, 
            paste0(gse_ID, "/metadata.tsv"), 
            sep = "\t")

# obtendo dados da plataforma
probe_table <- gse[[1]]@featureData@data
write.table(probe_table, 
            paste0(gse_ID, "/annotation_platform.tsv"),
            sep = "\t")



########## 2. Quality Control  --------------

# Go back to the parent directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



#### ----- a. Run AQM pre normalization -----####

# Find studies with CEL files (the file needs)
untar(tarfile = paste0(path_gse, gse_ID, '_RAW.tar'), 
      exdir = paste0(path_gse, gse_ID, '_RAW'))

cel_files = list.files(path = paste0(path_gse, gse_ID, '_RAW'),
                        pattern = ".cel", 
                        ignore.case = TRUE, 
                        recursive = TRUE, 
                        full.names = TRUE)

# Reading the CEL files (create Affybatch file) 
rawdata <- ReadAffy(filenames = cel_files)

# Create a directory to save AQM pre normalization (non-normalized)
dir_aqm_non_norm = paste0("intermediate/aqm_teste/", 
                           gse_ID, 
                           "_", 
                           gpl_ID, 
                           "_AQM_non_norm")

# Running AQM
arrayQualityMetrics(expressionset = rawdata, 
                    outdir = dir_aqm_non_norm, 
                    force = TRUE, 
                    do.logtransform = TRUE)

raw_expr = rawdata@assayData[["exprs"]]

# Check index.html file



### ----- b. Normalization -----####

# https://www.biostars.org/p/69570/
# https://stackoverflow.com/questions/25581769/creating-eset-object-from-preprocessed-expression-matrix

# Using RMA to normalize the data 
expr_norm = rma(rawdata)



#### ----- c. Run AQM pos normalization -----####
dir_aqm_norm = paste0("intermediate/aqm_teste/",
                       gse_ID, 
                       "_",
                       gpl_ID, 
                       "_AQM_norm")

# The expressionset needs to be the file that comes out ReadAffy (it's not just a simple dataframe, it's multiple lists together)
arrayQualityMetrics(expressionset = expr_norm, 
                    outdir = dir_aqm_norm, 
                    force = TRUE)

# Save the normalized expression file, if case:
expr_norm_final = expr_norm@assayData[["exprs"]]
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
dir_aqm_renorm = paste0("intermediate/aqm_teste/", 
                         gse_ID, 
                         "_",
                         gpl_ID, 
                         "_AQM_renorm")

# Creating a vector with the outlier samples (see the table in the index)
outlier_samples = c("GSE54992_RAW/XXXXXX_OUTLIER_SAMPLE_NAME.CEL.gz", 
                     "GSE54992_RAW/XXXXXX_OUTLIER_SAMPLE_NAME.CEL.gz", 
                     "GSE54992_RAW/XXXXXX_OUTLIER_SAMPLE_NAME.CEL.gz")

# Removing the outlier samples from cel_files (see the index and check the sample names)
cel_files2 = cel_files[!cel_files %in% outlier_samples]

# Get the AffyBatch from cel_files without outliers (pode pular)
rawdata2 = ReadAffy(filenames = cel_files2)

# Renorm without outliers
expr_renorm = rma(rawdata2)

# Check the AQM again
arrayQualityMetrics(expressionset = expr_renorm, 
                    outdir = dir_aqm_norm, 
                    force = TRUE)

# Save the renormalized expression data  
expr_renorm_final = expr_renorm@assayData[["exprs"]]
write.table(expr_renorm_final, 
            "intermediate/expr_final_norm_noOutlier.csv",
            sep = "\t",
            row.names = TRUE, 
            col.names = TRUE)

# Clean useless data
rm(list=ls())



########## 3. Annotation  --------------
# Get the saved probe_table again
probe_table = read.delim(paste0("data/", gseID, "/annotation_platform.tsv"))


# Select Normalized by authors or normalized by you:
# expr = fread("data/GSE54992/table_expression_array.tsv")
expr = read.delim("intermediate/expr_final_norm.csv")



# ----- a. Gene names from authors ----
# Turn probe_names into a column
expr = cbind(probeName = rownames(expr), expr)
rownames(expr) = NULL

# Filter probe_table to get the gene_ID
colnames(probe_table)
probe_table = probe_table[,c("ID", "Gene.Symbol")]

# Use left_join to get the names of the genes from probe_table
expr = expr %>% 
  left_join(probe_table, by = c("probeName" = "ID"))

# Separating using /// of "Gene.Symbol" e colocando em outra coluna ao lado (coluna "Delete")
# deletando a coluna "Delete" com a função select(-Delete)
expr = expr %>% 
  separate("Gene.Symbol", c("geneName_gse", "Delete"), " /// ") %>% 
  dplyr::select(-Delete)

expr = expr[,c(1, ncol(expr), 2:(ncol(expr)-1))]


#### ----- 4. Collapse   --------------

source("source.R")

# Collapsing
expr = collapse.rows(expr = expr, 
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
metadata = read.delim("data/", gseID, "/metadata.tsv")

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
expr = expr[, metadata$geo_accession]

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


# Comparing classes (Don't set p.value and lfc to get whole list)
degs.test1 = topTable(fit, coef = "test1", 
                      number = nrow(expr),
                      adjust.method = 'fdr', 
                      p.value=1, # get the whole list
                      lfc=0) # get the whole list
# deg.test2 = topTable(fit, coef = "test2", number = nrow(expr),
#                      adjust.method = 'fdr', p.value=0.05, lfc=log2(1))

#### Define DEGs
degs.test1['deg_status'] = 'no'
degs.test1$deg_status[degs.test1$logFC > 1 & degs.test1$adj.P.Val < 0.05] = "UP"
degs.test1$deg_status[degs.test1$logFC < -1 & degs.test1$adj.P.Val < 0.05] = "DOWN"

dim(degs.test1)[1]
dir.create('results/degs')
write.table(degs.test1, 
            'results/degs/degs_TB_vs_HC.tsv',
            sep = '\t', 
            row.names = T)






#### ----- 5. PLOTs   --------------

# ---- a. Volcano Plot ----

# Set colors
mycolors = c("royalblue3", "red2", "lightgrey")
names(mycolors) = c("DOWN", "UP", "no")

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
degs.test1$label = row.names(degs.test1)
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






# ---- b. Relative Expression of one Gene (Boxplot)----

gene_expr = expr["P2RX7",]
gene_expr = data.frame(
  Class=c(metadata$Class),
  Score=c(as.character(gene_expr[1,]))
)

# find means of sample scores
gene_expr$Score = as.numeric(gene_expr$Score)

class_means = vector()
for (j in unique(gene_expr$Class)) {
  class_means = c(class_means,
                  mean(gene_expr[gene_expr$Class == 'HD',
                                 "Score"]))
}
names(class_means) = unique(gene_expr$Class)
class_means = class_means[order(class_means)]


# make color for each class, with control class as light blue
groups = unique(gene_expr$Class)
palette = c("#86cce0", 
             "#4da566", 
             "#d67048", 
             "#b59519",
             "#FDBF6F",
             "#FF7F00",
             "#CAB2D6",
             "#6A3D9A")

if (length(groups) > length(palette)) {
  palette <- rep(palette, ceiling(length(groups)/length(palette)))
}
groups_coloured = palette[1:length(groups)]

control_lab = 'HC'

if (!missing(control_lab)) {
  if (!(control_lab %in% gene_expr$Class)) {
    
    stop("Please provide control label that features in the sample data")
  } else {
    
    groups_reordered = c(control_lab,
                          as.character(groups[-grep(control_lab,
                                                    groups)]))
    
    names(groups_coloured) = groups_reordered
  }
}


#### Plot scores as boxplot graphs
pdf("results/gene_boxplot.pdf")
ggplot2::ggplot(data = gene_expr,
                ggplot2::aes_string(y = "Score",
                                    x = "Class",
                                    fill = "Class")) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 6
  ) +
  ggplot2::labs(title = NULL,
                x = NULL,
                y = "Normalized expression") +
  ggplot2::theme(legend.position = "null") +
  ggplot2::scale_x_discrete(limits = names(class_means)) +
  ggplot2::geom_jitter(
    shape = 16,
    position = ggplot2::position_jitter(0.1),
    size = 2,
    color = "grey10",
    alpha = 0.7
  ) +
  ggplot2::theme(
    axis.line = ggplot2::element_line(size = 0.8,
                                      linetype = "solid"),
    axis.title =
      ggplot2::element_text(size = 14),
    panel.grid.major =
      ggplot2::element_line(linetype = "blank"),
    panel.grid.minor =
      ggplot2::element_line(linetype = "blank"),
    plot.title = ggplot2::element_text(size = 18),
    panel.background =
      ggplot2::element_rect(fill = "white"),
    axis.text.x =
      ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.2),
    axis.title.y = 
      ggplot2::element_text(vjust = 2),
    legend.text =
      ggplot2::element_text(size = 14),
    axis.text = ggplot2::element_text(size = 14)
  ) +
  ggplot2::scale_fill_manual(values = groups_coloured)
dev.off()






# ---- c. Heatmap samples for some genes (Z-score) ----

# Set the genes of interest
gene_list = c(
  'P2RX1',
  'P2RX2',
  'P2RX3',
  'P2RX4',
  'P2RX5',
  'P2RX6',
  'P2RX7'
)

# Check if your genes has expression signal
check = gene_list[!gene_list %in% row.names(expr)]
if (length(check) == 0) {
  print('All target genes are in expression data')
} else{
  print(c('Please check the following genes:', check))
}

# Check if your genes has expression signal
check = gene_list[!gene_list %in% row.names(degs.test1)]
if (length(check) == 0) {
  print('All target genes are in DGE data')
} else{
  print(c('Please check the following genes:', check))
}

# filter their expression to another table
expr_Heat = {}
for(i in gene_list){
  temp = expr[row.names(expr) == i,]
  expr_Heat = rbind(expr_Heat, temp)
}

# Calculate Z-score for each sample
expset_zscore = (expr_Heat-rowMeans(expr_Heat))/(rowSds(as.matrix(expr_Heat)))


# Set Heatmap colors for logFC 
col_logFC = 
  colorRamp2(c(-2, -1, 0, 1, 2),
             c('#217847', '#87dead', 'lightgrey', '#fec44f', '#d95f0e')
)

# Set Heatmap colors for significance
col_adjP =
  colorRamp2(c(0, 2, 4, 6, 8, 10),
             c(
               '#ffffff',
               '#d9d9d9',
               '#b3b3b3',
               '#8c8c8c',
               '#666666',
               '#404040'
             ))

# Set DEGs table
deg_Heat = {}
for(i in gene_list){
  temp = degs.test1[row.names(degs.test1) == i,]
  deg_Heat = rbind(deg_Heat, temp)
}

# Set DEG annotation
right_notes = rowAnnotation(
  logFC = deg_Heat$logFC,
  log10P.adjust = -log10(deg_Heat$adj.P.Val),
  col = list(logFC = col_logFC,
             log10P.adjust = col_adjP),
  simple_anno_size = unit(2, 'mm')
)

ht = Heatmap(
  as.matrix(expset_zscore),
  row_order = sort(gene_list, decreasing = T), #DEFINIR A ORDEM DO HEATMAP ROW
  row_km = 2, #CLUSTERIZAR O HEATMAP ROW
  column_gap = unit(1.5, 'mm'),
  row_gap = unit(1.5, 'mm'),
  rect_gp = gpar(col = "white", lwd = 1.9),
  row_names_gp = gpar(fontsize = 12, fontface = 'italic'),
  column_names_gp = gpar(fontsize = 12),
  name = 'sd from the mean',
  column_dend_reorder = F,
  column_split = factor(c(rep(1,6),rep(2,9)), levels = as.character(c(1,2))),
  right_annotation = right_notes,
  show_parent_dend_line = F
)

pdf("results/heatmap_by_samples.pdf")
draw(ht)
dev.off()

