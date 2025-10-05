################# AI-OMICS INTERNSHIP MODULE II ASSIGNMENT #############

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

BiocManager::install(c("GEOquery","affy","arrayQualityMetrics","hgu133plus2.db","AnnotationDbi"))
BiocManager::install("limma")

install.packages("dplyr")

library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(dplyr)
library(limma)
library(AnnotationDbi)
library(hgu133plus2.db)


gse_data <- getGEO("GSE118782", GSEMatrix = TRUE)

phenotype_data <-  pData(gse_data$GSE118782_series_matrix.txt.gz)

untar("GSE118782_RAW.tar", exdir = "CEL_Files")

raw_data <- ReadAffy(celfile.path = "CEL_Files")

raw_data   

arrayQualityMetrics(expressionset = raw_data,
                    outdir = "results_data",
                    force = TRUE,
                    do.logtransform = TRUE)

normalized_data <- rma(raw_data)
 
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "raw_data/normalized_data",
                    force = TRUE)

is.data.frame(normalized_data)
prepared_data <- as.data.frame(normalized_data)
processed_data <- prepared_data[-27, ]#indexing done to remove one outlier row
processed_data_trans <- t(processed_data)
processed_data <- as.data.frame(processed_data_trans)
dim(processed_data)   

row_median <- rowMedians(as.matrix(processed_data))

hist(row_median,
     breaks = 8, #reduced breaks to improve visualization
     freq = FALSE,
     main = "Median Intensity Distribution")


threshold <- 2.0 #based on histogram analysis
abline(v = threshold, col = "black", lwd = 2) 


indx <- row_median > threshold 
filtered_data <- processed_data[indx, ]

boxplot (processed_data,
         main = "Box Plot showing data after normalisation",
         xlab = "sample_name")

str(processed_data)
strdx.pca <- prcomp(processed_data[,c(1:39)],
                   center = TRUE,
                   scale. = TRUE)
summary(strdx.pca)
install.packages("ggfortify")
library(ggfortify)
strdx.pca.plot <- autoplot(strdx.pca,
                          data = processed_data)
strdx.pca.plot
var_coords <- as.data.frame(strdx.pca$rotation)
var_coords$Variable <- rownames(var_coords)
library(ggplot2)
ggplot(var_coords, aes(x = PC1, y = PC2, label = Variable)) +
  geom_point() +
  geom_text(vjust = 3, size = 1) +
  labs(title = "PCA Plot of data after normalization", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal() 

phenotype_data <- phenotype_data[-27, ] 
colnames(filtered_data) <- rownames(phenotype_data)

processed_data <- filtered_data 

class(phenotype_data$`sample type:ch1`) 

groups <- factor(phenotype_data$`sample type:ch1`,
                 levels = c("Control", "Breast Cancer"),
                 label = c("normal", "cancer"))

class(groups)
levels(groups)
