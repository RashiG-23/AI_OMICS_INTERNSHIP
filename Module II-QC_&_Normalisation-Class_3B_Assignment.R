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
                    outdir = "Results_data",
                    force = TRUE,
                    do.logtransform = TRUE)

normalized_data <- rma(raw_data)
 
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "Results_data/Normalized_data",
                    force = TRUE)

is.data.frame(normalized_data)
prepared_data <- as.data.frame(normalized_data)
processed_data <- prepared_data[-27, ] #indexing done to remove one outlier row

dim(processed_data)   

row_median <- rowMedians(as.matrix(processed_data))

hist(row_median,
     breaks = 7, #reduced breaks to improve visualization
     freq = FALSE,
     main = "Median Intensity Distribution")


threshold <- 1.52 #based on histogram analysis
abline(v = threshold, col = "black", lwd = 2) 


indx <- row_median > threshold 
filtered_data <- processed_data[indx, ]

transformed_filtered_data <- t(filtered_data) #transposed filtered_data to make sure that rows of phenotype_data match column of filtered_data
filtered_data <- as.data.frame(transformed_filtered_data)
phenotype_data <- phenotype_data[-27, ] 

colnames(filtered_data) <- rownames(phenotype_data)

processed_data <- filtered_data 

class(phenotype_data$`sample type:ch1`) 

groups <- factor(phenotype_data$`sample type:ch1`,
                 levels = c("Control", "Breast Cancer"),
                 label = c("normal", "cancer"))

class(groups)
levels(groups)
