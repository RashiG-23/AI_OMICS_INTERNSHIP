
input_dir <- "module2b" 
output_dir <- "Results"
files_to_process <- c("DEGS_Data_1.csv", "DEGs_Data_2.csv") 
result_list <- list()

classify_gene <- function(logFC, padj) {
  if (logFC > 1 && padj < 0.05) {
    return("Upregulated")
  } else if (logFC < -1 && padj < 0.05) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}

for (file in files_to_process) {
  
  data <- read.csv(file)
  
  cat("File imported. Checking for missing values...\n")
  
  data$padj[is.na(data$padj)] <- 1
  
  data$status <- mapply(classify_gene, data$logFC, data$padj)
  
  output_file_path <- file.path(output_dir, paste0("DEG_Results", basename(file)))

  write.csv(data, output_file_path, row.names = FALSE)
  
  cat("Processed and saved:", output_file_path, "\n")
}

print(table(data$status))
