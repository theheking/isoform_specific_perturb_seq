library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(ggpubr)
library(rstatix)
library(tidyr)
# Define directory and file paths
dir <- 'alt-prom-crispr-fiveprime/files/whippet_dedup'
files <- list.files(dir, pattern = ".gz", full.names = TRUE)
table_name <- file.path(dir, "all_data.txt")

# Print the list of files
print(files)



# Function to process each file
process_file <- function(file) {
  data <- read.table(file, header = TRUE, sep = "\t", row.names = NULL)
  colnames(data) <- c("Isoform", "Tpm", "Read_Counts")
  data <- data %>% select(Isoform, Tpm)
  colnames(data)[2] <- paste("Tpm", gsub(".deduplicated.isoform.tpm.gz", "", basename(file)), sep = "_")
  #reset rownames to Isoform 
  return(data)
}

# Read or create the table
if (file.exists(table_name)) {
  all_data <- read.table(table_name, header = TRUE, sep = "\t", row.names = NULL)
} else {
  #map the fiunction process to every element in list of files and merge the otuputs to simplest df use applyr
  all_data <- lapply(files, process_file)
  #flatten and merge on isforom column
  all_data <- Reduce(function(x, y) merge(x, y, by = "Isoform", all = TRUE), all_data)
  #instead of cbind use merge on Isoform
  write.table(all_data, table_name, sep = "\t", row.names = FALSE, quote = FALSE)
}
#differential expression usibng deseq
inDir = system.file("extdata", package="pasilla")
countFiles = list.files(inDir, pattern="fb.txt$", full.names=TRUE)
basename(countFiles)
flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
basename(flattenedFile)
#subsample the counts in Tpm_SLTM_MP column randomly to make Tpm_SLTM_MP_1 Tpm_SLTM_MP_2 Tpm_SLTM_MP_3
set.seed(123)
all_data <- all_data %>% mutate(Tpm_SLTM_MP_1 = Tpm_SLTM_MP + rnorm(n(), 0, 1))
#for every row in the dataframe randopmyl 



