library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(ggpubr)s
library(rstatix)
library(tidyr)

# Define directory and file paths
dir <- 'alt-prom-crispr-fiveprime/files/whippet_dedup'
files <- list.files(dir, pattern = ".gz", full.names = TRUE)
table_name <- file.path(dir, "all_data.txt")
process_file <- function(file) {
  data <- read.table(file, header = TRUE, sep = "\t", row.names = NULL)
  colnames(data) <- c("Isoform", "Tpm", "Read_Counts")
  data <- data %>% select(Isoform, Tpm)
  colnames(data)[2] <- paste("Tpm", gsub(".deduplicated.isoform.tpm.gz", "", basename(file)), sep = "_")
  #reset rownames to Isoform 
  return(data)
}

umi_sum <- function(file) {
  data <- read.table(file, header = TRUE, sep = "\t", row.names = NULL)
  colnames(data) <- c("Gene", "UMI Count")
  new_line <- c(basename(file), sum(data$`UMI Count`))
  #add new li
  return(new_line) }

cellcount_sum <- function(file) {
  data <- read.table(file, header = TRUE, sep = "\t", row.names = NULL)
  #create aline where they count the numbner of lines
  new_line <- c(basename(file), nrow(data))
  return(new_line) 
}
# Print the list of files
print(files)



#get the csv eery file in the folder check the number of lines
cellcount <- "alt-prom-crispr-fiveprime/files/pseudobulk_bam/cellbarcode_list_ideal"
#read every file in the folder and sum the second column tio create new df
files <- list.files(cellcount, pattern = ".csv", full.names = TRUE)
cellcount_total <-  t(as.data.frame(lapply( files,cellcount_sum)))
#change to number 
cellcount_total <- as.data.frame(cellcount_total) %>% rename("Cell Count" = "V2") %>% mutate("Cell Count" = as.numeric(`Cell Count`))
#removde the negative control
cellcount_total <- cellcount_total %>% filter(V1 != "Negative.csv")
ggboxplot(cellcount_total, y="Cell Count", fill = "blue", alpha = 0.6) +
  theme(legend.position = "none") +
  labs(x = " Cell Count", y = "Density") 

#read every file in the ffolder and sum the second column tio create new df
umi_folder <- "alt-prom-crispr-fiveprime/files/pseudobulk_bam/umi_count/"
files <- list.files(umi_folder, pattern = ".txt", full.names = TRUE)
umi_total <-  t(as.data.frame(sapply( files,umi_sum)))
#change to number
umi_total <- as.data.frame(umi_total) %>% rename("UMI Count" = "V2") %>% mutate("UMI Count" = as.numeric(`UMI Count`))
# Function to process each file
ggboxplot(umi_total, y="UMI Count", fill = "pink", alpha = 0.6)
  theme(legend.position = "none") +
  labs(x = " UMI Count", y = "Density") 

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

#sum every column and plot rowsums
rowsums <- all_data %>% select(-Isoform) %>% rowSums() %>% as.data.frame() %>%
  rename("Rowsums" = ".") %>% arrange(desc(Rowsums)) %>% #add 1 
  mutate(Rowsums = Rowsums + 1, logRowsums = log10(Rowsums)) %>% mutate(index_col = as.integer(rownames(rowsums)), logIndex= log10(index_col))
  
#plot the logrowsums
ggplot(rowsums, aes(x = logIndex, y = logRowsums)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Isoform", y = "log10(Rowsums)")
library(ggpubr)
ggline(data=rowsums, x="logIndex",y="logRowsums" ) +
  theme_minimal() +
  theme(legend.position = "none") + #plot index 
  labs(x = "Isoform", y = "log10(Rowsums)") 


