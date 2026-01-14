
###---------------
#LOAD LIBRARIES
###----------------
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
library(RColorBrewer)
library(RColorBrewer)
library(tidyr)
library(ggalluvial)

###---------------
#Brief 
###----------------
#specifically analyzing transcript-level TPM values and calculating knockdown efficiencies for different promoters (P1 and P2)


###---------------
#DEFINE LOCATION AND DIRECTIORY
###----------------
dir <- 'alt-prom-crispr-fiveprime/files/whippet_dedup'
files <- list.files(dir, pattern = ".gz", full.names = TRUE)
table_name <- file.path(dir, "all_data.txt")


# Print the list of files
print(files)

###---------------
#DEFINE FUNCTION
###----------------
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

#tx2gene has list of transcripts and their associated genes 
tx2gene_loc <- "alt-prom-crispr-fiveprime/files/kallisto/transcriptome/tx2gene.txt"
tx2gene <- read_delim(tx2gene_loc, delim = "\t", col_names = c("transcript_id", "gene_id","gene_name","transcript_name","chromosome","start","end","strand"))
print(sum((tx2gene$start - tx2gene$end) > 0))
tx2gene$transcript_index <- 1:nrow(tx2gene)
tx2gene$tss <- ifelse(tx2gene$strand == "+", tx2gene$start, tx2gene$end)


#read in the cells that are assigned, place into metadata find
#read in the bed file which conrtains p1 and p2 locations
bed_loc <- "alt-prom-crispr-fiveprime/cellranger_input/candidate_AP_MP.bed"
bed <- read_delim(bed_loc, delim = "\t", col_names = c("chromosome", "start", "end", "name", "strand"))
print(sum((bed$start - bed$end) > 0))
bed$tss <- ifelse(bed$strand == "+", bed$start, bed$end)
#do an example for esr1 mp and ap subpopulations
#split _ and get last 
bed$promoter_type <- str_split(bed$name, "_") %>% map_chr(~last(.x))
#get the second to last element
bed$gene_name <- str_split(bed$name, "_") %>% map_chr(~.x[length(.x)-1])
bed$promoter_gene <- paste0(bed$gene_name,"_",bed$promoter_type)
#check the number of gene name prior to merge
print(length(unique(bed$gene_name)))
print(length(unique(bed$promoter_gene)))
#merge the dataframe on gene name 

#merge the bed and tx2gene on gene name
bed <- merge(bed, tx2gene, by = "gene_name", suffixes = c("_target","_assigned"), how="inner")
print(length(unique(bed$gene_name)))
print(length(unique(bed$promoter_gene)))
#filter for for rows that only have 
#calculate column that has the distance from the tss
bed$distance <- abs(bed$tss_target - bed$tss_assigned)
#filter for rows that have a distance of less than 1000
bed <- bed %>% group_by(transcript_id) %>% mutate(min_distance = distance == min(distance))
print(length(unique(bed$transcript_id))==length(bed$transcript_id[bed$min_distance == TRUE]))
ggplot(bed) + 
  # geom_histogram(aes(x=distance, fill=min_distance),bins = 100) + theme_minimal() +
  labs(x = "Distance from TSS", y = "Count") + #add denistty plot
  geom_density(aes(x=distance, color=min_distance), alpha = 0.5)

#then find the numbers that assign to those transcripts
#get the number from the matrix number for those 
bed <- bed %>% filter(min_distance == TRUE)
#only keep genes with both P1 and P2
all_data_names <- stringr::str_split_i(colnames(all_data), "Tpm_", 2) 
#use these transcripts 
bed <- bed %>% group_by(gene_name) %>% filter(promoter_gene %in% unique(all_data_names)) %>% filter(n_distinct(promoter_type) == 2)
bed <- bed %>% filter(distance < 1000)
#split transcript ID by . and get the first element name column ISoform
bed$Isoform <- str_split_fixed(bed$transcript_id, "\\.", 2)[,1]
#merge bed and all_data
###!!!Load the conversion between transcript to gene !!!####
#read in the position --> transcript number conversion 

#merge the bed and all_data on Isoform
all_data <- merge(all_data, bed, by = "Isoform")
head(all_data)

#create empty table gene_table with columsn for both mp and ap
gene_table_ontarget <- data.frame(gene = character(), percent_knockdown_mp = numeric(), distance_mp = numeric(), Isoform_mp = character(), transcript_name_mp = character(), neg_ctrl_tpm_mp= character(), percent_knockdown_ap = numeric(), distance_ap = numeric(), Isoform_ap = character(), transcript_name_ap = character(), neg_ctrl_tpm_ap= character(), stringsAsFactors = FALSE)
gene_table_offtarget <- data.frame(gene = character(), percent_knockdown_mp = numeric(), distance_mp = numeric(), Isoform_mp = character(), transcript_name_mp = character(),  neg_ctrl_tpm_mp= character(),percent_knockdown_ap = numeric(), distance_ap = numeric(), Isoform_ap = character(), transcript_name_ap = character(),  neg_ctrl_tpm_ap= character(),stringsAsFactors = FALSE)
gene_table_total <- data.frame(gene = character(), percent_knockdown_mp_total = numeric(), percent_knockdown_ap = numeric(), stringsAsFactors = FALSE)
#loop through every gene 
#unique gene names
for (gene in unique(all_data$gene_name)){
  print(gene)
  gene_prom_mp <- paste0(gene, "_MP")
  gene_prom_mp_tpm <- paste0("Tpm_", gene_prom_mp)
  gene_prom_ap <- paste0(gene, "_AP")
  gene_prom_ap_tpm <- paste0("Tpm_", gene_prom_ap)
  negative_column <- "Tpm_Negative"
  #subset the data for the gene
  #firrst for main promoter
  gene_data_mp <- all_data %>% filter(promoter_gene == gene_prom_mp) 
  #then for the alternative promoter
  gene_data_ap <- all_data %>% filter(promoter_gene == gene_prom_ap)
  #check if all the sum(gene_data_mp[[gene_prom_mp_tpm]]) or sum(gene_data_mp[[gene_prom_ap_tpm]]) for every column is greater than the minimum filter tpm
  # if (( > filter) & (sum(gene_data_ap[[negative_column]]) > filter) ){
  #   print("skipped")
  #   next
  # }
  print(sum(gene_data_mp[[negative_column]]))
  print(sum(gene_data_ap[[negative_column]]))
  #calculate the % knockdown between negative control and ESR1_MP population by summing the Tpm and then dividing by the negative control
  #fill NA wuth 0 in gene_data_ap[[gene_prom_ap_tpm]] and gene_data_mp[[gene_prom_mp_tpm]] and negative control and gene_data_mp[[gene_prom_ap_tpm]]
  total_mp <- sum(gene_data_mp[[gene_prom_mp_tpm]]) + sum(gene_data_mp[[gene_prom_ap_tpm]])
  total_ap <- sum(gene_data_ap[[gene_prom_mp_tpm]]) + sum(gene_data_ap[[gene_prom_ap_tpm]])
  total_percentage_knockdown_mp <- (total_mp - sum(gene_data_mp[[negative_column]]))/sum(gene_data_mp[[negative_column]])*100
  total_percentage_knockdown_ap <- (total_ap - sum(gene_data_ap[[negative_column]]))/sum(gene_data_ap[[negative_column]])*100
  gene_data_mp_percent_knockdown_mp <- ( sum(gene_data_mp[[gene_prom_mp_tpm]])-sum(gene_data_mp[[negative_column]]))/sum(gene_data_mp[[negative_column]])*100
  gene_data_ap_percent_knockdown_ap <- ( sum(gene_data_ap[[gene_prom_ap_tpm]])-sum(gene_data_ap[[negative_column]]))/sum(gene_data_ap[[negative_column]])*100
  gene_data_mp_percent_knockdown_ap <- ( sum(gene_data_mp[[gene_prom_ap_tpm]])-sum(gene_data_mp[[negative_column]]))/sum(gene_data_mp[[negative_column]])*100
  gene_data_ap_percent_knockdown_mp <- ( sum(gene_data_ap[[gene_prom_mp_tpm]])-sum(gene_data_ap[[negative_column]]))/sum(gene_data_ap[[negative_column]])*100
  # whenever sum(gene_data_mp[[gene_prom_mp_tpm]]) == sum(gene_data_mp[[negative_column]])) make the value equal 0
  gene_data_mp_percent_knockdown_mp[sum(gene_data_mp[[gene_prom_mp_tpm]]) == sum(gene_data_mp[[negative_column]])] <- 0
  gene_data_ap_percent_knockdown_ap[sum(gene_data_ap[[gene_prom_ap_tpm]]) == sum(gene_data_ap[[negative_column]])] <- 0
  gene_data_mp_percent_knockdown_ap[sum(gene_data_mp[[gene_prom_ap_tpm]]) == sum(gene_data_mp[[negative_column]])] <- 0
  gene_data_ap_percent_knockdown_mp[sum(gene_data_ap[[gene_prom_mp_tpm]]) == sum(gene_data_ap[[negative_column]])] <- 0
  
 
  gene_list_kd <- list(gene,total_percentage_knockdown_mp,total_percentage_knockdown_ap)
  gene_list_ontarget <- list(gene, sum(gene_data_mp[[gene_prom_mp_tpm]]),
                             gene_data_mp_percent_knockdown_mp,median(gene_data_mp$distance),str_flatten_comma(gene_data_mp$Isoform),str_flatten_comma(gene_data_mp$transcript_name),sum(gene_data_mp[[negative_column]]), 
                             sum(gene_data_ap[[gene_prom_ap_tpm]]),
                             gene_data_ap_percent_knockdown_ap,median(gene_data_ap$distance),str_flatten_comma(gene_data_ap$Isoform),str_flatten_comma(gene_data_ap$transcript_name),sum(gene_data_ap[[negative_column]]))
  gene_list_offtarget <- list(gene, sum(gene_data_mp[[gene_prom_ap_tpm]]),
                              gene_data_mp_percent_knockdown_ap,median(gene_data_mp$distance),str_flatten_comma(gene_data_mp$Isoform),str_flatten_comma(gene_data_mp$transcript_name),sum(gene_data_mp[[negative_column]]), 
                              sum(gene_data_ap[[gene_prom_mp_tpm]]),
                              gene_data_ap_percent_knockdown_mp,median(gene_data_ap$distance),str_flatten_comma(gene_data_ap$Isoform),str_flatten_comma(gene_data_ap$transcript_name),sum(gene_data_ap[[negative_column]]))
  #add list to table keep column names
  print(gene_prom_ap_tpm)
  gene_table_ontarget <- rbind(gene_table_ontarget, gene_list_ontarget)
  gene_table_offtarget <- rbind(gene_table_offtarget, gene_list_offtarget)
  gene_table_total <- rbind(gene_table_total, gene_list_kd)
}

#change gene_table names 
column_names <- c("gene", "sum_MP","percent_knockdown_MP", "distance_MP", "Isoform_MP", "transcript_name_MP", "neg_ctrl_tpm_MP",
                  "sum_AP","percent_knockdown_AP", "distance_AP", "Isoform_AP", "transcript_name_AP","neg_ctrl_tpm_AP")
colnames(gene_table_ontarget) <-column_names
colnames(gene_table_offtarget) <- column_names
colnames(gene_table_total) <- c("gene", "percent_knockdown_MP_total", "percent_knockdown_AP_total")
gene_table <- merge(gene_table_ontarget, gene_table_offtarget, by = "gene", suffixes = c("_ontarget", "_offtarget"))
gene_table <- merge(gene_table, gene_table_total, by = "gene")
head(gene_table)
#plot the data
gene_table <- as.data.frame(gene_table)
dim(gene_table)
gene_table$successfulKD <- ifelse(gene_table$percent_knockdown_AP_ontarget < 40,TRUE,FALSE)

#save gene_table as a file
write.csv(gene_table,paste0(dir,"/gene_table.csv"))

gene_table_heatmap <-gene_table %>%
 dplyr::select(gene,percent_knockdown_MP_ontarget, percent_knockdown_AP_ontarget, percent_knockdown_MP_offtarget, percent_knockdown_AP_offtarget,percent_knockdown_MP_total,percent_knockdown_AP_total) 
#put column to rowname
rownames(gene_table_heatmap) <- gene_table_heatmap$gene
gene_table_heatmap <- gene_table_heatmap %>% dplyr::select(-gene)
breaks2 <- seq(-100, 100, length.out = 150)


pheatmap(gene_table_heatmap %>% 
           dplyr::select(percent_knockdown_MP_ontarget,percent_knockdown_MP_offtarget,percent_knockdown_MP_total)%>% drop_na() ,
         cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = TRUE, show_colnames = TRUE, main = "P1", 
         fontsize = 8, cellwidth = 8, cellheight = 8 , scale ="none" ,cutree_rows=5 ,colorRampPalette(brewer.pal(9,"Spectral"))(150),
         breaks = breaks2)

breaks2 <- seq(-100, 100, length.out = 150)
#repeat to P2
pheatmap(gene_table_heatmap %>% 
           dplyr::select(percent_knockdown_AP_ontarget,percent_knockdown_AP_offtarget,percent_knockdown_AP_total) %>% drop_na(),
         cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = TRUE, show_colnames = TRUE, main = "P2", 
         fontsize = 8, cellwidth = 8, cellheight = 10 , scale ="none" ,cutree_rows=5 ,colorRampPalette(brewer.pal(9,"Spectral"))(150),
         breaks = breaks2)

gene_table <- gene_table %>%
  drop_na() %>% filter(neg_ctrl_tpm_MP_ontarget > 0.5 & neg_ctrl_tpm_AP_ontarget > 0.5)  

dim(gene_table)


gene_table_melt <-  gene_table %>% 
  filter(percent_knockdown_AP_ontarget < 40) %>% 
  melt(id.vars = c("gene","neg_ctrl_tpm_AP_ontarget","neg_ctrl_tpm_MP_ontarget"),
                                       measure.vars = c("percent_knockdown_MP_ontarget", "percent_knockdown_AP_ontarget", 
                                                        "percent_knockdown_MP_offtarget", "percent_knockdown_AP_offtarget"))
gene_table_melt$target <- ifelse(grepl("ontarget", gene_table_melt$variable), "ontarget", "offtarget")
gene_table_melt$promoter <- ifelse(grepl("MP", gene_table_melt$variable), "P1", "P2")
ggplot(gene_table_melt) +
  geom_boxplot(aes(x=variable, y=value, fill=variable)) + 
  theme_minimal() + labs(x = "Promoter", y = "Percent Knockdown") + #rotate 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #order MP_ontarget MP_offtarget just filter for 20 to -100
  scale_y_continuous(limits = c(-100, 20)) 


ggpaired(gene_table_melt %>% filter(promoter=="P1"), x = "target", y = "value",
         color = "target", line.color = "gray", line.size = 0.4,
         palette = "npg")+
  stat_compare_means(paired = TRUE) + 
  ggtitle("P1") + theme_minimal() + labs(y = "Percent Knockdown") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggpaired(gene_table_melt %>% filter(promoter=="P2"), x = "target", y = "value",
         color = "target", line.color = "gray", line.size = 0.4,
         palette = "npg")+
  stat_compare_means(paired = TRUE) + 
  ggtitle("P2") + theme_minimal() + labs(y = "Percent Knockdown") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



palletee1 <- c("#00AFBB", "#0052bb")
palletee <- c ("#FC4E07", "#fc8b07")

ggboxplot(gene_table_melt, x = "target", y = "value", color = "target", add = "jitter", facet.by = "promoter", ylab = "Percent Knockdown", xlab = "Promoter", 
          palette = c(palletee[1],palletee1[1])) +
  stat_compare_means( paired = FALSE) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


gene_table_melt %>% filter(promoter=="P1") %>%
  ggviolin( x = "target", y = "value",fill = "target",add = "boxplot",    ylab = "Percent Knockdown", xlab = "Promoter", 
           palette = c(palletee1[1],palletee[1]),add.params = list(fill = "white") ) +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means( paired = FALSE) 
gene_table_melt %>% filter(promoter=="P2") %>%
  ggviolin( x = "target", y = "value",fill = "target",add = "boxplot",    ylab = "Percent Knockdown", xlab = "Promoter", 
            palette = c(palletee[1],palletee1[1]) ,add.params = list(fill = "white") ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means( paired = FALSE) 



p<-ggscatter(gene_table_melt %>% filter(target=="ontarget"), x = "neg_ctrl_tpm_MP_ontarget", y = "value",
         color = "promoter", line.color = "gray", line.size = 0.4,
         palette = "npg")+
  ggtitle("P1") + theme_minimal() + labs(y = "Percent Knockdown") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
a <- ggscatter(gene_table_melt %>% filter(target=="offtarget"), x = "neg_ctrl_tpm_AP_ontarget", y = "value",
          color = "promoter", line.color = "gray", line.size = 0.4,
          palette = "npg")+
  ggtitle("P2") + theme_minimal() + labs(y = "Percent Knockdown") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

a+ geom_text(data = gene_table_melt %>% filter(promoter=="P1"), aes(label = gene), nudge_y = 5, check_overlap = TRUE) +
  geom_text(data = gene_table_melt %>% filter(promoter=="P2"), aes(label = gene), nudge_y = 5, check_overlap = TRUE) +
  theme_minimal() + labs(y = "Percent Knockdown") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(limits = c(-100, 20)) + scale_fill_manual(values = palletee)
a + geom_text(data = gene_table_melt %>% filter(promoter=="P1"), aes(label = gene), nudge_y = 5, check_overlap = TRUE) +
  geom_text(data = gene_table_melt %>% filter(promoter=="P2"), aes(label = gene), nudge_y = 5, check_overlap = TRUE) +
  theme_minimal() + labs(y = "Percent Knockdown") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(limits = c(-100, 20)) + scale_fill_manual(values = palletee)
esr1 <- gene_table_melt %>% filter(gene == "ESR1") %>% 
  filter(variable %in% c("percent_knockdown_MP_ontarget","percent_knockdown_MP_offtarget",
                         "percent_knockdown_AP_ontarget","percent_knockdown_AP_offtarget"))


gene_table_melt %>% #filter for value between -100 and 100
  filter(value >= -100 & value <= 100) %>% #filter for tpm
  filter((neg_ctrl_tpm_AP_ontarget > 1) & (neg_ctrl_tpm_MP_ontarget> 1)) %>% #plot density plot of negative ctrl 
  melt(id.vars = c("gene"),measure.vars=c("neg_ctrl_tpm_AP_ontarget","neg_ctrl_tpm_MP_ontarget")) %>%
    unique() %>%
  ggplot(aes(x=value, fill=variable,  alpha=0.4)) + #log scale
  geom_density() +
  labs(x = "Negative Control TPM", y = "Density") + theme_minimal() +
    scale_x_continuous(trans='log10') + scale_fill_manual(values = c("#00AFBB", "#E7B800"))


gene_table_melt <- gene_table_melt %>% drop_na()
gene_table_melt %>% filter(variable %in% c("percent_knockdown_MP_ontarget","percent_knockdown_MP_offtarget")) %>%
  ggplot(aes(x=variable, y=value, group=gene,colour =gene )) +
  geom_line() 


gene_table %>% filter(variable %in% c("neg_ctrl_tpm_MP","neg_ctrl_tpm_AP")) %>%
  ggplot(aes(x=variable, y=value, group=gene,colour =gene )) +
  geom_line()

t.test(gene_table$sum_MP_ontarget, gene_table$neg_ctrl_tpm_MP_ontarget, paired = TRUE)
t.test(gene_table$sum_AP_ontarget, gene_table$neg_ctrl_tpm_AP_ontarget, paired = TRUE)

gene_table_melt %>% filter(target=="ontarget") %>%
   ggviolin( x = "promoter", y = "value",fill = "promoter",add = "boxplot",    ylab = "Percent Knockdown", xlab = "Promoter", 
                    palette = c(palletee[1],palletee1[1]),add.params = list(fill = "white") ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

