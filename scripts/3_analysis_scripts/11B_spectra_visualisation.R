library(ggplot2)
library(dplyr)
library(pheatmap)

#load two files alt-prom-crispr-fiveprime/files/spectra_boxplot.csv
file1 <- read.csv("alt-prom-crispr-fiveprime/files/spectra_boxplot.csv")
g1m_data <- read.csv("~/Desktop/Weatheritt_Lab_Y2/alt-prom-crispr-fiveprime/files/mean_df_melt_10_pivot_sub.csv")
g2_data <- read.csv("~/Desktop/Weatheritt_Lab_Y2/alt-prom-crispr-fiveprime/files/mean_df_melt_70_pivot_sub.csv")

#get mean values per group file1$variable 
file1_mean_variable <- file1 %>% 
  group_by(variable) %>% 
  summarise(mean_value = mean(value)) 
#sory by the highest to lowest mean value
file1_mean_variable <- file1_mean_variable %>% 
  arrange(mean_value)
#set levels of variable to be in order of mean value 
file1$variable <- factor(file1$variable, levels = file1_mean_variable$variable)

file <- file1 %>% 
  ggplot( aes(y=variable, x=value, color=variable)) + 
  geom_boxplot() + 
  theme_minimal() +  #sort y variable by the mean value 
  #remove lines 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "Mean Spectra Score by Promoter", x = "Promoter", y = "Mean Spectra Score") + 
  scale_color_manual(values = c("#FBD0D1", "#F7A5A7", "#801E1E", "#E55352", "#F37374", "#B54242", "#0F4B82", "#5A9BD3", "#81B9E5", "#3E7BBF")) + 
  theme(legend.position = "none")
file


#two heatmaps 
#place the gene column as row names
rownames(g1m_data) <- g1m_data$gene
rownames(g2_data) <- g2_data$gene
#remove gene column
g1m_data <- g1m_data %>% select(-gene)
g2_data <- g2_data %>% select(-gene)
a <- pheatmap(g1m_data, cluster_rows = FALSE, cluster_cols = FALSE,
              show_rownames = TRUE, show_colnames = TRUE , #set size of vox
              cellwidth = 10, cellheight = 10 , #same range for both 
              #create breaks fro the 
              breaks =  seq(min(g1m_data$P1,g1m_data$P2,g2_data$p1, g2_data$p2), max(g1m_data$P1,g1m_data$P2,g2_data$p1, g2_data$p2), length.out = 100),
              #colour rbrewer range RdYlBlu
              color= rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100))) #set the range of the heatmap


b <- pheatmap(g2_data, cluster_rows = FALSE, cluster_cols = FALSE,
              show_rownames = TRUE, show_colnames = TRUE , #set size of vox
              cellwidth = 10, cellheight = 10 , #same range for both 
              #create breaks fro the 
              breaks =  seq(min(g1m_data$P1,g1m_data$P2,g2_data$p1, g2_data$p2), max(g1m_data$P1,g1m_data$P2,g2_data$p1, g2_data$p2), length.out = 100),
              #colour rbrewer range RdYlBlu
              color= rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100))) #set the range of the heatmap



