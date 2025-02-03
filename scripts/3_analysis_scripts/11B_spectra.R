#library read in dataframe_spectra_boxplot
library(ggplot2)
library(reshape2)
library(ggpubr)

#readin file
spectra_boxplot <- read.csv("~/Desktop/Weatheritt_Lab_Y2/alt-prom-crispr-fiveprime/files/spectra_boxplot.csv", row.names=1)
spectra_boxplot <- as.data.frame(spectra_boxplot)
#sort other way
#remove legend
#change colour

#plot using ggpubr boxplot
ggboxplot(spectra_boxplot, x = "variable",
          y = "value", color = "variable", 
           ylab = "Spectra", xlab = "Group", 
          orientation="horizontal",palette =c("#397BC5","#73BDFF","#529CE6",
                                              "#104982","#B44141","#FF7373",
                                              "#E65252","#7F1D1D","#FFA6A6",
                                              "#FFD1D1","#FFA6A6","#FF7373"),
          legend.title="Group", legend.position = "right")


#import other spectra files mean_df_melt_10_pivot_sub
mean_df_melt_10_pivot_sub <- read.csv("~/Desktop/Weatheritt_Lab_Y2/alt-prom-crispr-fiveprime/files/mean_df_melt_10_pivot_sub.csv", row.names=1)
mean_df_melt_10_pivot_sub <- as.data.frame(mean_df_melt_10_pivot_sub)
mean_df_melt_70_pivot_sub <- read.csv("~/Desktop/Weatheritt_Lab_Y2/alt-prom-crispr-fiveprime/files/mean_df_melt_70_pivot_sub.csv", row.names=1)
mean_df_melt_70_pivot_sub <- as.data.frame(mean_df_melt_70_pivot_sub)

mean_df_melt_10_pivot_sub
#form heatmap
library("RColorBrewer")
library(pheatmap)
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
pheatmap(mean_df_melt_10_pivot_sub, cluster_rows = F, cluster_cols = F, color = col, breaks = seq(0.016,0.007, length.out = 256), 
         fontsize = 8, cellwidth = 10, cellheight = 10)
pheatmap(mean_df_melt_70_pivot_sub, cluster_rows = F, cluster_cols = F, color = col, breaks = seq(0.016, 0.007, length.out = 256), 
         fontsize = 8, cellwidth = 10, cellheight = 10)

  
