#read the files/singlecell_shortread_analysis/differential_exp_list_full.csv
library(ggpubr)
library(reshape2)
library(dplyr)
dfm <- read.csv("alt-prom-crispr-fiveprime/files/singlecell_shortread_analysis/differential_exp_list.csv")
#create a dot chart that shows per the P1 and P2 for the number of 
#group the target_gene, promoter MP_Gene_num, AP_gene_num , Overlap_Gene_Num
dfm <- dfm %>% 
  select(Target_Gene, MP_Gene_Num, AP_Gene_Num, Overlap_Gene_Num) %>%
  unique()
#melt the dfm 
dfm <- melt(dfm, id.vars = "Target_Gene")
#set every AP Gene number as negative
dfm$log2value <- log2(dfm$value)

# dfm$value[dfm$variable == "AP_Gene_Num"] <- -dfm$value[dfm$variable == "AP_Gene_Num"]
dfm$log2value[dfm$variable == "AP_Gene_Num"] <- -dfm$log2value[dfm$variable == "AP_Gene_Num"]

#dropna
dfm <- dfm[!is.na(dfm$log2value),]
#create a new column that set the log2value 
dfm$log2value_disp <- ifelse(dfm$variable=="Overlap_Gene_Num",0,dfm$log2value)
#create the another column with only the overlap number 
dfm$log2value_label <- ifelse(dfm$variable=="Overlap_Gene_Num",dfm$value,"")
ggdotchart(dfm, x = "Target_Gene", y = "log2value_disp",
           color = "variable",                                # Color by groups
           palette = c("#00AFBB", "#E7B800", "#6F69AF"), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           rotate = TRUE,                                # Rotate vertically
           dot.size = 5,                                 # Large dot size
           y.text.col = TRUE,
           label = dfm$log2value_label,                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 8, 
                             vjust = 0.5), 
           ggtheme = theme_pubr()                        # ggplot2 theme
)+ scale_y_continuous(breaks = c(-10,-6,-2,0,2,6,10))+
  theme_cleveland()                              # Add dashed grids

#set xlabels 

