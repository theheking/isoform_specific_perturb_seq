library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggExtra)
library(ggrepel)
library(tidyr)
library(stringr)

palette <- c("#50acb9","#de6239","#deb947")
cnv_score="alt-prom-crispr-fiveprime/files/singlecell_shortread_analysis/cnv_score.csv"
cnv_score=read.csv(cnv_score)
#plot 
cnv_score=as.data.frame(cnv_score %>% filter(binding_regions %in% c("esr1") & (elements=="T")))
dim(cnv_score)
#count the number of promoter_type column equivlauent of value_counts
cnv_score %>% count(promoter_type)


#vilin plot wilcox.test()	kruskal.test()	
ggviolin(cnv_score, x = "promoter_type", y = "cnv_score", fill = "promoter_type",  adjust=2,
         add = "boxplot",add.params = list(fill = "white")) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs( x = "Cell cnv_score", y = "CNV score") +
  theme(legend.position = "none") + #choose a a palette 
  scale_fill_manual(values = palette) + #add statistics
  stat_compare_means( label.y = c(0.8,0.75), label.x = 1, method = "wilcox.test", comparisons = list(c("Control","AP"),c("Control","MP"))   , 
                      hide.ns = FALSE, size = 3, tip.length = 0.01, method.args = list(alternative = "two.sided")) +
  annotate("text", x=1, y=0.8, label= "n=897") 



#order cnv to be Control MP AP
cnv_score$promoter_type <- factor(cnv_score$promoter_type, levels = c("Control", "MP", "AP"))
#reoder the levels
# cnv_score$promoter_type <- factor(cnv_score$promoter_type, levels = c("AP", "Control","MP"))
##conver this to ggplotgeom_vuiolin
ggplot(cnv_score, aes(x=promoter_type, y=cnv_score, fill=promoter_type)) +
  geom_violin(adjust=7) +
  geom_boxplot(width=0.1) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs( x = "Cell cnv_score", y = "CNV score") +
  theme(legend.position = "none") + #choose a a palette 
  scale_fill_manual(values = palette)  +
  geom_signif(  test = "wilcox.test",y_position = c(0.715, 0.73,0.748),
 comparisons = list(c("Control", "MP"),c("Control", "AP"),c("MP","AP")))
  
##import cnv_score_list_total
cnv_score="alt-prom-crispr-fiveprime/files/singlecell_shortread_analysis/cnv_score_allgene.csv"
cnv_score_all=read.csv(cnv_score)
dim(cnv_score_all)
cnv_score_all <- unique(cnv_score_all)
dim(cnv_score_all)
cnv_score_neg <- cnv_score_all %>% filter(promoter_type=="Control")
#describe cnv_score_neg


#select only for promoter_type MP and AP
cnv_score_all <- cnv_score_all %>% filter(promoter_type!="Control")
#create pivot table for mean of MP A
cnv_score_all_pivot <- cnv_score_all %>% 
  dplyr::select(promoter_type, cnv_score,perturbation) %>%
  pivot_wider(names_from = promoter_type, values_from = cnv_score,values_fn = mean) %>% #make a new column with minimum
  mutate(pval_min= pmin(AP,MP)) %>%
  mutate(logneg_pval= -log10(pval_min), significant =pval_min<0.05 ) #get the minimum pvalue


cnv_score_all_pivot <- drop_na(cnv_score_all_pivot)
#rread in the nterm 
nterm <- "alt-prom-crispr-fiveprime/files/reference/all_pivot_simple_nterm.txt"
nterm <- read.table(nterm, header=TRUE)
#merge nterm with cnv_score_all_pivot
cnv_score_all_pivot <- merge(cnv_score_all_pivot,nterm, by.x="perturbation", by.y="Gene_symbol")
#read in the significant distance edistance
#plot the mean of the cnv_score for MP and AP
p <- ggplot() +
  geom_point(data=cnv_score_all_pivot,aes(y = AP, x = MP, color=Nterminus_Change), alpha=0.6 )+ 
  scale_color_brewer(palette = "Dark2")+
  scale_size_manual(values = c(1.2, 2.8))
  
model <- lm(AP ~ MP, data = cnv_score_all_pivot)

cnv_score_all_pivot <- cnv_score_all_pivot %>%
  mutate(residuals = resid(model))

residual_threshold <- 2 * sd(cnv_score_all_pivot$residuals)

cnv_score_all_pivot_resid <- cnv_score_all_pivot %>%
  mutate(outlier = ifelse(abs(residuals) > residual_threshold, TRUE, FALSE))

#create the cnv_score_neg onto a line
ggMarginal(p,type = "densigram", groupColour = TRUE, groupFill = TRUE)

p <- p+
  geom_line(data=cnv_score_neg , aes(x=cnv_score ,y=cnv_score) , linetype=1 ) +
   geom_text_repel(data =unique(cnv_score_all_pivot_resid[cnv_score_all_pivot_resid$outlier==TRUE,]),
    aes(y = AP, x = MP,label = perturbation), max.overlaps =getOption("ggrepel.max.overlaps", default = 80) ,
    size = 3,min.segment.length = unit(0, 'lines'), 
  ) + theme_minimal() + scale_color_brewer(palette = "Dark2") + 
  #plot line across x=y
   xlim(0.255,0.3) + ylim(0.255,0.3) +
  labs(x = "AP CNV score", y = "MP CNV score")
p
  # geom_abline(intercept = 0, slope = 1, linetype = "dashed") + xlim(0.1225,0.131) + ylim(0.1225,0.131) +
ggMarginal(p,type = "densigram", groupColour = TRUE, groupFill = TRUE)



##import individual gene cnv score 
cnv_score_allgene="alt-prom-crispr-fiveprime/files/cnv_mean_score_esr1.csv"
cnv_score_allgene=read.csv(cnv_score_allgene, row.names = 1)
###
library(dplyr)
#import tsv and then get the chr6 amp
head(cnv_score_allgene)
#mutate cnv_score_allgene and get largest variation
cnv_score_allgene <- cnv_score_allgene %>% 
  mutate(variance = apply(cnv_score_allgene[c("ESR1_AP",   "ESR1_MP","non.targeting_Control")], 1, var))  %>% 
  arrange(desc(variance))

cnv_score_allgene_6 <- cnv_score_allgene %>% 
  filter(chromosome=="chr6") %>% #sortby
  filter(str_detect(gene_id, "HIST")) %>%
  arrange(start) %>% #filter after 26031594
  filter(start > 26031594) %>%
  head(n=4)
    # head(n=80) %>% tail (n =35) %>% head(n=8) %>%
  #filter for genes that begin with HIST
# %>%
  # tail(n=4)

cnv_score_allgene_filtered <- cnv_score_allgene %>% filter(chromosome=="chr10"|chromosome=="chr11") %>% #sortby
  arrange(start) %>% filter(gene_id %in% c("CDK1","CCND1"))

#append cnv_score_allgene_filtered and cnv_score_allgene_6
cnv_score_allgene_6 <- rbind(cnv_score_allgene_6,cnv_score_allgene_filtered)
#order to chromosome 6 and then start 
# cnv_score_allgene_6 <- cnv_score_allgene_6 %>% filter(str_detect(gene_id, "HIST"))
rownames(cnv_score_allgene_6) <- cnv_score_allgene_6$gene_id

#filter for genes that begin with HIST
#plot the three values ESR1_AP      ESR1_MP non.targeting_Control
#and gene_id
library(RColorBrewer)
max(cnv_score_allgene_6[c("ESR1_AP","ESR1_MP","non.targeting_Control")])
min(cnv_score_allgene_6[c("ESR1_AP","ESR1_MP","non.targeting_Control")])
#center zero around white
paletteLength <- 20
test <- cnv_score_allgene_6[c("ESR1_AP","ESR1_MP","non.targeting_Control")]

myColor <- rev(colorRampPalette(brewer.pal(10, "Spectral"))(paletteLength))
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))

pheatmap::pheatmap(test, cluster_rows=FALSE,
                   cluster_cols=TRUE, show_rownames=TRUE, show_colnames=TRUE, fontsize=8, border_color = NA, 
                   fontsize_row = 8, fontsize_col = 8, cellwidth = 10, cellheight = 10, scale = "row")
###cnv_score_allgene_6 CDK1 and CCND1


rownames(cnv_score_allgene_filtered) <- cnv_score_allgene_filtered$gene_id
#filter for genes that begin with HIST
#plot the three values ESR1_AP      ESR1_MP non.targeting_Control
#and gene_id
library(RColorBrewer)
max(cnv_score_allgene_filtered[c("ESR1_AP","ESR1_MP","non.targeting_Control")])
min(cnv_score_allgene_filtered[c("ESR1_AP","ESR1_MP","non.targeting_Control")])
#center zero around white
paletteLength <- 20
test <- cnv_score_allgene_filtered[c("ESR1_AP","ESR1_MP","non.targeting_Control")]

myColor <- rev(colorRampPalette(brewer.pal(10, "Spectral"))(paletteLength))
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))


pheatmap::pheatmap(test, cluster_rows=FALSE,
                   cluster_cols=FALSE, show_rownames=TRUE, show_colnames=TRUE, fontsize=8, border_color = NA, 
                   fontsize_row = 8, fontsize_col = 8, cellwidth = 10, cellheight = 10, scale = "row")
###cnv_score_allgene_6 CDK1 and CCND1

pheatmap::pheatmap(test, cluster_rows=FALSE)



