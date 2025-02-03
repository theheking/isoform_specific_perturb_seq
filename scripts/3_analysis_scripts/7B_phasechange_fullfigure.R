library(ggpubr)
library(ggExtra)
library(ggrepel)
library(dplyr)
library(tidyr)
library(magrittr)

phase_assignment_phase <- "alt-prom-crispr-fiveprime/files/phase_assignment_phase.csv"
phase_assignment <- read.csv(phase_assignment_phase,header = TRUE)
##show the correlation across the 
g_one_corr <- phase_assignment %>% 
  dplyr::select(target_gene,promoter_type,corr) %>% #drop duplicate rows
  dplyr::distinct() %>% #remove non-targeting
  dplyr::filter(target_gene  != "non-targeting") %>%
  drop_na() 
#plot the average density 
ggdensity(g_one_corr,x = "corr", fill = "promoter_type",add = "mean") +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(title = "Correlation of guide assignment across both promoters",
       x = "Correlation",
       y = "Density") 


phase_assignment$label <- paste0(phase_assignment$target_gene, " ", phase_assignment$phase)
phase_assignment$label_final <- paste0(phase_assignment$target_gene, " ", phase_assignment$phase," ",phase_assignment$guide_num)
phase_assignment <- phase_assignment %>% filter(pvalue != 0)
phase_assignment

#plot the tstat 
g_one_tstat <- phase_assignment %>% 
  dplyr::select(target_gene,promoter_type,tstat) %>% #drop duplicate rows
  dplyr::distinct() %>% #remove non-targeting
  dplyr::filter(target_gene  != "non-targeting") %>%
  drop_na() %>% 
  ggdensity(x = "tstat", fill = "promoter_type",add = "mean")
g_one_tstat

g_one_pivot_nt <- phase_assignment %>% filter(target_gene == "non-targeting") %>%
  group_by(phase) %>%
  mutate(mean_percentage = mean(percentage)) %>%
  select(phase,mean_percentage) %>% #dropduplicates
  distinct() %>% ungroup() 


#for every gene and guide_num and phase get the mean % per promoter
# g_one_pivot <- phase_assignment %>% filter(gene != "non-targeting") %>%
#   left_join(g_one_pivot_nt, by = "phase") %>%
#   dplyr::select(gene,promoter,phase,percentage,pvalue) %>%
#   pivot_wider(names_from = promoter, values_from = c("percentage","pvalue"),values_fn = mean) %>% #make a new column with minimum
#   mutate(pval_min= pmin(pvalue_AP,pvalue_MP)) %>%
#   mutate(logneg_pval= -log10(pval_min), significant =pval_min<0.05 )  

g_one_pivot <- phase_assignment %>% filter(target_gene != "non-targeting") %>%
  left_join(g_one_pivot_nt, by = "phase") %>%
  dplyr::select(target_gene,promoter_type,phase,percentage) %>%
  pivot_wider(names_from = promoter_type, values_from = c("percentage"), values_fn = mean) 
# %>%
  # mutate(fdry_n= ifelse((fdr_num_MP > 0.5)|(fdr_num_AP > 0.5),TRUE,FALSE))

phase_assignment_minpvalue <- phase_assignment %>% select(target_gene,phase,adj_pvalue) %>%
  group_by(target_gene,phase) %>% #extract minimum
  mutate(min_pvalue = min(adj_pvalue)) %>% 
  select(target_gene,phase,min_pvalue) %>%
  unique() %>%
  mutate(significant = min_pvalue < 0.05) %>% 
  ungroup() 

#merge phase_assignment_minpvalue g_one_pivot
g_one_pivot <- g_one_pivot %>% left_join(phase_assignment_minpvalue,
                                         by = c("target_gene","phase")) 

color1="#BD4526"
color2="#7671B3"
color3="#347238"
#create palette
color_palette <- c(color1,color2,color3)
p <- ggplot() +
  geom_point(data=g_one_pivot,aes(y = AP, x = MP, color = phase, size=significant) )+ 
  # geom_line(data=as.data.frame(g_one_pivot_nt),aes(x=mean_percentage,y=mean_percentage,color='black',alpha=1))  + 
  scale_color_manual(values = color_palette)+
  scale_size_manual(values = c(1.2, 2.8))+ 
  labs(title = "Phase assignment of guides",
       x = "Mean percentage MP",
       y = "Mean percentage AP") +
  theme(legend.position = "none")


p <-   p+ geom_text_repel(
  data =g_one_pivot[(g_one_pivot$significant==TRUE),],
  aes(y = AP, x = MP,label = target_gene), max.overlaps =getOption("ggrepel.max.overlaps", default = 80) ,
  size = 3,min.segment.length = unit(0, 'lines'), 
) + theme_minimal()  + #place range
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + #range xlim ylim
  xlim(20, 50) + ylim(20, 50) 
ggMarginal(p, type = "densigram", groupColour = TRUE, groupFill = TRUE) 



