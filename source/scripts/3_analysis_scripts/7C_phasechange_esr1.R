library(ggpubr)
library(ggExtra)
library(ggrepel)
library(dplyr)
library(tidyr)
library(magrittr)

g_one <- "alt-prom-crispr-fiveprime/files/phase_assignment_phase.csv"
#plot the correlations of the guide asignment across both
g_one <- read.csv(g_one,header = TRUE)
g_one$label <- paste0(g_one$target_gene, " ", g_one$phase)
g_one$label_final <- paste0(g_one$target_gene, " ", g_one$phase," ",g_one$guide_num)
g_one_ntc <- g_one %>% filter(target_gene == "non-targeting")

#plot the distribution of ESR1 
p <- g_one %>% filter((target_gene == "ESR1")  ) %>% 
  ggboxplot(x = "promoter_type", y = "percentage", color = "promoter_type",                add = "jitter",
            ,palette = "jco", facet.by = "phase", short.panel.labs = FALSE) 
p  +  stat_compare_means( aes(label = paste0("p = ", after_stat(p.format))))


p <- ggboxplot(ToothGrowth, x = "supp", y = "len",
                color = "supp", palette = "npg",
                add = "jitter",
                facet.by = "dose", short.panel.labs = FALSE)
# Use only p.format as label. Remove method name.
p + stat_compare_means(
  aes(label = paste0("p = ", after_stat(p.format)))
)
# g_one$percentage[g_one$percentage == 0] <- 0.1
##calculate the relative percentage of the target_gene to percentage_nontargeting
percentage_nontargeting <- g_one %>% filter(target_gene  == "ESR1", promoter_type =="MP") %>% 
  group_by(phase) %>% 
  summarise(percentage = mean(percentage)) %>% 
  ungroup() %>% 
  mutate(target_gene = "ESR1")
g_one_esr1 <- g_one %>% filter(target_gene  == "ESR1", promoter_type =="AP")

g_one <- g_one_esr1 %>% 
  left_join(percentage_nontargeting, by = c("phase")) %>% 
  mutate(percentage_change = (percentage.y-percentage.x)/percentage.x) %>% 
  dplyr::select(-percentage.x, -percentage.y) 
##subselect esr1 for gene
my_comparisons <- list(c("G1", "G2M"), c("S", "G1"), c("S", "G2M"))
g_one_esr1 <- g_one %>% filter(gene.x == "ESR1", promoter_type =="AP") 

ggboxplot(g_one_esr1,x = "phase", y = "percentage_change", 
          color = "phase",palette = "jco")+ 
stat_compare_means(comparisons = my_comparisons,method="t.test") 



#plot the percentage in each promoter_type 
plot_esr1 <-g_one %>% filter(gene.x == "ESR1", promoter_type =="MP") %>% select(promoter_type , percentage_change, phase)
my_comparisons <- list( c("G1", "G2M"), c("S", "G1"), c("S", "G2M") )
ggboxplot( plot_esr1, x = "phase", y = "percentage_change",palette =c("#00AFBB", "#E7B800", "#FC4E07"),
           color = "phase")+
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Phase", y = "Percentage Change Relative ") +
  stat_compare_means(comparisons = my_comparisons,method="t.test")
##subselect esr1 for gene
