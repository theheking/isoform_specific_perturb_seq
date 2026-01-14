

library(shiny)
library(stringr)
library(dplyr)
library(pcr)  ## For qPCR analyis
library(ggplot2)  ##For graphs
library(readxl)
library(rlang)
library(tidyr)
library(DT)


#read in the final biological replicate and then calculate the delta delta ct
excel_filename="alt-prom-crispr-fiveprime/files/esr1/2024-07-12ESR1KO.xls"

zscore_thresh=2
reference_group="NTC"
sample_column_name="sample_id"
target_name="Target Name"
ct_column_name="CT" 
return_data="normal"
mode_tube="separate_tube"
method_choice="delta_delta_ct" 
reference_gene="hGAPDH"
test_reference_group="NTC"


radiostat="t.test"
test_gene="hGAPDH"
refgroup="NTC"
test_group="ESR1 KO MP"


results_readin<- readin_quantstudio(excel_filename=excel_filename)
#trim white space in Sample Name 
results_readin$`Sample Name` <- trimws(results_readin$`Sample Name`)
#create two new columns one that contains everything before Rep and one that contains Rep
results_readin$sample_id <- trimws(stringr::str_extract(results_readin$`Sample Name`, ".*(?=Rep)"))
#extract only the final rep number in a column
results_readin$`replicate_num` <- ifelse(grepl("Rep II" ,results_readin$`Sample Name`), "Rep II", "Rep I")

#split into two different dataframes, run the reference on one and test the other 
results_readin_rep2 <- results_readin %>% filter(`replicate_num` == "Rep II")
#get the CT value with column sample_id
#for each group at number 
results_readin_rep2 <- sort_by(results_readin_rep2, results_readin_rep2$`Target Name`)
results_readin_rep2$number <- rep(c(1:18),nrow(results_readin_rep2)/18)
print(results_readin_rep2$number)
results_readin_rep2_simple <- as.data.frame(results_readin_rep2 %>% 
                                                select("number",sample_column_name,target_name,ct_column_name) %>%
  pivot_wider(names_from=target_name,values_from=ct_column_name,values_fn = list ) %>% select(-c("number"))) %>% 
  select(c("sample_id","ESR1_Common","ESR1_P1",'hGAPDH')) %>%
  filter(sample_id  %notin% c("Wt dcas9") )
group_var <- ifelse(results_readin_rep2_simple$sample_id=="ESR1 KO MP1","ESR1 KO MP",results_readin_rep2_simple$sample_id)
group_var <- ifelse(group_var=="ESR1 KO MP2","ESR1 KO MP",group_var)
group_var <- ifelse(group_var=="ESR1 KO AP1","ESR1 KO AP",group_var)
group_var_2 <- ifelse(group_var=="ESR1 KO AP2","ESR1 KO AP",group_var)

results_readin_rep2_simple <- results_readin_rep2_simple %>% select(-c(sample_id))
results_summary_2 <- mutate_all(results_readin_rep2_simple, function(x) as.numeric(as.character(x)))
# View(results_summary)
#delta delta CT 


results_summary_new_p2 <- pcr_ddct(
  df = results_summary_2, ## Name of the table loaded in line 13
  group_var = group_var_2, ## Select the groups, generated in line 17
  reference_gene = reference_gene, ## Define the reference gene of your table
  reference_group = reference_group, ## Define the control group, from line 17
  plot=TRUE)  ##to indicate that GAPDH and PTBP1 were measured in the same tube
results_summary_new



results_readin_rep1 <- results_readin %>% filter(`replicate_num` == "Rep I")
#get the CT value with column sample_id
#for each group at number 
results_readin_rep1 <- sort_by(results_readin_rep1, results_readin_rep1$`Target Name`)
results_readin_rep1$number <- rep(c(1:18),nrow(results_readin_rep1)/18)
print(results_readin_rep1$number)
results_readin_rep1_simple <- as.data.frame(results_readin_rep1 %>% 
                                              select("number",sample_column_name,target_name,ct_column_name) %>%
                                              pivot_wider(names_from=target_name,values_from=ct_column_name,values_fn = list ) %>% select(-c("number"))) %>% 
  select(c("sample_id","ESR1_Common","ESR1_P1",'hGAPDH')) %>%
  filter(sample_id  %notin% c("Wt dcas9") )
group_var <- ifelse(results_readin_rep1_simple$sample_id=="ESR1 KO MP1","ESR1 KO MP",results_readin_rep1_simple$sample_id)
group_var <- ifelse(group_var=="ESR1 KO MP2","ESR1 KO MP",group_var)
group_var <- ifelse(group_var=="ESR1 KO AP1","ESR1 KO AP",group_var)
group_var_1 <- ifelse(group_var=="ESR1 KO AP2","ESR1 KO AP",group_var)

results_readin_rep1_simple <- results_readin_rep1_simple %>% select(-c(sample_id))
results_summary_1 <- mutate_all(results_readin_rep1_simple, function(x) as.numeric(as.character(x)))
# View(results_summary)
#delta delta CT 


results_summary_new <- pcr_ddct(
  df = results_summary_1, ## Name of the table loaded in line 13
  group_var = group_var_1, ## Select the groups, generated in line 17
  reference_gene = reference_gene, ## Define the reference gene of your table
  reference_group = reference_group, ## Define the control group, from line 17
  plot=TRUE)  ##to indicate that GAPDH and PTBP1 were measured in the same tube
results_summary_new

#concatenate the two p1 and p2 
one <- results_summary_new[[1]]
one$sample <- "Rep I"
two <- results_summary_new_p2[[1]]
two$sample <- "Rep II"
results_summary_all <- rbind(one,two)

#just do it only 
results_summary_all <- pcr_ddct(
  df = rbind(results_summary_1,results_summary_2), ## Name of the table loaded in line 13
  group_var = c(group_var_1,group_var_2), ## Select the groups, generated in line 17
  reference_gene = reference_gene, ## Define the reference gene of your table
  reference_group = reference_group, ## Define the control group, from line 17
  plot=TRUE)  ##to indicate that GAPDH and PTBP1 were measured in the same tube

results_summary_all <- results_summary_all[[1]]


results_summary_all$log_relative_expression <- log10(results_summary_all$relative_expression)
results_summary_all$log_lower <- log10(results_summary_all$lower)
results_summary_all$log_upper <- log10(results_summary_all$upper)

results_summary_all <- results_summary_all %>% group_by(group,gene) %>% 
  mutate(median=median(relative_expression),log_median=median(log_relative_expression))


my_comparisons <- list( c("ESR1 KO AP","NTC"),  c("NTC" ,"ESR1 KO MP") )

#plot only the ESR1 Common 
library(ggpubr)

p <- results_summary_all %>% filter(gene=="ESR1_P1") %>% select(c(group,relative_expression,
                                                                   lower,upper)) %>%
  unique() %>% 
  ggplot + geom_col(aes(x=group,y=relative_expression)) + yscale("sqrt")
  
p + geom_errorbar(aes(x=group,ymin = lower, ymax = upper), width = 0.2)
#add stats to the P1 P2 


p <- results_summary_all %>% filter(gene=="ESR1_Common") %>% select(c(group,relative_expression,
                                                                  lower,upper)) %>%
  unique() %>% 
  ggplot + geom_col(aes(x=group,y=relative_expression)) + yscale("log10")
p + geom_errorbar(aes(x=group,ymin = lower, ymax = upper), width = 0.2)


results_summary_all_p1 <- results_summary_all %>% filter(gene=="ESR1_P1") %>% select(c(group,relative_expression,
                                                             lower,upper)) %>% unique() 
#
t.test()
p <- ggbarplot(results_summary_all_p1,x="group","relative_expression", color = "group", palette = "Paired", title="ESR1_P1")+
  labs(fill="Promoter KD") + xlab("KD of Interest") + ylab("Relative Expression") 
p + geom_errorbar(aes(x=group,ymin = lower, ymax = upper), width = 0.2)


