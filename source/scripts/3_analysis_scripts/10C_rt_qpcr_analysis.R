library(shiny)
library(stringr)
library(dplyr)
library(pcr)  ## For qPCR analyis
library(ggplot2)  ##For graphs
library(readxl)
library(rlang)
library(tidyr)
library(DT)
###!!!!Functions!!!!#####
#This provides the column names used in the uiOutput
colnames_readin_quantstudio <- function(excel_filename){
  tryCatch(
    { results<-suppressMessages(readxl::read_excel(excel_filename,sheet="Results", col_names = FALSE))
    ##extract the column names and 
    result_columnnames <- results %>% slice(which(results[1] == "Well"))
    #change the name of the columns
    names(result_columnnames) <- NULL
    result_columnnames <- unlist(c(result_columnnames))
    
    #assign them to the new dataframe
    results <-data.frame(results[-(1:(which(results[1] == "Well"))),])
    colnames(results) <-result_columnnames
    colnames_readin_list <- colnames(results)
    # print(paste0("colnames readin list as follow:",colnames_readin_list))
    },
    error = function(e) {
      # return a safeError if a parsing error occurs
      stop(safeError(e))
      colnames_readin_list <- "Attach file"
    })
  
  return(colnames_readin_list)
}

#This creates the results frame 
pcr_prep_grazi <- function(results,zscore_thresh, reference_gene, reference_group, sample_column_name,target_name, ct_column_name, return_data="normal") {   #read in excel file as input
  #Read in results frame
  # print(paste("pcrprep",sample_column_name,target_name, ct_column_name))
  results <- results %>% drop_na(ct_column_name) %>% drop_na((!!sym(sample_column_name))) %>% drop_na((!!sym(target_name))) %>% drop_na(ct_column_name)
  results$CT_VALUE <- suppressWarnings(sapply(results[ct_column_name], as.numeric)[,1])
  results <- results %>% drop_na(CT_VALUE)
  
  ##Create a z-score
  results_zscore <- results %>% group_by((!!sym(sample_column_name)), (!!sym(target_name))) %>% 
    mutate(z_score = unlist(scale(CT_VALUE)))
  
  ##Create a df showing the stupid zscore
  results_zscore <- results_zscore %>% filter(abs(z_score)< zscore_thresh)
  
  
  ##Create a new ID with Repeats ect.
  results_zscore["Sample_ID"] <- results_zscore %>% ungroup()%>%  select((!!sym(sample_column_name)))
  
  results_zscore <- results_zscore %>%
    group_by(Sample_ID, (!!sym(target_name))) %>%
    dplyr::summarize(CT_VALUE = mean(CT_VALUE, na.rm=TRUE),.groups = 'drop')
  # results_zscore$Sample_ID<- unlist(lapply(results_zscore$Sample_ID, function(x) str_remove(x, "_RepIII$|_RepI$|_RepII$|_RepIV$|_RepVI$| Rep$| RepVII$| RepIII$| RepI$| RepII$| RepIV$| RepVI$| Rep$| RepVII$") ))
  results_zscore["Sample_ID"] <-unlist(lapply(results_zscore["Sample_ID"], function(x) str_remove_all(x, "_All$|_All_FACS$|_Ex9_FACS$|_Ex9$|_Ex10$|_Ex10_FACS$|_Ex11$|_Ex12$|_FACS$|_RepIII$|_RepI$|_RepII$|_RepIV$|_RepVI$| Rep$| RepVII$| RepIII$| RepI$| RepII$| RepIV$| RepVI$| Rep$| RepVII$")))
  
  #restucture the df using pvot wider to make multiple column of the Ct value
  # select only the CT values, target name and sample column name
  #initalise columns
  results_df <- results_zscore %>% ungroup() %>% select( c(Sample_ID,(!!sym(target_name)), CT_VALUE))  %>% group_by(Sample_ID, (!!sym(target_name))) %>%
    mutate(row = row_number()) %>%
    pivot_wider( names_from = (!!sym(target_name)), values_from = CT_VALUE ) %>%
    select(-row)
  
  
  ##Define group_var as the sample id and then remove that column
  group_var <- results_df$Sample_ID
  names(group_var) <- NULL
  
  results_df <- results_df %>% ungroup() %>% select(-Sample_ID)
  if (return_data =="normal") {
    results_print <- results_df} else if (return_data =="gene"){
      results_print <- unique(colnames(results_df))
      
    } else if (return_data =="group") {
      results_print <- unique(group_var)
    }
  
  return(results_print)}

#This provides the input of the excel filename and output of the dataframe
readin_quantstudio<- function(excel_filename){
  results<-suppressWarnings(readxl::read_excel(excel_filename, sheet="Results", col_names = FALSE, trim_ws=TRUE))
  ##extract the column names and 
  result_columnnames <- results %>% slice(which(results[1] == "Well"))
  
  #change the name of the columns
  names(result_columnnames) <- NULL
  result_columnnames <- unlist(c(result_columnnames))
  
  #assign them to the new dataframe
  results <-data.frame(results[-(1:(which(results[1] == "Well"))),])
  colnames(results) <-result_columnnames
  # print(head(results))
  # #drop na
  results <- results[!is.na(results$CT),]
  return(results)}
readin <- function(excel_filename,sheet_name,input_filetype){
  if (input_filetype==1){
    results<-suppressMessages(readxl::read_excel(excel_filename,sheet = sheet_name))
    colnames_readin_list <- colnames(results)}
  if(input_filetype==2){
    results<-suppressWarnings(readxl::read_excel(excel_filename, sheet="Results", col_names = FALSE, trim_ws=TRUE))
    ##extract the column names and 
    result_columnnames <- results %>% slice(which(results[1] == "Well"))
    #change the name of the columns
    names(result_columnnames) <- NULL
    result_columnnames <- unlist(c(result_columnnames))
    
    #assign them to the new dataframe
    results <-data.frame(results[-(1:(which(results[1] == "Well"))),])
    colnames(results) <-result_columnnames
    # print(results)
    # #drop na
    results <- results[!is.na(results$CT),]}
  
  return(results)
}

#This creates the group var as output
group_var_func <- function(results,zscore_thresh, reference_gene, reference_group, sample_column_name,target_name, ct_column_name, return_data="normal") {   #read in excel file as input
  #Read in results frame
  
  
  results <- results %>% drop_na(ct_column_name) %>% drop_na((!!sym(sample_column_name))) %>% drop_na((!!sym(target_name))) %>% drop_na(ct_column_name)
  results$CT_VALUE <- suppressWarnings(sapply(results[ct_column_name], as.numeric)[,1])
  results <- results %>% drop_na(CT_VALUE)
  # print(head(results))
  
  ##Create a z-score
  results_zscore <- results %>% group_by((!!sym(sample_column_name)), (!!sym(target_name))) %>% 
    mutate(z_score = unlist(scale(CT_VALUE)))
  # print("  ##Create a z-score")
  
  ##Create a df showing the stupid zscore
  results_zscore <- results_zscore %>% filter(abs(z_score)< zscore_thresh)
  print("  ##Create a df showing the stupid zscore")
  
  
  ##Create a new ID with Repeats ect.
  results_zscore["Sample_ID"] <- results_zscore %>% ungroup()%>%  select((!!sym(sample_column_name)))
  
  results_zscore <- results_zscore %>%
    group_by(Sample_ID, (!!sym(target_name))) %>%
    dplyr::summarize(CT_VALUE = mean(CT_VALUE, na.rm=TRUE), .groups="drop")
  
  # print("    ##Create a new ID with Repeats ect.")
  
  # results_zscore$Sample_ID<- unlist(lapply(results_zscore$Sample_ID, function(x) str_remove(x, "_RepIII$|_RepI$|_RepII$|_RepIV$|_RepVI$| Rep$| RepVII$| RepIII$| RepI$| RepII$| RepIV$| RepVI$| Rep$| RepVII$") ))
  results_zscore["Sample_ID"] <-unlist(lapply(results_zscore["Sample_ID"], function(x) str_remove_all(x, "_All$|_All_FACS$|_Ex9_FACS$|_Ex9$|_Ex10$|_Ex10_FACS$|_Ex11$|_Ex12$|_FACS$|_RepIII$|_RepI$|_RepII$|_RepIV$|_RepVI$| Rep$| RepVII$| RepIII$| RepI$| RepII$| RepIV$| RepVI$| Rep$| RepVII$")))
  
  #restucture the df using pvot wider to make multiple column of the Ct value
  # select only the CT values, target name and sample column name
  #initalise columns
  results_df <- results_zscore %>% ungroup() %>% select( c(Sample_ID,(!!sym(target_name)), CT_VALUE))  %>% group_by(Sample_ID, (!!sym(target_name))) %>%
    mutate(row = row_number()) %>%
    pivot_wider( names_from = (!!sym(target_name)), values_from = CT_VALUE ) %>%
    select(-row) 
  # print("      #initalise columns.")
  
  ##Define group_var as the sample id and then remove that column
  group_var <- results_df$Sample_ID
  # print("      #  ##Define group_var as the sample id and then remove that columns.")
  
  names(group_var) <- NULL
  # print(paste0(group_var,"final group var"))
  return(group_var)
}

#This creates the pcr analysis
pcr_analysis <- function(results_summary, mode_tube,method_choice, reference_gene, reference_group, group_var) {   #read in excel file as input
  #make all input to results_summary numeric 
  results_summary <- mutate_all(results_summary, function(x) as.numeric(as.character(x)))
  results_summary <- pcr_analyze(
    df = results_summary, ## Name of the table loaded in line 13
    method = method_choice, ## Default method for analysis expression between conditions
    group_var = group_var, ## Select the groups, generated in line 17
    reference_gene = reference_gene, ## Define the reference gene of your table
    reference_group = reference_group, ## Define the control group, from line 17
    mode = mode_tube,
    plot=TRUE)  ##to indicate that GAPDH and PTBP1 were measured in the same tube
  
  print("notsuccess yet")
  
  return(results_summary)}

plot_pcr <- function(results_summary, reference_gene, reference_group, method_choice ){
  tst_Graph <- ggplot(`results_summary`, aes(x= group, y= relative_expression, fill = gene)) + ##Load the data and define X and Y from the summary (line 22)
    geom_bar(stat = "identity", position = "dodge") +
    ggtitle(paste("Using",method_choice, "to Find Expression \n Relative to ",reference_gene, "and Reference Group",reference_group )) +
    xlab("pre-gRNA") + ylab("Relative expression") + ##Title and axes, change accordingly
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
    geom_errorbar(aes(x= group, ymin=lower, ymax=upper), width=.6, position = position_dodge(.9)) #+  ## add error bars from summary
  tst_Graph
  # facet_wrap(~gene, scales = "free") ## If you want the facets. If not, remove the final + in line38 and add # at the beginning of line 30
}

pcr_testing <- function(results_readin, test_group,test_gene,stat_test, test_reference_group, sample_column_name, target_name, ct_column_name, zscore_thresh){
  ##Read in excel sheet
  results <- results_readin
  results <- results %>% drop_na(ct_column_name) %>% drop_na((!!sym(sample_column_name))) %>% drop_na((!!sym(target_name))) %>% drop_na(ct_column_name)
  results$CT_VALUE <- suppressWarnings(sapply(results[ct_column_name], as.numeric)[,1])
  results <- results %>% drop_na(CT_VALUE)
  
  ##Create a z-score
  results_zscore <- results %>% group_by((!!sym(sample_column_name)), (!!sym(target_name))) %>% 
    mutate(z_score = unlist(scale(CT_VALUE)))
  
  ##Create a df showing the stupid zscore
  results_zscore <- results_zscore %>% filter(abs(z_score)< zscore_thresh)
  
  
  ##Create a new ID with Repeats ect.
  results_zscore["Sample_ID"] <- results_zscore %>% ungroup()%>%  select((!!sym(sample_column_name)))
  
  results_zscore <- results_zscore %>%
    group_by(Sample_ID, (!!sym(target_name))) %>%
    dplyr::summarize(CT_VALUE = mean(CT_VALUE, na.rm=TRUE), .groups="drop") 
  
  results_zscore["Sample_ID"] <-unlist(lapply(results_zscore["Sample_ID"], function(x) str_remove(x, "_All$|_All_FACS$|_Ex9_FACS$|_Ex9$|_Ex10$|_Ex10_FACS$|_Ex11$|_Ex12$|_FACS$|_RepIII$|_RepI$|_RepII$|_RepIV$|_RepVI$| Rep$| RepVII$| RepIII$| RepI$| RepII$| RepIV$| RepVI$| Rep$| RepVII$")))
  
  #restucture the df using pvot wider to make multiple column of the Ct value
  # select only the CT values, target name and sample column name
  #initalise columns
  results_df <- results_zscore %>% ungroup() %>% select( c(Sample_ID,(!!sym(target_name)), CT_VALUE))  %>% group_by(Sample_ID, (!!sym(target_name))) %>%
    mutate(row = row_number()) %>%
    pivot_wider( names_from = (!!sym(target_name)), values_from = CT_VALUE ) %>%
    select(-row) 
  
  ##Define group_var as the sample id and then remove that column
  # group_var <- results_df$Sample_ID
  # names(group_var) <- NULL
  # 
  # results_df <- results_df %>% ungroup() %>% select(-Sample_ID)
  # results_df["blank"] <- 1
  # 
  # #subset based on the test group and reference group
  # results_df_test <- results_df %>% 
  #   slice(which(group_var %in% c(test_reference_group,test_group)))  %>%  ##filter dataframe for rows to do with groups
  #   select(c("blank",test_gene)) #filter for columns with control and test gnees
  # 
  #
  # group_var_test <- group_var[group_var %in% c(test_reference_group,test_group)]
  
  # print(results_df_test)
  # test_output <- pcr_test(
  #   df = results_df_test, ## Name of the table loaded in line 13
  #   method = stat_test, ## Choose the statistical test
  #   group_var = group_var_test, ## Select the groups, generated in line 17
  #   reference_gene = "blank", ## Define the reference gene of your table
  #   reference_group = test_reference_group) ## Define the control group, from line 17  
  return(results_df)}

create_results <- function(excel_filename,ct_column_name){
  results_readin<- readin_quantstudio(excel_filename=excel_filename)
  #trim white space in Sample Name 
  results_readin$`Sample Name` <- trimws(results_readin$`Sample Name`)
  #create two new columns one that contains everything before Rep and one that contains Rep
  results_readin$sample_id <- trimws(stringr::str_extract(results_readin$`Sample Name`, ".*(?=Rep)"))
  #extract only the final rep number in a column
  results_readin$`replicate_num` <- ifelse(grepl("Rep II" ,results_readin$`Sample Name`), "Rep II", "Rep I")
  
  #split into two different dataframes, run the reference on one and test the other 
  results_readin_rep1 <- results_readin %>% filter(`replicate_num` == "Rep I")
  group_var <- group_var_func(results=results_readin_rep1,  zscore_thresh=zscore_thresh, reference_gene=reference_gene, reference_group=reference_group, sample_column_name=sample_column_name,target_name=target_name, ct_column_name=ct_column_name)
  results_summary <- pcr_prep_grazi(results_readin_rep1, zscore_thresh=zscore_thresh, reference_gene=reference_gene, reference_group=reference_group, sample_column_name=sample_column_name,target_name=target_name, ct_column_name=ct_column_name)
  pcr_output_rep1 <- pcr_analysis(results_summary, group_var=group_var, mode_tube=mode_tube, reference_gene=reference_gene, reference_group=reference_group, method_choice=method_choice)
  pcr_output_rep1$replicate_num <- "Rep I"
  ##
  results_readin_rep2 <- results_readin %>% filter(`replicate_num` == "Rep II")
  group_var <- group_var_func(results=results_readin_rep2,  zscore_thresh=zscore_thresh, reference_gene=reference_gene, reference_group=reference_group, sample_column_name=sample_column_name,target_name=target_name, ct_column_name=ct_column_name)
  results_summary <- pcr_prep_grazi(results_readin_rep2, zscore_thresh=zscore_thresh, reference_gene=reference_gene, reference_group=reference_group, sample_column_name=sample_column_name,target_name=target_name, ct_column_name=ct_column_name)
  pcr_output_rep2 <- pcr_analysis(results_summary, group_var=group_var, mode_tube=mode_tube, reference_gene=reference_gene, reference_group=reference_group, method_choice=method_choice)
  pcr_output_rep2$replicate_num <- "Rep II"
  #append results
  results_summary <- rbind(pcr_output_rep1, pcr_output_rep2)
  results_summary$`promoter_kd` <- ifelse(grepl("AP" ,results_summary$group), "ESR1 KO AP", results_summary$group)
  results_summary$`promoter_kd` <- ifelse(grepl("MP" ,results_summary$promoter_kd), "ESR1 KO MP", results_summary$promoter_kd) 
  return(results_summary)} 

###!!!!Run!!!!#####

excel_filename_1="alt-prom-crispr-fiveprime/files/esr1/2024-04-29_ESR1_KO.xls"
excel_filename_2="alt-prom-crispr-fiveprime/files/esr1/2024-07-12ESR1KO.xls"

# reference_group="Wt dcas9"
# test_reference_group="Wt dCas9"
# refgroup="Wt dCas9"
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
test_group="ESR1 KO MP1"

#apply create_resuilts to c(excel_filename_1, excel_filename_2)
results_summary <- create_results(excel_filename = excel_filename_1,ct_column_name=ct_column_name)
#add column bio_rep
results_summary <- results_summary[[1]]
results_summary$bio_rep <- "Rep I"
results_summary_2 <- create_results(excel_filename = excel_filename_2,ct_column_name=ct_column_name)
results_summary_2 <- results_summary_2[[1]]
results_summary_2$bio_rep <- "Rep II"
results_summary <- rbind(results_summary, results_summary_2)
#change ESR1 KO AP1 and ESR1 KO AP2 rt ESR1 KO AP same with MP
results_summary$promoter_kd <- results_summary$group
results_summary$promoter_kd <- ifelse(results_summary$group=="ESR1 KO AP1","ESR1 KO AP" , results_summary$promoter_kd)
results_summary$promoter_kd <- ifelse(results_summary$group=="ESR1 KO AP2","ESR1 KO AP" , results_summary$promoter_kd)
results_summary$promoter_kd <- ifelse(results_summary$group=="ESR1 KO MP1","ESR1 KO MP" , results_summary$promoter_kd)
results_summary$promoter_kd <- ifelse(results_summary$group=="ESR1 KO MP2","ESR1 KO MP" , results_summary$promoter_kd)
results_summary$promoter_kd <- factor(results_summary$promoter_kd, levels=c("Wt dcas9","NTC" , "ESR1 KO MP", "ESR1 KO AP"))

# results_summary <- create_results(excel_filename = excel_filename_1,ct_column_name=ct_column_name)
##create a per gene 
tst_Graph <- `results_summary`  %>% filter(bio_rep=="Rep I")  %>% filter(gene != "ESR1_206") %>% filter(promoter_kd != "Wt dcas9") %>% ggplot() + ##Load the data and define X and Y from the summary (line 22)
  geom_bar(aes(x= promoter_kd, y= relative_expression, fill = gene), stat = "identity", position = "dodge") + 
  ggtitle(paste("Delta Delta CT to Find Expression \n Relative to ",reference_gene, "and Reference Group NTC" )) +
  theme_bw()+ scale_fill_brewer(palette = "Paired") +
  xlab("KD of Interest") + ylab("Relative Expression") + ##Title and axes, change accordingly
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_errorbar(aes(x= promoter_kd, ymin=lower, ymax=upper), width=.6, position = position_dodge(.9))+
  facet_wrap(~gene, scales = "free")  
  ##change palette ggplot
tst_Graph 
##create a per gene 
tst_Graph <- ggplot(`results_summary`) + ##Load the data and define X and Y from the summary (line 22)
  geom_bar(aes(x= promoter_kd, y= relative_expression, fill = bio_rep), stat = "identity", position = "dodge") + 
  ggtitle(paste("Delta Delta CT to Find Expression \n Relative to ",reference_gene, "and Reference Group NTC" )) +
  theme_bw()+ scale_fill_brewer(palette = "Paired") +
  xlab("KD of Interest") + ylab("Relative Expression") + ##Title and axes, change accordingly
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_errorbar(aes(x= promoter_kd, ymin=lower, ymax=upper), width=.6, position = position_dodge(.9))+
  facet_wrap(~gene, scales = "free")  
##change palette ggplot
tst_Graph 
tst_Graph <- ggplot(`results_summary`, aes(x= promoter_kd, y= relative_expression, fill = gene)) + ##Load the data and define X and Y from the summary (line 22)
geom_boxplot()+scale_fill_brewer(palette = "Paired") +  xlab("KD of Interest") + ylab("Relative Expression") +
  facet_wrap(~bio_rep, scales = "free") ##  ##Title and axes, change accordingly
tst_Graph

tst_Graph <- `results_summary` %>% filter(bio_rep=="Rep II") %>% filter(gene != "ESR1_206") %>% filter(promoter_kd != "Wt dcas9") %>% ggplot( aes(x= promoter_kd, y= relative_expression, fill = gene)) + ##Load the data and define X and Y from the summary (line 22)
  geom_boxplot()+scale_fill_brewer(palette = "Paired") +  xlab("KD of Interest") + ylab("Relative Expression") 
tst_Graph


library(ggpubr)
my_comparisons <- list( c("ESR1 KO AP","Wt dcas9"),  c("NTC" ,"Wt dcas9") ,c("Wt dcas9", "ESR1 KO MP"))

results_summary$promoter_kd <- factor(results_summary$promoter_kd, levels=c("Wt dcas9","NTC" , "ESR1 KO MP", "ESR1 KO AP"))

p1 <-  results_summary %>% filter(promoter_kd!= "Wt dcas9") %>% filter(gene=="ESR1_P1") %>% 
  ggboxplot( x = "promoter_kd", y = "relative_expression",
          color = "promoter_kd", palette = "Paired", title="ESR1_P1")+
  labs(fill="Promoter KD") + xlab("KD of Interest") + ylab("Relative Expression") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") #log yaxis
  #change order of the boxplot
#one sample t-test
t.test(results_summary$relative_expression[results_summary$promoter_kd=="ESR1 KO MP" & results_summary$gene=="ESR1_P1"], mu=1, alternative="two.sided")
t.test(results_summary$relative_expression[results_summary$promoter_kd=="ESR1 KO AP" & results_summary$gene=="ESR1_P1"], mu=1, alternative="two.sided")
p1

p2 <-  results_summary %>% filter(promoter_kd!= "Wt dcas9")  %>% filter(gene=="ESR1_206") %>%
  ggboxplot( x = "promoter_kd", y = "relative_expression",
             color = "promoter_kd", palette = "Paired", title="ESR1_206") +
  labs(fill="Promoter KD") + xlab("KD of Interest") + ylab("Relative Expression") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")  #change order of the boxplot
t.test(results_summary$relative_expression[results_summary$promoter_kd=="ESR1 KO MP" & results_summary$gene=="ESR1_206"], mu=1, alternative="two.sided")
t.test(results_summary$relative_expression[results_summary$promoter_kd=="ESR1 KO AP" & results_summary$gene=="ESR1_206"], mu=1, alternative="two.sided")



p2
p3 <-  results_summary %>% filter(gene=="ESR1_Common") %>% filter(promoter_kd != "Wt dcas9") %>%
  ggboxplot( x = "promoter_kd", y = "relative_expression",
             color = "promoter_kd", palette = "Paired", title="ESR1_Common") +
  labs(fill="Promoter KD") + xlab("KD of Interest") + ylab("Relative Expression") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")  #change order of the boxplot
p3
t.test(results_summary$relative_expression[results_summary$promoter_kd=="ESR1 KO MP" & results_summary$gene=="ESR1_Common"], mu=1, alternative="two.sided")
t.test(results_summary$relative_expression[results_summary$promoter_kd=="ESR1 KO AP" & results_summary$gene=="ESR1_Common"], mu=1, alternative="two.sided")




my_comparisons <- list( c("NTC","ESR1 KO AP"),c("NTC", "ESR1 KO MP"))
#make a my_coparisons list include Rep I and Rep II
# my_comparisons <- list(c(c("Rep I","ESR1 KO AP"),c("Rep I","NTC")), c(c("NTC","Rep I"),c("ESR1 KO MP","Rep I")), c(c("ESR1 KO AP","Rep II"),c("NTC","Rep II")), c(c("NTC","Rep II"),c("ESR1 KO MP","Rep II")))
p1 <-  results_summary %>% filter(gene=="ESR1_P1") %>% filter(bio_rep=="Rep II") %>% filter(promoter_kd != "Wt dcas9") %>% 
  ggboxplot( x = "promoter_kd", y = "relative_expression", 
             color = "promoter_kd",add = "jitter", palette = "Paired", title="ESR1_P1")+
  labs(fill="Promoter KD") +yscale("log2") + xlab("KD of Interest") + ylab("Relative Expression") +
  stat_compare_means(ref.group = "NTC",aes(label = after_stat(p.signif)), method = "t.test",paired=FALSE )  #change order of the boxplot

p2 <-  results_summary %>% filter(gene=="ESR1_206") %>%filter(bio_rep=="Rep II") %>% filter(promoter_kd != "Wt dcas9") %>%
  ggboxplot( x = "promoter_kd", y = "relative_expression",
             color = "promoter_kd",add = "jitter",palette = "Paired", title="ESR1_206") +
  labs(fill="Promoter KD") + xlab("KD of Interest") + ylab("Relative Expression") +
  stat_compare_means(ref.group = "NTC",aes(label = after_stat(p.signif)), method = "t.test",paired=FALSE) #change order of the boxplot

p3 <-  results_summary %>% filter(gene=="ESR1_Common") %>% filter(bio_rep=="Rep II") %>% filter(promoter_kd != "Wt dcas9") %>%
  ggboxplot( x = "promoter_kd", y = "relative_expression",
             color = "promoter_kd", add = "jitter",palette = "Paired", title="ESR1_Common") +
  labs(fill="Promoter KD") + xlab("KD of Interest") + ylab("Relative Expression") +
  stat_compare_means(ref.group = "NTC",aes(label = after_stat(p.signif)), method = "t.test",paired=FALSE ) #change order of the boxplot
p1
p2
p3


#get the input resultsfrom the excel file
results_readin<- readin_quantstudio(excel_filename=excel_filename_2)

#melt df to have target name as columns 
results_readin_wider <- results_readin %>% select(c(CT, `Sample Name`, `Target Name`)) %>%
  pivot_wider(names_from = `Target Name`, values_from = CT) %>%
  tidyr::unnest(cols = -`Sample Name`) %>% #split a column into multiple columns one with the rep
  mutate(`Sample Name` = trimws(`Sample Name`)) #remove white space)
#split the string befo
#select the MP1 
results_readin_wider_filt <- results_readin_wider %>% 
  filter(!`Sample Name` %in% c("Wt dcas9 Rep I","Wt dcas9 Rep II")) 
# %>%
results_readin_wider_filt$rep <- unlist(trimws(str_split_fixed(results_readin_wider_filt$`Sample Name`, "Rep", 2)[,2]))
#filter for II
results_readin_wider_filt <- results_readin_wider_filt %>% filter(rep=="II")
  # select(c(`Sample Name`,ESR1_Common,hGAPDH)) #remove the sample name column
group_var <- unlist(trimws(str_split_fixed(results_readin_wider_filt$`Sample Name`, "Rep", 2)[,1]))
#replace every ESR1 KO MP2 with ESR1 KO MP1
# group_var <- ifelse(group_var=="ESR1 KO AP2", "ESR1 KO AP1", group_var)
# group_var <- ifelse(group_var=="ESR1 KO MP2", "ESR1 KO MP1", group_var)
#change the type of data column 
results_readin_wider_filt$ESR1_Common <- as.numeric(results_readin_wider_filt$ESR1_Common)
results_readin_wider_filt$hGAPDH <- as.numeric(results_readin_wider_filt$hGAPDH)
results_readin_wider_filt$ESR1_P1 <- as.numeric(results_readin_wider_filt$ESR1_P1)
results_readin_wider_filt$ESR1_206 <- as.numeric(results_readin_wider_filt$ESR1_206)
res <- pcr_analyze(results_readin_wider_filt  %>% select(-c(`Sample Name`,rep)) ,
                   group_var = group_var,
                   reference_group = 'NTC',
                   method = 'delta_delta_ct',
                   reference_gene="hGAPDH",
                   plot = TRUE)
res

res <- pcr_analyze(results_readin_wider_filt  %>% select(-c(`Sample Name`,rep)) ,
                   group_var = group_var,
                   reference_group = 'NTC',
                   method = 'delta_delta_ct',
                   reference_gene="hGAPDH")
res

