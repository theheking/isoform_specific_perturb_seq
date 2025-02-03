##Read in the library
library(ggpubr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(lmtest)
library(rstatix)
library(lme4)

# Set order
order <- c("Control", "MP", "AP")
hue_order <- c('T', 'F')
order <- c("P1", "P2")

# Set the color for the histogram
color1 <- '#4d00c7'
palecolor1 <- "#b366ff"
color2 <- '#da3c07'
palecolor2 <- "#ff8954"
color3 <- '#05d3d3'
color4 <- '#c6c7c5'
color5 <- 'black'
p1 <- "#e6b821"
p2 <- "#f05223"
control <- "#03afbb"
# Create the color palette
palette <- c(palecolor1, palecolor2, color3)
palette <- c(p1, p2, control)
new_palette <- c(color1, color2, color3, palecolor1, palecolor2, color4, color5, color2, color1, color2, color1, color2, color3, color4)
slope <- function(x1, y1, x2, y2) {
  if(x2 - x1 != 0) {
    return ((y2-y1) / (x2-x1))
  } else {
    return (Inf)
  }
}
#convert this to R
##import the data 
esr1_folder <- 'alt-prom-crispr-fiveprime/files/esr1/'
#read in esr1 
esr1 <- read.table(paste0(esr1_folder, '20242003_ERS1_Prolassay.txt'), header = T, sep = '\t')
esr1



# drop the date time column
esr1 <- esr1 %>% select(-`Date.Time`)

# melt the dataframe
esr1 <- esr1 %>% pivot_longer(cols = -Elapsed, names_to = "variable", values_to = "value")

# change the experimental conditions
esr1 <- esr1 %>% mutate(
  treatment = ifelse(str_detect(variable, "Tamoxifen"), "Tamoxifen", "Untreated"),
  perturbation = case_when(
    str_detect(variable, "wt") ~ "WT",
    str_detect(variable, "dCas9") ~ "dCas9",
    str_detect(variable, "NTC") ~ "NTC",
    str_detect(variable, "MP1") ~ "MP1",
    str_detect(variable, "MP2") ~ "MP2",
    str_detect(variable, "AP1") ~ "AP1",
    str_detect(variable, "AP2") ~ "AP2",
    str_detect(variable, "BFP") ~ "BFP",
    TRUE ~ variable
  ),
  prom_perturbation = case_when(
    str_detect(variable, "wt") ~ "Control",
    str_detect(variable, "dCas9") ~ "Control",
    str_detect(variable, "NTC") ~ "Control",
    str_detect(variable, "MP1") ~ "MP",
    str_detect(variable, "MP2") ~ "MP",
    str_detect(variable, "AP1") ~ "AP",
    str_detect(variable, "AP2") ~ "AP",
    str_detect(variable, "BFP") ~ "Control",
    TRUE ~ variable
  ),
  standarderror = ifelse(str_detect(variable, "Std.Err.Img"), "Std.Err.Img", "Value"),
  prom_perturbation_treatment = paste(prom_perturbation, treatment, sep = "_"),
  perturbation_treatment = paste(perturbation, treatment, sep = "_")
)


# recreate the original plot
esr1_value <- esr1 %>% filter(standarderror == "Value")

#create a new column with well take the final three values from variable  remove .
esr1_value$well <- esr1_value$variable %>% str_sub(-3) %>% str_remove("\\.")

#group by the well and subtract the minimum value for value per well
esr1_value$mean <- esr1_value %>% mutate(mean = mean(value)) 
##calculate the per group doubling time 
# Fit linear model for each group and summarize results
esr1_value$Elapsed <- as.numeric(esr1_value$Elapsed)
esr1_value$well <- as.factor(esr1_value$well)
esr1_value$value <- as.numeric(esr1_value$value)


# Define the function to find time when cell confluence is double the minimum value
find_time_double_min_confluence <- function(data, time_col, confluence_col) {
  # Compute the minimum cell confluence value
  min_confluence <- min(data[[confluence_col]], na.rm = TRUE)
  # Compute double the minimum value
  double_min_confluence <- 2 * min_confluence
  # Find the rows where cell confluence crosses the doubled value
  crossing_points <- data %>%
    arrange(!!sym(time_col)) %>%
    mutate(next_confluence = lead(!!sym(confluence_col)),  # Get the cell confluence value of the next row
           next_time = lead(!!sym(time_col))) %>%       # Get the time value of the next row
    filter(!!sym(confluence_col) <= double_min_confluence & next_confluence >= double_min_confluence) %>%
    select(!!sym(time_col), next_time, !!sym(confluence_col), next_confluence)
  
  # Interpolate the time where cell confluence is double the minimum value
  result <- crossing_points %>%
    rowwise() %>%
    mutate(interpolated_time = !!sym(time_col) + (double_min_confluence - !!sym(confluence_col)) / 
             (next_confluence - !!sym(confluence_col)) * 
             (next_time - !!sym(time_col))) %>%
    select(interpolated_time)
  
  # Return the result
  return(result)
}



# Define the function to find time when cell confluence is double the minimum value
double_min_confluence <- function(data, time_col, confluence_col) {
  # Compute the minimum and maximum cell confluence values
  min_confluence <- min(data[[confluence_col]], na.rm = TRUE)
  max_confluence <- max(data[[confluence_col]], na.rm = TRUE)
  
  # Compute double the minimum value
  double_min_confluence <- 2 * min_confluence
  
  # Find the time and confluence values corresponding to the minimum and maximum
  min_time <- data %>% filter(!!sym(confluence_col) == min_confluence) %>% pull(!!sym(time_col))
  max_time <- data %>% filter(!!sym(confluence_col) == max_confluence) %>% pull(!!sym(time_col))
  
  # Calculate the time when cell confluence is double the minimum value using the logarithmic formula
  interpolated_time <- min_time + (log(2) * (max_time - min_time)) / 
    log(max_confluence / min_confluence)
  
  # Return the result as a data frame
  result <- data.frame(interpolated_time = interpolated_time)
  
  return(result)
}

# Define a function to apply to each group
doubling_time <- esr1_value %>%
    group_by(well,perturbation_treatment,prom_perturbation_treatment,treatment,perturbation) %>% # Group by all columns except for time and confluence
    do(double_min_confluence(., "Elapsed", "value")) %>%
    ungroup()
#create another perturbation column splitting prom_perturbation_treatment by _
doubling_time$perturbation_prom <- doubling_time$prom_perturbation_treatment %>%
  str_split_fixed("_", n = 2) %>% #select the second
  as.data.frame() %>%
  pull(1)

#summarise the doubling time per perturbation_treatment
#plot the doubling time  doubling_time boxplot

#change the type of the interpolated time to numeric
doubling_time$interpolated_time <- as.numeric(doubling_time$interpolated_time)
p <- doubling_time %>% filter(perturbation!="BFP") %>% 
  ggboxplot(x = "perturbation_prom", y = "interpolated_time", 
               color = "perturbation_prom", palette = "npg",
               add = "jitter",
               facet.by = "treatment",  short.panel.labs = FALSE)
# p
# Use only p.format as label. Remove method name.
p + stat_compare_means( ref.group = "Control",
  aes(label = paste0("p = ", after_stat(p.format)))
)

##simplfy tbelow

doubling_time <- doubling_time %>%
  mutate(perturbation_treatment = case_when(
    perturbation_treatment %in% c("AP1_Tamoxifen", "AP2_Tamoxifen") ~ "AP_Tamoxifen",
    perturbation_treatment %in% c("MP1_Tamoxifen", "MP2_Tamoxifen") ~ "MP_Tamoxifen",
    perturbation_treatment %in% c("AP1_Untreated", "AP2_Untreated") ~ "AP_Untreated",
    perturbation_treatment %in% c("MP1_Untreated", "MP2_Untreated") ~ "MP_Untreated",
    TRUE ~ perturbation_treatment
  ))

doubling_time %>% 
  group_by(perturbation_treatment)%>% filter(perturbation!="BFP")  %>% 
  summarise(mean_doubling_time = mean(interpolated_time, na.rm = TRUE), 
            sd_doubling_time = sd(interpolated_time, na.rm = TRUE))
# Create a lm plot
ggplot(esr1_value, aes(x = Elapsed, y = value, color = perturbation_treatment)) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  labs(x = "Time (h)", y = "Phase Object Confluence (%)") 

esr1_value <- esr1_value %>% group_by(well) %>% mutate(value = value - min(value))
esr1_value <- esr1_value %>%group_by(perturbation_treatment,Elapsed) %>% mutate(mean = mean(value))

# Create a lm plot
ggplot(esr1_value, aes(x = Elapsed, y = value, color = perturbation_treatment)) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  labs(x = "Time (h)", y = "Phase Object Confluence (%)") 

formula <- y ~ poly(x, 2, raw = TRUE)

#repeated the Untreated
t <- esr1_value %>% filter(treatment == "Untreated") %>% 
  ggplot( ) +
  geom_point(aes(x = Elapsed, y = mean, color = prom_perturbation)) +
  geom_smooth(aes(x = Elapsed, y = value, color = prom_perturbation),method = "lm", formula = y ~ poly(x, 2), se = TRUE) +
  labs(x = "Time (h)", y = "Phase Object Confluence (%)") + scale_color_manual(values=palette) +
  stat_regline_equation(
    aes(x = Elapsed, y = value, color = prom_perturbation, label =  paste(..eq.label..,..adj.rr.label.., sep = "~~")),
    formula = formula
  )
t+ theme_light()

esr1_value %>% filter(treatment == "Tamoxifen") %>% 
ggline( x = "Elapsed", y = "value", 
       add = c("mean_se"),
       color = "prom_perturbation") + rotate_x_text(45)

#remove esr1_value perutrbation BFP 
esr1_value <- esr1_value %>% filter(perturbation != "BFP")

#repeated the Untreated
t <- esr1_value %>% filter(treatment == "Untreated") %>% 
  ggplot( ) +
  geom_point(aes(x = Elapsed, y = mean, color = prom_perturbation)) +
  geom_smooth(aes(x = Elapsed, y = value, color = prom_perturbation),method = "lm", formula = y ~ poly(x, 2), se = TRUE) +
  labs(x = "Time (h)", y = "Phase Object Confluence (%)") + scale_color_manual(values=palette) 
t+ theme_light()


t <- esr1_value %>% filter(treatment == "Tamoxifen") %>% 
  ggplot( ) +
  geom_point(aes(x = Elapsed, y = mean, color = prom_perturbation)) +
  geom_smooth(aes(x = Elapsed, y = value, color = prom_perturbation),method = "lm", formula = y ~ poly(x, 2), se = TRUE) +
  labs(x = "Time (h)", y = "Phase Object Confluence (%)") + scale_color_manual(values=palette) 
t+ theme_light()


##calculate the doubling time 
# esr1_value %>% mutate(doubling_time = 1/slope)



# esr1_value$prom_perturbation <- as.factor(esr1_value$prom_perturbation, levels=c("Control","MP","AP"))
#change order of levels esr1_va
# levels(esr1_value$prom_perturbation) <- c("AP", "Control","MP")
#reoder the levels 
#place first factr of esr1_value as Controlprom_perturbation #reorder
esr1_value$prom_perturbation <- relevel(factor(esr1_value$prom_perturbation), ref = "Control")

##use a LMM with and without prom_perutrnation term 
#fit the model and chec
model1 <- lmer(value ~ prom_perturbation + (1|well), data=esr1_value%>% filter(treatment == "Tamoxifen"))
#repeat with lm 
model1 <- lm(value ~ prom_perturbation, data=esr1_value%>% filter(treatment == "Tamoxifen"))
summary(model1)


model2 <- lmer(value ~ prom_perturbation + (1|well), data=esr1_value%>% filter(treatment == "Untreated"))
summary(model2)
#repeat with lm 
model2 <- lm(value ~ prom_perturbation, data=esr1_value%>% filter(treatment == "Untreated"))
#set base factor as control
summary(model2)

#calculate the confidence intervanl
confint(model1)
confint(model2)
#calculate the slope for every min and max value per group
slope_df <- esr1_value %>% group_by(perturbation_treatment,treatment,perturbation,prom_perturbation,well) %>% summarise(
  min_time = min(Elapsed),
  max_time = max(Elapsed),
  min_value = min(value),
  max_value = max(value),
  slope = slope(min_time, min_value, max_time, max_value)
)

#create a new ciolumn with prom_perturbation and treatment
slope_df$prom_perturbation_treatment <- paste(slope_df$prom_perturbation, slope_df$treatment, sep = "_")

#create a displot
#add more complimentary colours to palette
palette <- c(  p1,  p1,p2,  p2 ,control,control)
#order rhe data AP untreated tamoxifen, MP untreated tamox
levels(slope_df$prom_perturbation_treatment) <- c("AP_Untreated", "MP_Untreated", "Control_Untreated", "AP_Tamoxifen", "MP_Tamoxifen", "Control_Tamoxifen")
slope_df %>% 
  filter((perturbation!="dCas9") & (perturbation!="WT" )) %>%
  ggplot(aes(x = slope, fill = prom_perturbation_treatment)) +
  # geom_histogram(aes(y = ..density..), bins = 8, alpha = 0.1, color = "black") +
  geom_density(alpha = 0.7,adjust = 2) +
  labs(x = "Slope", y = "Density")  +
  theme_minimal() 
#calculate the significance between different groups
#compare sloped between intreated and tampxofen 


# slope_df %>% filter((perturbation!="dCas9") & (perturbation!="WT")) %>%
  # group_by(prom_perturbation_treatment) %>% get_summary_stats(slope, type = "mean_sd")

pwc <- slope_df %>% ungroup() %>%
  filter((perturbation!="dCas9") & (perturbation!="WT" )) %>%
  select(slope, prom_perturbation_treatment)  %>%
  pairwise_t_test(slope ~ prom_perturbation_treatment, p.adjust.method = "bonferroni", comparisons= list(c("AP_Untreated", "AP_Tamoxifen"), c("MP_Untreated", "MP_Tamoxifen"), c("Control_Untreated", "Control_Tamoxifen")))

res.aov <- slope_df %>% ungroup() %>%
  filter((perturbation!="dCas9") & (perturbation!="WT" )) %>%
  select(slope, prom_perturbation_treatment)  %>% 
  anova_test(slope ~ prom_perturbation_treatment)

res.aov

