#import data
upset_plot <- "../../files/guide_calling/upsetplot.csv"
upset_plot
upset_data <- read.csv(upset_plot, header = TRUE, row.names = 1)
#plot the data
#change the True to 1 and False to 0
upset_data[upset_data == "True"] <- 1
upset_data[upset_data == "False"] <- 0
#change the datatypes to integer
upset_data <- as.data.frame(lapply(upset_data, as.integer))
upset(upset_data,order.by ="freq")
head(upset_data)