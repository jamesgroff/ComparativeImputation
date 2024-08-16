# Notes: How would you deal with NA values? Why is it not in method order?
# Overall, doesn't look good as FDR do not go under 0.05 (unusable for biologists)
# Save the confusion matrices in a pdf file

# Goal: From the table of metrics and average CM given by ImputationMV.R, make figures for visualization in papers or presentations

# Disclaimer: You should run ImputationMV.R and all the required steps before running this code script so that you have all the necessary files

# First, set session or go to the directory where all the MaxQuant results are available (normally in combined/txt folder)

# General setup of packages
library(ggplot2)
library(reshape2)
library(gridExtra)

# Load the average confusion matrices for visualization in heatmap formats
load("average_CM_proteinGroups.RData")
load("average_CM_Phospho(STY)Sites.RData")
load("average_CM_Phospho(STY)Sites_normed.RData")

# Get all the metrics of proteins, phosphosites (STY and normed)
metrics_prot <- read.table("proteinGroups_avg_metrics.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)
metrics_phos <- read.table("Phospho (STY)Sites_avg_metrics.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)
metrics_norm <- read.table("Phospho (STY)Sites_normed_avg_metrics.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)

# Add method and category group vectors in those dataframes for plotting simplicity
add_methods_and_categories <- function(metrics_df){
  # Get the methods and categories from the rownames of the data itself
  methods <- unique(sub("_.*", "", rownames(metrics_df)))
  categories <- unique(sub(".*?_", "", rownames(metrics_df)))
  categories <- sub(".[A-z]+$", "",categories)
  
  metrics_df$Method <- rep(methods, each=3)
  metrics_df$Categories <- rep(categories, times=5)
  
  return(metrics_df)
}

metrics_prot <- add_methods_and_categories(metrics_prot)
metrics_phos <- add_methods_and_categories(metrics_phos)
metrics_norm <- add_methods_and_categories(metrics_norm)

# Heatmaps for confusion matrices

plot_CM <- function(confusion_matrix, method, category){
  # If the confusion matrix does not exist ("NA"), make an empty plot and break out of the function
  if (all(is.na(confusion_matrix))){
    ggplot() + theme_void()
    return()
  }
  
  # Plot the heatmap
  ggplot(confusion_matrix, aes(x=factor(Reference), y=factor(Prediction, levels=rev(levels(as.factor(Prediction)))), fill=value)) + 
    geom_tile() +
    geom_text(aes(label=format(round(value,4), nsmall=2)), color="white", size=5) +
    scale_fill_gradient(low="lightblue", high="darkblue") +
    scale_x_discrete(labels=c("Non-DE", "DE")) + 
    scale_y_discrete(labels=c("DE", "Non-DE")) + 
    labs(title=paste(category, "with", method, "method"), 
         x="Reference", y="Prediction") +
    theme_minimal() +
    theme(legend.position="none")
}

plot_list <- function(CM_list, main_title){
  # Call the heatmap plot function to plot multiple confusion matrices
  # with reshaped long format confusion matrix, method name and category name in mind
  plots <- lapply(seq_along(CM_list), function(i){
    lapply(seq_along(CM_list[[i]]), function(j){
      plot_CM(melt(CM_list[[i]][[j]]), names(CM_list)[i], names(CM_list[[i]])[j])
    })
  })
    
  # Flatten the list of plots
  plots <- unlist(plots, recursive=FALSE) 
  
  # Arrange and display plots
  do.call(grid.arrange, c(plots, ncol=3, top=paste("Average Confusion Matrices in", main_title)))
}

pdf(file="Average_CM_by_datatype.pdf", width=14, height=8)
plot_list(average_CM_overall_proteinGroups, "Protein Groups")
plot_list(average_CM_overall_Phophosites, "Phosphosites (STY)")
plot_list(average_CM_overall_normedPhosphosites, "Normed ph - ft")
dev.off()

# Barplots for metrics

barplot_metrics <- function(metrics_df){
  # Make a vector of metric names from the dataframe itself
  metric_names <- c("FDR", "Accuracy", "MCC", "NRMSE")
  sd_names <- c("FDR_sd", "Accuracy_sd", "MCC_sd", "NRMSE_sd")
  
  # Create a list to store individual plots
  plot_list <- vector("list", length=length(metric_names))
  
  # Get the categories from the dataframe
  cats <- rownames(metrics_df)
  
  # Get the categories names when doing the x-axis labeling
  categories <- unique(sub(".*?_", "", rownames(metrics_df)))
  categories <- sub(".[A-z]+$", "",categories)
  categories <- rep(categories, times=5)
  
  for (i in 1:length(metric_names))
  plot_list[[i]] <- ggplot(metrics_df, aes(x=cats, y=!!sym(metric_names[i]), fill=factor(Method))) +
    geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=!!sym(metric_names[i]) - !!sym(sd_names[i]), ymax=!!sym(metric_names[i]) + !!sym(sd_names[i])), width=.2, position=position_dodge(.9)) +
    labs(title=paste("Average", metric_names[i]), x="Imputation Method and Category Ratio", y= metric_names[i]) +
    scale_x_discrete(labels=categories) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none") +
    geom_text(aes(label=Method, y=!!sym(metric_names[i]) + !!sym(sd_names[i])), vjust=-0.3, color="black", size=2)
  
  return(plot_list)
}

# Plot those barplots as arguments
protein_barplots <- barplot_metrics(metrics_prot)
Phosphosites_barplots <- barplot_metrics(metrics_phos)
normedPhosphosites_barplots <- barplot_metrics(metrics_norm)

# Save them in a pdf file (ignore warning messages, as they are for missing values to be represented)
pdf(file="Metrics_by_Methods.pdf", width=14, height=8)
grid.arrange(grobs=protein_barplots, nrow=2, ncol=2, top="Comparison of Average Metrics in Protein Groups")
grid.arrange(grobs=Phosphosites_barplots, nrow=2, ncol=2, top="Comparison of Average Metrics in Phosphosites (STY)")
grid.arrange(grobs=normedPhosphosites_barplots, nrow=2, ncol=2, top="Comparison of Average Metrics in Normed ph - ft")
dev.off()
