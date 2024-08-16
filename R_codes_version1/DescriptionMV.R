# Important Note: From the ratios and intensities histograms given, which of the two intensities from the ratio are missing (color them by both, one, none),
# and of those intensities that have only one missing, which ones are top hits and which ones are not.

# Goal: Describe through a matrix or plots the distribution of missing values from MaxQuant datas in proteomes and phosphosites by labels, conditions and intensities

# Required packages to install: stringr, ggplot2, gridExtra, reshape2

# First, set session or go to the directory where all the MaxQuant results are available (normally in combined/txt folder)

#rm(list=ls())

# For conditions, create a txt/excel file named "overviewConditions.txt" where you indicate your groups/columns 
# where the rows are the replicates or experiments and the columns are the labels from SILAC experiments (see overview file as an example)
conditions_table_orig <- read.table("overviewConditions.txt", header=TRUE, sep="\t")

# If there is a column named experiment and not replicates (like in MB_03 or MB_04), copy each row and paste repeatedly with a given number of replicates
if ("Experiment" %in% colnames(conditions_table_orig)){
  # Prompt the user to enter the number of replicates per experiment
  replicate_prompt <- function(){
    num_replicate <- readline("How many replicates do you have per experiment? ")
    num_replicate <- as.numeric(num_replicate)
    return(num_replicate)
  }
  n_rep <- replicate_prompt()
  # Repeat each row of the conditions table n_rep times
  conditions_table_orig <- conditions_table_orig[rep(row.names(conditions_table_orig), each=n_rep),]
  # Reset row names
  row.names(conditions_table_orig) <- NULL
}

# Remove the replicate/experiment column as it is not needed for futher usage
conditions_table <- conditions_table_orig[, -1]

# Get all the datasets of proteins and phosphosites
proteinGroups <- read.table("proteinGroups.txt", header=TRUE, sep="\t", fill=TRUE, quote="")
Phosphosites <- read.table("Phospho (STY)Sites.txt", header=TRUE, sep="\t", fill=TRUE, quote="")

# If you have the Perseus combined dataset (after normalization by log2), you can also get the dataset for combined data
# It should already be filtered and ridden of undesirable rows.
normedPhosphosites <- read.table("Phospho (STY)Sites_normed.txt", header=TRUE, sep="\t", fill=TRUE, quote="")

# Get column names for naming verification and interested columns
#colnames(conditions_table)
#colnames(proteinGroups)
#colnames(Phosphosites)
#colnames(normedPhosphosites)

# Remove rows that have potential contaminants, reversed sequences and that were only identified by site
remove_rows_protein <- function(df){
  df_tmp <- subset(df, Only.identified.by.site!="+")
  df_tmp <- subset(df_tmp, Reverse!="+")
  df_tmp <- subset(df_tmp, Potential.contaminant!="+")
  
  return (df_tmp)
}

remove_rows_phos <- function(df){
  df_tmp <- subset(df, Reverse!="+")
  df_tmp <- subset(df_tmp, Potential.contaminant!="+")
  
  return (df_tmp)
}

proteinGroups <- remove_rows_protein(proteinGroups)
Phosphosites <- remove_rows_phos(Phosphosites)

# Labels bias

# Filter only normalized ratios as columns and only as separate experiments
desired_columns_lab <- names(proteinGroups)[grep("^Ratio.[H|M].[M|L].normalized.[A-z]+[0-9]+$", names(proteinGroups))]
labels_proteinGroups <- proteinGroups[,desired_columns_lab]
labels_Phosphosites <- Phosphosites[,desired_columns_lab]

# Rename them with a function to reduce the column names to its essential parts and as ratio naming
reduce_colnames <- function(df){
  colnames(df) <- gsub("normalized.", "", as.character(colnames(df)))
  colnames(df) <- gsub("Ratio.", "", as.character(colnames(df)))
  colnames(df) <- sub(".", "/", as.character(colnames(df)), fixed=TRUE)
}
colnames(labels_proteinGroups) <- reduce_colnames(labels_proteinGroups)
colnames(labels_Phosphosites) <- reduce_colnames(labels_Phosphosites)

# For the normed phosphosites, you would need to rename the columns for the function accessibility
desired_columns_norm <- names(normedPhosphosites)[grep("^Ratio.[H|M].[M|L].normalized.[A-z]+[0-9]+_ph_x.y", names(normedPhosphosites))]
labels_normedPhosphosites <- normedPhosphosites[,desired_columns_norm]
colnames(labels_normedPhosphosites) <- desired_columns_lab
colnames(labels_normedPhosphosites) <- reduce_colnames(labels_normedPhosphosites)

# Intensities bias

# Filter only intensities as columns and only as separate experiments
desired_columns_int <- names(proteinGroups)[grep("^Intensity.[H|M|L].[A-z]+[0-9]+$", names(proteinGroups))]
intensities_proteinGroups <- proteinGroups[,desired_columns_int]
intensities_Phosphosites <- Phosphosites[,desired_columns_int]

# Conditions bias based on intensity

# Function to rearrange the intensities matrix with the conditions table in mind
rearrange_matrix_intensity <- function(conditions_matrix, numeric_matrix){
  # Initialize rearranged matrix that should be the same as the numeric matrix
  rearranged_matrix <- numeric_matrix
  
  # Loop through each conditions from the conditions matrix
  for (i in 1:nrow(conditions_matrix)){
    for (j in 1:ncol(conditions_matrix)){
      # Get the condition value from the conditions table
      condition_value <- conditions_matrix[i,j]
      # Calculate the appropriate column index for replacing it in the numeric matrix by using row-major ordering
      col_index <- (i-1)*ncol(conditions_matrix) + (j-1) + 1
      # Rename the column with the condition value and the replicate number
      colnames(rearranged_matrix)[col_index] <- paste0(condition_value,".Exp",i)
    }
  }
  # Reorder columns to have ordered conditions and replicate
  rearranged_matrix <- rearranged_matrix[, order(colnames(rearranged_matrix))]
  
  return(rearranged_matrix)
}

conditions_by_intensity_proteinGroups <- rearrange_matrix_intensity(conditions_table, intensities_proteinGroups)
conditions_by_intensity_Phosphosites <- rearrange_matrix_intensity(conditions_table, intensities_Phosphosites)

# Conditions bias based on labels (Replace the conditions with "conditions ratios" and rearrange them with labels)
library(stringr)

# Function to generate pair combinations for a given row
get_pair_combinations <- function(row_number, row_data){
  pairs_row <- combn(row_data, 2)
  sorted_pairs <- t(apply(pairs_row, 2, str_sort)) # Sort each pair alphabetically
  combined_pairs <- apply(sorted_pairs, 1, paste, collapse="/")
  complete_pairs <- paste0(combined_pairs,".Exp", row_number)
  return(complete_pairs)
}
  
# Function to rearrange the labels matrix with the conditions table in mind
rearrange_matrix_label <- function(conditions_matrix, numeric_matrix){
  # Initialize rearranged matrix that should be the same as the numeric matrix
  rearranged_matrix <- numeric_matrix
  
  # Initialize a matrix to store pair combinations
  num_combinations <- choose(ncol(conditions_matrix),2)
  pair_combinations_matrix <- matrix("", nrow=nrow(conditions_matrix), ncol=num_combinations)
  
  for (i in 1:nrow(conditions_matrix)){
    pair_combinations_matrix[i,] <- get_pair_combinations(i, conditions_matrix[i,])
  }
  # Combine pair combinations into a single vector
  conditions_vector <- as.vector(pair_combinations_matrix)
  
  # Replace the column names from the rearranged matrix with the conditions ratio name vector
  colnames(rearranged_matrix) <- conditions_vector
  
  # Reorder columns to have ordered conditions and replicates
  rearranged_matrix <- rearranged_matrix[, order(colnames(rearranged_matrix))]
  
  return(rearranged_matrix)
}

conditions_by_label_proteinGroups <- rearrange_matrix_label(conditions_table, labels_proteinGroups)
conditions_by_label_Phosphosites <- rearrange_matrix_label(conditions_table, labels_Phosphosites)
conditions_by_label_normedPhosphosites <-rearrange_matrix_label(conditions_table, labels_normedPhosphosites)

# Bins calculations for MV distribution

# Function to calculate n+1 bins of missing value ratios for n replicates with multiple categories (like M/L, H/L or H/M for labels):
get_missing_value_bins_by_category <- function(df){
  # Replace zeroes with NA
  df[df==0] <- NA
  
  # Extract categories from column names
  categories <- unique(gsub("\\d+$","", colnames(df)))
  
  # Create a list to store the results for each category
  result_list <- list()
  
  for (category in categories){
    # Subset dataframe for the current category
    cat_cols <- grep(paste0("^", category), colnames(df))
    cat_df <- df[, cat_cols, drop = FALSE]
    
    # Calculate the number of missing values per row
    missing_counts <- rowSums(is.na(cat_df))
    
    # Determine the maximum number of missing values
    max_missing <- max(missing_counts)
    
    # Create a matrix to store the counts
    bins <- matrix(0, nrow = max_missing + 1, ncol=1)
    
    # Count the frequencies of rows with different numbers of missing values
    for (i in 0:max_missing) {
      bins[i+1,1] <- sum(missing_counts==i)
    }
    
    colnames(bins) <- category
    rownames(bins) <- paste(0:max_missing, "MV", sep=" ")
    
    # Store the results for the current category
    result_list[[category]] <- bins / nrow(cat_df)
  }
  # Convert the list into a proper matrix
  final_result <- do.call(cbind,result_list)
  
  return(final_result)
}

# If there are different labels of ratios, run the following lines
bins_by_labels_proteinGroups <- get_missing_value_bins_by_category(labels_proteinGroups)
bins_by_labels_Phosphosites <- get_missing_value_bins_by_category(labels_Phosphosites)
bins_by_labels_normedPhosphosites <- get_missing_value_bins_by_category(labels_normedPhosphosites)

# The same for intensities
bins_by_intensities_proteinGroups <- get_missing_value_bins_by_category(intensities_proteinGroups)
bins_by_intensities_Phosphosites <- get_missing_value_bins_by_category(intensities_Phosphosites)

# The same for conditions by intensities
bins_by_conditions_intensities_proteinGroups <- get_missing_value_bins_by_category(conditions_by_intensity_proteinGroups)
bins_by_conditions_intensities_Phosphosites <- get_missing_value_bins_by_category(conditions_by_intensity_Phosphosites)

# The same for conditions by labels
bins_by_conditions_labels_proteinGroups <- get_missing_value_bins_by_category(conditions_by_label_proteinGroups)
bins_by_conditions_labels_Phosphosites <- get_missing_value_bins_by_category(conditions_by_label_Phosphosites)
bins_by_conditions_labels_normedPhosphosites <- get_missing_value_bins_by_category(conditions_by_label_normedPhosphosites)

# Print the distributions values/ratios for verification
print(bins_by_labels_proteinGroups)
print(bins_by_labels_Phosphosites)
print(bins_by_labels_normedPhosphosites)

print(bins_by_intensities_proteinGroups)
print(bins_by_intensities_Phosphosites)

print(bins_by_conditions_intensities_proteinGroups)
print(bins_by_conditions_intensities_Phosphosites)

print(bins_by_conditions_labels_proteinGroups)
print(bins_by_conditions_labels_Phosphosites)
print(bins_by_conditions_labels_normedPhosphosites)

# Take notes of the MV distributions with the last two variables from your experiments for the generation of artificial matrices,
# there should not be many significant difference of ratios between the categories

# Save those MV distributions through txt files
bins_overall_proteinGroups <- cbind(bins_by_labels_proteinGroups, bins_by_conditions_labels_proteinGroups, bins_by_intensities_proteinGroups, bins_by_conditions_intensities_proteinGroups)
bins_overall_Phosphosites <- cbind(bins_by_labels_Phosphosites, bins_by_conditions_labels_Phosphosites, bins_by_intensities_Phosphosites, bins_by_conditions_intensities_Phosphosites)
bins_overall_normedPhosphosites <- cbind(bins_by_labels_normedPhosphosites, bins_by_conditions_labels_normedPhosphosites)

write.table(bins_overall_proteinGroups, file="proteinsMVdistribution.txt", sep="\t", col.names=NA)
write.table(bins_overall_Phosphosites, file="phosphositesMVdistribution.txt", sep="\t", col.names=NA)
write.table(bins_overall_normedPhosphosites, file="phosphosites_normedMVdistribution.txt", sep="\t", col.names=NA)

# There are many ways to visualize those distributions after matrix/df implementation
library(ggplot2)
library(gridExtra)
library(reshape2)

# Barplots format for each type of matrix or dataframe

plot_barplots <- function(matrix_data){
  # Get number of columns in the matrix
  num_cols <- ncol(matrix_data)
  
  # Create a list to store individual plots
  plot_list <- vector("list", length=num_cols)
  
  # Iterate over each column and create barplots
  for (i in 1:num_cols){
    plot_list[[i]] <- ggplot(data.frame(Columns=rownames(matrix_data), Proportion=matrix_data[,i]), aes(x=Columns, y=Proportion)) +
      geom_bar(stat="identity", fill="skyblue", color="black") +
      ylim(0,1) +
      labs(x=" Number of missing values", y="Proportion of MV", title=colnames(matrix_data)[i]) +
      theme_minimal()
  }
  
  return(plot_list)
}

#Create barplots
protein_barplots <- plot_barplots(bins_overall_proteinGroups)
Phosphosites_barplots <- plot_barplots(bins_overall_Phosphosites)
normedPhosphosites_barplots <- plot_barplots(bins_overall_normedPhosphosites)

# Save them in a pdf file (you could skip this if it does not work)
pdf(file="Barplot_by_data_type_MV.pdf", width=11, height=8)
grid.arrange(grobs=protein_barplots, nrow=4, ncol=ncol(bins_by_intensities_proteinGroups), top="Overall Missing Value Distributions for Protein Groups")
grid.arrange(grobs=Phosphosites_barplots, nrow=4, ncol=ncol(bins_by_intensities_Phosphosites), top="Overall Missing Value Distributions for Phosphosites (STY)")
grid.arrange(grobs=normedPhosphosites_barplots, nrow=2, ncol=ncol(bins_by_labels_normedPhosphosites), top="Overall Missing Value Distributions for Normed ph - ft")
dev.off()

# Format to compare between types (first average out the bins by category matrix), all the matrices must be in a list first

# Function to compute the average of matrices
compute_matrix_average <- function(matrix_list){
  averages <- sapply(matrix_list, function(single_matrix){
    rowMeans(single_matrix, na.rm=TRUE)
  })
  return(averages)
}

# Function to compare a given list of matrices by average of rows/replicates
make_bins_list_average <- function(bins_list){
  matrix_averages <- compute_matrix_average(bins_list)
  if (ncol(matrix_averages) == 2){
    colnames(matrix_averages) <- c("Protein Groups", "Phosphosites (STY)")
  }
  else{
    colnames(matrix_averages) <- c("Protein Groups", "Phosphosites (STY)", "Normed ph - ft")
  }
  return(matrix_averages)
}

# Make a list of all the matrices per bias type and compute its average of matrices
bins_by_labels_list <- list(bins_by_labels_proteinGroups, bins_by_labels_Phosphosites, bins_by_labels_normedPhosphosites)
bins_by_intensities_list <- list(bins_by_intensities_proteinGroups, bins_by_intensities_Phosphosites)
bins_by_conditions_labels_list <- list(bins_by_conditions_labels_proteinGroups, bins_by_conditions_labels_Phosphosites, bins_by_conditions_labels_normedPhosphosites)
bins_by_conditions_intensities_list <- list(bins_by_conditions_intensities_proteinGroups, bins_by_conditions_intensities_Phosphosites)

matrix_averages_labels <- make_bins_list_average(bins_by_labels_list)
matrix_averages_intensites <- make_bins_list_average(bins_by_intensities_list)
matrix_averages_conditions_labels <- make_bins_list_average(bins_by_conditions_labels_list)
matrix_averages_conditions_intensities <- make_bins_list_average(bins_by_conditions_intensities_list)

# Barplot visualization
average_barplots_labels <- plot_barplots(matrix_averages_labels)
average_barplots_intensities <- plot_barplots(matrix_averages_intensites)
average_barplots_conditions_labels <- plot_barplots(matrix_averages_conditions_labels)
average_barplots_conditions_intensities <- plot_barplots(matrix_averages_conditions_intensities)

pdf(file="Barplot_by_bias_type_MV.pdf", width=11, height=8)
grid.arrange(grobs=average_barplots_labels, ncol=length(average_barplots_labels), top="Missing Value Distributions for average SILAC labels ratio")
grid.arrange(grobs=average_barplots_intensities, ncol=length(average_barplots_intensities), top="Missing Value Distributions for average SILAC intensities")
grid.arrange(grobs=average_barplots_conditions_labels, ncol=length(average_barplots_conditions_labels), top="Missing Value Distributions for average SILAC conditions by labels ratio")
grid.arrange(grobs=average_barplots_conditions_intensities, ncol=length(average_barplots_conditions_intensities), top="Missing Value Distributions for average SILAC conditions by intensities")
dev.off()

# Heatmap format overall

plot_heatmap <- function(matrix_data, main_title){
  # Reshape the matrix into a long format
  missing_data <- melt(matrix_data)
  
  # Plot the heatmap of the matrix
  ggplot(missing_data, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile() +
    geom_text(aes(label=format(round(value, 4), nsmall=4)), color="white", size=5) +
    scale_fill_gradient(low="lightblue", high="darkblue") +
    guides(fill=guide_colorbar(title="Proportion", barwidth=0.5, barheight=20)) +
    labs(x="Data type", y="Number of missing values", title=paste("Heatmap of Missing Value Distribution for", main_title)) +
    theme_minimal()
}

# Plot heatmaps per type of data between biases overall
pdf(file="Heatmap_by_data_type_MV.pdf", width=13, height=8)
plot_heatmap(bins_overall_proteinGroups, "Protein Groups overall")
plot_heatmap(bins_overall_Phosphosites, "Phosphosites (STY) overall")
plot_heatmap(bins_overall_normedPhosphosites, "Normed ph - ft overall")
dev.off()

# Plot heatmaps per type of bias between datasets
pdf(file="Heatmap_by_bias_type_MV.pdf", width=11, height=8)
plot_heatmap(matrix_averages_labels, "labels average")
plot_heatmap(matrix_averages_intensites, "intensities average")
plot_heatmap(matrix_averages_conditions_labels, "conditions by labels average")
plot_heatmap(matrix_averages_conditions_intensities, "conditions by intensities average")
dev.off()

# Plot histograms per dataframe (once for global in blue, then with averaged missing values in red)

plot_histogram <- function(dataframes, measure_type, pdf_name){
  # Start the saving of pdf files
  pdf(file=pdf_name, width=11, height=8)
  
  # Set up the dataframes index for later iterations
  df_index <- 0
  
  for (df in dataframes){
    # Add 1 to the dataframes index
    df_index <- df_index + 1
    
    # Replace zeroes with NA
    df[df==0] <- NA
    
    # Extract categories from column names and adjust mfrow as such
    categories <- unique(gsub("\\d+$","", colnames(df)))
    par(mfrow=c(length(categories), 1), mar=rep(2,4))
    
    # Choose datatypes depending on the length of dataframes list
    if (length(dataframes) == 2){
      data_type <- c("Protein Groups", "Phosphosites (STY)")
    }
    else{
      data_type <- c("Protein Groups", "Phosphosites (STY)", "Normed ph - ft")
    }
    
    for (category in categories){
      # Subset dataframe for the current category
      cat_cols <- grep(paste0("^", category), colnames(df))
      cat_df <- df[, cat_cols, drop = FALSE]
      
      # Rename the category without the ".Exp"
      new_category <- gsub("(.*)\\..*","\\1", category)
      
      # Log2 transform the data before anything (except for normed data as it is already done)
      if (data_type[df_index] != "Normed ph - ft"){
        log2_df <- sapply(cat_df, log2)
      }
      else{
        log2_df <- cat_df
      }
      
      # Calculate means per row
      mean_df <- apply(log2_df, 1, function(x) mean(na.omit(x)))
      
      # Put everything in one summary dataframe
      summary_df <- data.frame(cat_df, average=mean_df, MV=rowSums(is.na(cat_df)))
      
      # Extract data that have missing values
      MV_df <- summary_df$average[summary_df$MV>0]
      
      # Build histogram for the current category global distribution
      hist_df <- hist(mean_df, breaks="Scott", col="blue",
                      main=paste("Histogram of", new_category, "for", data_type[df_index]), xlab=paste(measure_type, new_category, "measure"))
      
      # Build histogram for missing values with the same bins from the global distribution
      hist(MV_df, breaks=hist_df$breaks, col="red", add=TRUE)
      
      # Add a legend to indicate distribution difference
      legend("topright", legend=c("Global", "MV only"), fill=c("blue", "red"))
    }
  }
  
  # Close the pdf device
  dev.off()
}

# Create list of dataframes by labels, intensities and conditions
df_by_labels_list <- list(labels_proteinGroups, labels_Phosphosites, labels_normedPhosphosites)
df_by_intensities_list <- list(intensities_proteinGroups, intensities_Phosphosites)
df_by_conditions_labels_list <- list(conditions_by_label_proteinGroups, conditions_by_label_Phosphosites, conditions_by_label_normedPhosphosites)
df_by_conditions_intensities_list <- list(conditions_by_intensity_proteinGroups, conditions_by_intensity_Phosphosites)

# Plot histograms for every possible dataset that is applicable
plot_histogram(df_by_labels_list, "Ratio", "Histogram_by_labels_MV.pdf")
plot_histogram(df_by_intensities_list, "", "Histogram_by_intensities_MV.pdf")
plot_histogram(df_by_conditions_labels_list, "Ratio", "Histogram_by_conditions_labels_MV.pdf")
plot_histogram(df_by_conditions_intensities_list, "", "Histogram_by_conditions_intensities_MV.pdf")
