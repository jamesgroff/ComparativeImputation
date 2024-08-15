# Notes: Evaluation of consistencies are done one-sided only (from intensities to ratios)

# Remark: The inconsistencies you will observe in the intensities to ratios assignments can be due to either:
# - peptides that have two ratios considered in proteins cases when MaxQuant analysis is done
# - isoforms that have only the highest conformation are considered in Phosphosites cases

# Goal: Describe through a matrix or plots the distribution of missing values from MaxQuant datas in proteomes and phosphosites by labels, conditions and intensities
# Additionally, give a glimpse at which of the intensities are missing in the ratios through histograms and its basic metrics

# Required packages to install: stringr, ggplot2, gridExtra, reshape2, parameters

# Before starting this script, please create a txt/excel file named "overviewConditions.txt" where you indicate your groups/columns for your conditions,
# You can find examples of such condition tables in "overview.docx"

### --- Parameters---
# Please modify the parameters below to personalize your outputs needed with either TRUE or FALSE statements

# Clean up the current environment, recommended if you want to start over
#rm(list=ls())

# Indicate if you have these MaxQuant outputs from your directory
get_proteinGroups <- TRUE
get_Phosphosites <- TRUE
# Indicate if you also have the Perseus combined dataset of normed phosphosites 
# (only possible when proteinGroups and Phosphosites data are there, otherwise the script will stop with an error)
get_normedPhosphosites <- FALSE

# Does your conditions table have replicates as rows?
have_replicates <- TRUE
# If not, then indicate the number of replicates here
nb_replicates <- 5

# What kind of figures do you want in your final outputs of MV distribution?
make_barplots <- TRUE
make_heatmaps <- TRUE
make_histograms <- TRUE

# Do you also want to have details on the intensities to ratios state plotted as histograms? 
# (combined or separated dependent on the class: Both, None, Numerator and Denominator)
# (only possible when proteinGroups and Phosphosites data are there)
make_combined_state_histograms <- TRUE
make_separate_state_histograms <- TRUE

# Do you want to make a consistence evaluation when classifying your intensities to ratios data?
make_consistency_eval <- TRUE

### ---End of parameters---

# First, set session or go to the directory where all the MaxQuant results are available (normally in combined/txt folder)

# From this point on, DO NOT MODIFY THE CODE!!!

# For conditions, create a txt/excel file named "overviewConditions.txt" where you indicate your groups/columns 
# where the rows are the replicates or experiments and the columns are the labels from SILAC experiments (see "overview" file as an example)
conditions_table_orig <- read.table("overviewConditions.txt", header=TRUE, sep="\t")

# If there are no replicates specified in the conditions table, copy each row and paste repeatedly with the given number of replicates
if (!have_replicates){
  n_rep <- nb_replicates
  # Repeat each row of the conditions table n_rep times
  conditions_table_orig <- conditions_table_orig[rep(row.names(conditions_table_orig), each=n_rep),]
  # Reset row names
  row.names(conditions_table_orig) <- NULL
}

# Remove the replicate/experiment column as it is not needed for futher usage
conditions_table <- conditions_table_orig[, -1]

# Get all the datasets of proteins and phosphosites
if (get_proteinGroups){
  proteinGroups <- read.table("proteinGroups.txt", header=TRUE, sep="\t", fill=TRUE, quote="")
}
if (get_Phosphosites){
  Phosphosites <- read.table("Phospho (STY)Sites.txt", header=TRUE, sep="\t", fill=TRUE, quote="")
}

# If you have the Perseus combined dataset (after normalization by log2), you can also get the dataset for combined data
# It should already be filtered and ridden of undesirable rows.
if (get_normedPhosphosites){
  normedPhosphosites <- read.table("Phospho (STY)Sites_normed.txt", header=TRUE, sep="\t", fill=TRUE, quote="")
}

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

if (get_proteinGroups){
  proteinGroups <- remove_rows_protein(proteinGroups)
}
if (get_Phosphosites){
  Phosphosites <- remove_rows_phos(Phosphosites)
}

# Labels bias

# Filter only normalized ratios as columns and only as separate experiments
if (get_proteinGroups){
  desired_columns_lab <- names(proteinGroups)[grep("^Ratio.[H|M].[M|L].normalized.[A-z]+[0-9]+$", names(proteinGroups))]
  labels_proteinGroups <- proteinGroups[,desired_columns_lab]
}
if (get_Phosphosites){
  desired_columns_lab <- names(Phosphosites)[grep("^Ratio.[H|M].[M|L].normalized.[A-z]+[0-9]+$", names(Phosphosites))]
  labels_Phosphosites <- Phosphosites[,desired_columns_lab]
}

# Rename them with a function to reduce the column names to its essential parts and as ratio naming
reduce_colnames <- function(df){
  colnames(df) <- gsub("normalized.", "", as.character(colnames(df)))
  colnames(df) <- gsub("Ratio.", "", as.character(colnames(df)))
  colnames(df) <- sub(".", "/", as.character(colnames(df)), fixed=TRUE)
}
if (get_proteinGroups){
  colnames(labels_proteinGroups) <- reduce_colnames(labels_proteinGroups)
}
if (get_Phosphosites){
  colnames(labels_Phosphosites) <- reduce_colnames(labels_Phosphosites)
}

# For the normed phosphosites, you would need to rename the columns for the function accessibility
if (get_normedPhosphosites){
  desired_columns_norm <- names(normedPhosphosites)[grep("^Ratio.[H|M].[M|L].normalized.[A-z]+[0-9]+_ph_x.y", names(normedPhosphosites))]
  labels_normedPhosphosites <- normedPhosphosites[,desired_columns_norm]
  colnames(labels_normedPhosphosites) <- desired_columns_lab
  colnames(labels_normedPhosphosites) <- reduce_colnames(labels_normedPhosphosites)
}

# Intensities bias

# Filter only intensities as columns and only as separate experiments
if (get_proteinGroups){
  desired_columns_int <- names(proteinGroups)[grep("^Intensity.[H|M|L].[A-z]+[0-9]+$", names(proteinGroups))]
  intensities_proteinGroups <- proteinGroups[,desired_columns_int]
}
if (get_Phosphosites){
  desired_columns_int <- names(Phosphosites)[grep("^Intensity.[H|M|L].[A-z]+[0-9]+$", names(Phosphosites))]
  intensities_Phosphosites <- Phosphosites[,desired_columns_int]
}

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
if (get_proteinGroups){
  conditions_by_intensity_proteinGroups <- rearrange_matrix_intensity(conditions_table, intensities_proteinGroups)
}
if (get_Phosphosites){
  conditions_by_intensity_Phosphosites <- rearrange_matrix_intensity(conditions_table, intensities_Phosphosites)
}
# Conditions bias based on labels (Replace the conditions with "conditions ratios" and rearrange them with labels)
library(stringr)

# Make correct labels based on the first row of the overview table

# Function to generate pair combinations for a given row
get_pair_combinations <- function(row_number, row_data){
  pairs_row <- combn(row_data, 2)
  sorted_pairs <- t(apply(pairs_row, 2, str_sort)) # Sort each pair alphabetically
  combined_pairs <- apply(sorted_pairs, 1, paste, collapse="/")
  complete_pairs <- paste0(combined_pairs,".Exp", row_number)
  return(complete_pairs)
}

# Function to rearrange the labels matrix with the conditions table in mind and possible inversions
rearrange_matrix_label <- function(conditions_matrix, numeric_matrix){
  # Initialize rearranged matrix that should be the same as the numeric matrix
  rearranged_matrix <- numeric_matrix
  
  # Pairing rearrangement part
  
  # Initialize a matrix to store pair combinations
  num_combinations <- choose(ncol(conditions_matrix),2)
  pair_combinations_matrix <- matrix("", nrow=nrow(conditions_matrix), ncol=num_combinations)
  
  for (i in 1:nrow(conditions_matrix)){
    row_pairs <- get_pair_combinations(i, conditions_matrix[i,])
    pair_combinations_matrix[i,] <- row_pairs
  }
  # Combine pair combinations into a single vector
  conditions_vector <- as.vector(pair_combinations_matrix)
  
  # Replace the column names from the rearranged matrix with the conditions ratio name vector
  colnames(rearranged_matrix) <- conditions_vector
  
  # Inversion of ratios part
  
  # Get the number of rows and columns from the conditions matrix
  num_conditions <- ncol(conditions_matrix)
  num_replicates <- nrow(conditions_matrix)
  
  # Initialize the ratio matrix
  ratio_matrix <- matrix(NA, nrow=num_replicates, ncol=num_conditions*(num_conditions-1) / 2)
  
  # Get all the possible ratios from the conditions matrix
  for (i in 1:num_replicates){
    ratios <- character(length = num_conditions)
    count <- 1
    for (j in 1:(num_conditions-1)){
      for (k in (j+1):num_conditions){
        ratios[count] <- paste0(conditions_matrix[i,j], "/", conditions_matrix[i,k])
        count <- count + 1
      }
    }
    ratio_matrix[i,] <- ratios
  }
  # Get the first row of the ratios matrix as basis for invert test
  base_ratios <- as.vector(ratio_matrix[1,])
  
  # Initialize an invert binary matrix to apply a possible inversion to the columns of the rearranged matrix
  invert_matrix <- matrix(NA, nrow=num_replicates, ncol=num_conditions)
  
  # Check at each value of the ratio matrix if an inversion should be applied or not based on its first row
  for (i in 1:num_replicates){
    for (j in 1:ncol(ratio_matrix)){
      if (ratio_matrix[i,j] %in% base_ratios){
        invert_matrix[i,j] <- FALSE
      }else{
        invert_matrix[i,j] <- TRUE
      }
    }
  }
  
  # Make the matrix into a vector and get column indexes where the value is set to TRUE
  invert_vector <- as.vector(invert_matrix)
  invert_indexes <- which(invert_vector)
  
  # If it is already normalized, apply another transformation where it's just a sign change
  if (deparse(substitute(numeric_matrix))=="labels_normedPhosphosites"){
    rearranged_matrix[,invert_indexes] <- (-1) * rearranged_matrix[,invert_indexes]
  }else{
    # Apply an inversion transformation to the columns if the invert vector is TRUE
    rearranged_matrix[,invert_indexes] <- 1 / rearranged_matrix[,invert_indexes]
  }
  # Reorder columns to have ordered conditions and replicates
  rearranged_matrix <- rearranged_matrix[, order(colnames(rearranged_matrix))]
  
  # Renaming columns part
  
  # Get the first row of the conditions table that would be the basis of label to condition assignment
  conditions_row <- conditions_matrix[1,]
  
  # Get the names of the columns to be renamed
  renamed_columns <- colnames(rearranged_matrix)
  
  # Iterate to each case that depends on the number of conditions (2-plex or 3-plex)
  for (l in 1:length(renamed_columns)){
    column <- renamed_columns[l]
    
    # Get the experiment string
    experiment <- unlist(str_extract_all(column, "(?<=\\.).*"))
    # Get the numerator string from the ratios
    numerator <- str_extract(column, "[^/]+")
    # Get the denominator string from the ratios (a bit more complicated)
    denominator <- unlist(str_extract_all(column, "(?<=/).*?(?=\\.)"))
    
    # Check for identities of those numerators and denominators in the first row of conditions
    num_index <- which(conditions_row == numerator)
    denom_index <- which(conditions_row == denominator)
    
    # Swap the ratios if denominator index is higher than numerator
    if (num_index < denom_index){
      correct_column <- paste0(denominator,"/",numerator,".",experiment)
    } else {
      correct_column <- column
    }
    
    # Rename the column appropriately
    renamed_columns[l] <- correct_column
  }
  
  # Generate pair combinations based on the first row's arrangement of the condition table
  first_row_pairs <- get_pair_combinations(1, conditions_matrix[1,])
  
  colnames(rearranged_matrix) <- renamed_columns
  
  return(rearranged_matrix)
}

if (get_proteinGroups){
  conditions_by_label_proteinGroups <- rearrange_matrix_label(conditions_table, labels_proteinGroups)
}
if (get_Phosphosites){
  conditions_by_label_Phosphosites <- rearrange_matrix_label(conditions_table, labels_Phosphosites)
}
if (get_normedPhosphosites){
  conditions_by_label_normedPhosphosites <-rearrange_matrix_label(conditions_table, labels_normedPhosphosites)
}

# Additional classification state hits from intensities to ratios (Phosphosites may take longer, around 5 minutes)

# Function to classify the ratios per column from intensities to labels (not the conditions)
classify_ratios_int_to_lab <- function(intensities_df, ratios_df){
  # Replace zeroes with NA
  intensities_df[intensities_df==0] <- NA
  ratios_df[ratios_df==0] <- NA
  
  # Check for missing values in both datasets (Binary df of TRUE/FALSE)
  missing_intensities <- is.na(intensities_df)
  missing_ratios <- is.na(ratios_df)
  
  # Initialize a dataframe to store classification results
  classification_df <- data.frame(matrix(nrow=nrow(ratios_df), ncol=ncol(ratios_df)))
  colnames(classification_df) <- colnames(ratios_df)
  rownames(classification_df) <- rownames(ratios_df)
  
  # Iterate through each cell in the ratio dataframe
  for (j in 1:ncol(ratios_df)){
    # Extract the column name from the ratios df
    column <- colnames(ratios_df)[j]
    numerator_col <- substr(column, 1, 1) # Extract numerator
    denominator_col <- substr(column, 3, 3) # Extract denominator
    experiment_number <- substr(column, 5, nchar(column))# Extract experiment number
    
    # Get the equivalent from the intensities data names
    numerator_name <- paste0("Intensity.",numerator_col,".",experiment_number)
    denominator_name <- paste0("Intensity.",denominator_col,".",experiment_number)
    
    for (i in 1:nrow(ratios_df)){
      # Check for MVs in numerator and denominator columns of intensities
      numerator_missing <- missing_intensities[i, numerator_name]
      denominator_missing <- missing_intensities[i, denominator_name]
      
      # Classify ratios based on presence of values (Both, none, denominator or numerator values present)
      if (!numerator_missing && !denominator_missing){
        classification_df[i,j] <- "Both"
      } else if (numerator_missing && denominator_missing){
        classification_df[i,j] <- "None"
      } else if (numerator_missing){
        classification_df[i,j] <- "Denom"
      } else {
        classification_df[i,j] <- "Num"
      }
    }
  }
  return(classification_df)
}

if (get_proteinGroups){
  classification_int_to_lab_proteinGroups <- classify_ratios_int_to_lab(intensities_proteinGroups, labels_proteinGroups)
}
if (get_Phosphosites){
  classification_int_to_lab_Phosphosites <- classify_ratios_int_to_lab(intensities_Phosphosites, labels_Phosphosites)
}

# Function to classify the ratios per column from conditions intensities to conditions labels
classify_ratios_cond_int_to_lab <- function(intensities_df, ratios_df){
  # Replace zeroes with NA
  intensities_df[intensities_df==0] <- NA
  ratios_df[ratios_df==0] <- NA
  
  # Check for missing values in both datasets (Binary df of TRUE/FALSE)
  missing_intensities <- is.na(intensities_df)
  missing_ratios <- is.na(ratios_df)
  
  # Initialize a dataframe to store classification results
  classification_df <- data.frame(matrix(nrow=nrow(ratios_df), ncol=ncol(ratios_df)))
  colnames(classification_df) <- colnames(ratios_df)
  rownames(classification_df) <- rownames(ratios_df)
  
  # Iterate through each cell in the ratio dataframe
  for (j in 1:ncol(ratios_df)){
    # Extract the column name from the ratios df
    column <- colnames(ratios_df)[j]
    numerator_col <- str_extract(column, "[^/]+") # Extract numerator
    denominator_col <- unlist(str_extract_all(column, "(?<=/).*?(?=\\.)")) # Extract denominator
    experiment_number <- unlist(str_extract_all(column, "(?<=\\.).*"))# Extract experiment number
    
    # Get the equivalent from the intensities data names
    numerator_name <- paste0(numerator_col,".",experiment_number)
    denominator_name <- paste0(denominator_col,".",experiment_number)
    
    for (i in 1:nrow(ratios_df)){
      # Check for MVs in numerator and denominator columns of intensities
      numerator_missing <- missing_intensities[i, numerator_name]
      denominator_missing <- missing_intensities[i, denominator_name]
      
      # Classify ratios based on presence of values (Both, none, denominator or numerator values present)
      if (!numerator_missing && !denominator_missing){
        classification_df[i,j] <- "Both"
      } else if (numerator_missing && denominator_missing){
        classification_df[i,j] <- "None"
      } else if (numerator_missing){
        classification_df[i,j] <- "Denom"
      } else {
        classification_df[i,j] <- "Num"
      }
    }
  }
  return(classification_df)
}

if (get_proteinGroups){
  classification_cond_int_to_lab_proteinGroups <- classify_ratios_cond_int_to_lab(conditions_by_intensity_proteinGroups, conditions_by_label_proteinGroups)
}
if (get_Phosphosites){
  classification_cond_int_to_lab_Phosphosites <- classify_ratios_cond_int_to_lab(conditions_by_intensity_Phosphosites, conditions_by_label_Phosphosites)
}
# Store the row ids from original data of the numerators or denominators that are present in those classifications (for verification)
if (get_proteinGroups){
  numerator_row_ids_intlab_proteinGroups <- names(which(rowSums("Num" == classification_int_to_lab_proteinGroups) > 0))
  denominator_row_ids_intlab_proteinGroups <- names(which(rowSums("Denom" == classification_int_to_lab_proteinGroups) > 0))
  
  numerator_row_ids_condintlab_proteinGroups <- names(which(rowSums("Num" == classification_cond_int_to_lab_proteinGroups) > 0))
  denominator_row_ids_condintlab_proteinGroups <- names(which(rowSums("Denom" == classification_cond_int_to_lab_proteinGroups) > 0))
}

if (get_Phosphosites){
  numerator_row_ids_intlab_Phosphosites <- names(which(rowSums("Num" == classification_int_to_lab_Phosphosites) > 0))
  denominator_row_ids_intlab_Phosphosites <- names(which(rowSums("Denom" == classification_int_to_lab_Phosphosites) > 0))
  
  numerator_row_ids_condintlab_Phosphosites <-names(which(rowSums("Num" == classification_cond_int_to_lab_Phosphosites) > 0))
  denominator_row_ids_condintlab_Phosphosites <- names(which(rowSums("Denom" == classification_cond_int_to_lab_Phosphosites) > 0))
}

# (Optional) Evaluation of consistency of those classifications (makes sure that all assignations of None, Num and Denom are only in NaN values)
consistence_evaluation <- function(original_df, classification_df){
  # Make sure that the two dataframes have the same length, othwerwise stop the evaluation
  if (!all(dim(original_df) == dim(classification_df))){
    stop("Error: The dimensions of the value and classified dataframes must match")
  }
  # Initiate an exception variable that will tell if there are inconsitencies in the dataframes
  exception <- FALSE
  # Initiate an inconsistencies vector that stores all the inconsistencies positions
  inconsistencies <- vector(mode="character")
  # Let's evaluate evey single value for consistency
  for (i in seq_len(nrow(original_df))){
    for (j in seq_len(ncol(original_df))){
      # In the case where the value is missing
      if (is.na(original_df[i,j])) {
        if(any(c("None", "Num", "Denom") == classification_df[i,j])){
          next
        } else {
          inconsistencies <- c(inconsistencies, paste0("(",i,",",j,")"))
          exception <- TRUE
        }
      }
      # In the case where the value is there (there should be no inconsistencies)
      if (!is.na(original_df[i,j])){
        if(classification_df[i,j] == "Both"){
          next
        } else {
          print(paste0("There are inconsistencies of full value from intensities to ratios in this position: (",i,",",j,")"))
          exception <- TRUE
        }
      }
    }
  }
  if (!exception){
    print("The classifications are consistent through the whole dataframe")
  }
  return(inconsistencies)
}

# Check if there are inconsistencies in all the possible dataframes
if (make_consistency_eval){
  if (get_proteinGroups){
    inconsistencies_lab_proteinGroups <- consistence_evaluation(labels_proteinGroups, classification_int_to_lab_proteinGroups)
  }
  if (get_Phosphosites){
    inconsistencies_lab_Phosphosites <- consistence_evaluation(labels_Phosphosites, classification_int_to_lab_Phosphosites)
  }
}
# The inconsistencies will be the same when looking at conditions (more on proteinGroups than Phosphosites, see remark on top of the script)
# Next step would be to visualize those states through histograms (at the end of this script)

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

# Run the following code depending on your parameter choices
if (get_proteinGroups){
  bins_by_labels_proteinGroups <- get_missing_value_bins_by_category(labels_proteinGroups)
  bins_by_intensities_proteinGroups <- get_missing_value_bins_by_category(intensities_proteinGroups)
  bins_by_conditions_labels_proteinGroups <- get_missing_value_bins_by_category(conditions_by_label_proteinGroups)
  bins_by_conditions_intensities_proteinGroups <- get_missing_value_bins_by_category(conditions_by_intensity_proteinGroups)
}

if (get_Phosphosites){
  bins_by_labels_Phosphosites <- get_missing_value_bins_by_category(labels_Phosphosites)
  bins_by_intensities_Phosphosites <- get_missing_value_bins_by_category(intensities_Phosphosites)
  bins_by_conditions_labels_Phosphosites <- get_missing_value_bins_by_category(conditions_by_label_Phosphosites)
  bins_by_conditions_intensities_Phosphosites <- get_missing_value_bins_by_category(conditions_by_intensity_Phosphosites)

}

if (get_normedPhosphosites){
  bins_by_labels_normedPhosphosites <- get_missing_value_bins_by_category(labels_normedPhosphosites)
  bins_by_conditions_labels_normedPhosphosites <- get_missing_value_bins_by_category(conditions_by_label_normedPhosphosites)
}

# Take notes of the MV distributions from your experiments for the generation of artificial matrices,
# there should not be many significant difference of proportions between the categories

# Save those MV distributions through txt files
if (get_proteinGroups){
  bins_overall_proteinGroups <- cbind(bins_by_labels_proteinGroups, bins_by_conditions_labels_proteinGroups, bins_by_intensities_proteinGroups, bins_by_conditions_intensities_proteinGroups)
  write.table(bins_overall_proteinGroups, file="proteinsMVdistribution.txt", sep="\t", col.names=NA)
}

if (get_Phosphosites){
  bins_overall_Phosphosites <- cbind(bins_by_labels_Phosphosites, bins_by_conditions_labels_Phosphosites, bins_by_intensities_Phosphosites, bins_by_conditions_intensities_Phosphosites)
  write.table(bins_overall_Phosphosites, file="phosphositesMVdistribution.txt", sep="\t", col.names=NA)
}

if (get_normedPhosphosites){
  bins_overall_normedPhosphosites <- cbind(bins_by_labels_normedPhosphosites, bins_by_conditions_labels_normedPhosphosites)
  write.table(bins_overall_normedPhosphosites, file="phosphosites_normedMVdistribution.txt", sep="\t", col.names=NA)
}

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

# Create barplots if wanted
if (make_barplots){
  if (get_proteinGroups){
    protein_barplots <- plot_barplots(bins_overall_proteinGroups)
  }
  if (get_Phosphosites){
    Phosphosites_barplots <- plot_barplots(bins_overall_Phosphosites)
  }
  if (get_normedPhosphosites){
    normedPhosphosites_barplots <- plot_barplots(bins_overall_normedPhosphosites)
  }
  
  # Save them in a pdf file (you could skip this if it does not work)
  pdf(file="Barplot_by_data_type_MV.pdf", width=11, height=8)
  if (get_proteinGroups){
    grid.arrange(grobs=protein_barplots, nrow=4, ncol=ncol(bins_by_intensities_proteinGroups), top="Overall Missing Value Distributions for Protein Groups")
  }
  if (get_Phosphosites){
    grid.arrange(grobs=Phosphosites_barplots, nrow=4, ncol=ncol(bins_by_intensities_Phosphosites), top="Overall Missing Value Distributions for Phosphosites (STY)")
  }
  if (get_normedPhosphosites){
    grid.arrange(grobs=normedPhosphosites_barplots, nrow=2, ncol=ncol(bins_by_labels_normedPhosphosites), top="Overall Missing Value Distributions for Normed ph - ft")
  }
  dev.off()
}

# Format to compare between types (first average out the bins by category matrix), all the matrices must be in a list first

# Function to compute the average of matrices
compute_matrix_average <- function(matrix_list){
  averages <- sapply(matrix_list, function(single_matrix){
    rowMeans(single_matrix, na.rm=TRUE)
  })
  return(averages)
}

# Function to compare a given list of matrices by average of rows/replicates
make_bins_list_average <- function(bins_list, prot){
  matrix_averages <- compute_matrix_average(bins_list)
  if (ncol(matrix_averages) == 1 & prot){
    colnames(matrix_averages) <- c("Protein Groups")
  } else if (ncol(matrix_averages) == 1 & !prot){
    colnames(matrix_averages) <- c("Phosphosites (STY)")
  } else if (ncol(matrix_averages) == 2){
    colnames(matrix_averages) <- c("Protein Groups", "Phosphosites (STY)")
  } else {
    colnames(matrix_averages) <- c("Protein Groups", "Phosphosites (STY)", "Normed ph - ft")
  }
  return(matrix_averages)
}

# Make a list of all the matrices per bias type and according to parameters and compute its average of matrices
if (make_barplots){
  have_prot <- TRUE
  if (get_proteinGroups & get_Phosphosites & get_normedPhosphosites){
    bins_by_labels_list <- list(bins_by_labels_proteinGroups, bins_by_labels_Phosphosites, bins_by_labels_normedPhosphosites)
    bins_by_intensities_list <- list(bins_by_intensities_proteinGroups, bins_by_intensities_Phosphosites)
    bins_by_conditions_labels_list <- list(bins_by_conditions_labels_proteinGroups, bins_by_conditions_labels_Phosphosites, bins_by_conditions_labels_normedPhosphosites)
    bins_by_conditions_intensities_list <- list(bins_by_conditions_intensities_proteinGroups, bins_by_conditions_intensities_Phosphosites)
  } else if (get_proteinGroups & get_Phosphosites){
    bins_by_labels_list <- list(bins_by_labels_proteinGroups, bins_by_labels_Phosphosites)
    bins_by_intensities_list <- list(bins_by_intensities_proteinGroups, bins_by_intensities_Phosphosites)
    bins_by_conditions_labels_list <- list(bins_by_conditions_labels_proteinGroups, bins_by_conditions_labels_Phosphosites)
    bins_by_conditions_intensities_list <- list(bins_by_conditions_intensities_proteinGroups, bins_by_conditions_intensities_Phosphosites)
  } else if (get_proteinGroups){
    bins_by_labels_list <- list(bins_by_labels_proteinGroups)
    bins_by_intensities_list <- list(bins_by_intensities_proteinGroups)
    bins_by_conditions_labels_list <- list(bins_by_conditions_labels_proteinGroups)
    bins_by_conditions_intensities_list <- list(bins_by_conditions_intensities_proteinGroups)
  } else if (get_Phosphosites){
    have_prot <- FALSE
    bins_by_labels_list <- list(bins_by_labels_Phosphosites)
    bins_by_intensities_list <- list(bins_by_intensities_Phosphosites)
    bins_by_conditions_labels_list <- list(bins_by_conditions_labels_Phosphosites)
    bins_by_conditions_intensities_list <- list(bins_by_conditions_intensities_Phosphosites)
  }
  
  matrix_averages_labels <- make_bins_list_average(bins_by_labels_list, have_prot)
  matrix_averages_intensites <- make_bins_list_average(bins_by_intensities_list, have_prot)
  matrix_averages_conditions_labels <- make_bins_list_average(bins_by_conditions_labels_list, have_prot)
  matrix_averages_conditions_intensities <- make_bins_list_average(bins_by_conditions_intensities_list, have_prot)
  
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
}

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
if (make_heatmaps){
  pdf(file="Heatmap_by_data_type_MV.pdf", width=13, height=8)
  if (get_proteinGroups){
    plot_heatmap(bins_overall_proteinGroups, "Protein Groups overall")
  }
  if (get_Phosphosites){
    plot_heatmap(bins_overall_Phosphosites, "Phosphosites (STY) overall")
  }
  if (get_normedPhosphosites){
    plot_heatmap(bins_overall_normedPhosphosites, "Normed ph - ft overall")
  }
  dev.off()
}

# Plot heatmaps per type of bias between datasets only when the barplots are also made
if (make_barplots & make_heatmaps){
  pdf(file="Heatmap_by_bias_type_MV.pdf", width=11, height=8)
  plot_heatmap(matrix_averages_labels, "labels average")
  plot_heatmap(matrix_averages_intensites, "intensities average")
  plot_heatmap(matrix_averages_conditions_labels, "conditions by labels average")
  plot_heatmap(matrix_averages_conditions_intensities, "conditions by intensities average")
  dev.off()
}

# Plot histograms per dataframe (once for global in blue, then with averaged missing values in red)
# Make a table of means, variance and skewness of those histograms too
library(parameters)

plot_histogram <- function(dataframes, prot, measure_type, pdf_name){
  # Start the saving of pdf files
  pdf(file=pdf_name, width=11, height=8)
  
  # Set up the dataframes index depending on the list's length for later iterations
  df_index <- 0
  
  # Initiate a dataframe table that summarizes basic metrics of histograms
  summary_table_df <- data.frame(Mean=numeric(), Variance=numeric(), Skewness=numeric(), stringsAsFactors=FALSE)
  
  # Set up the summary index as 1 initially (row indexing)
  summary_df_index <- 1
  
  for (df in dataframes){
    # Add 1 to the dataframes index
    df_index <- df_index + 1
    
    # Replace zeroes with NA
    df[df==0] <- NA
    
    # Extract categories from column names and adjust mfrow as such
    categories <- unique(gsub("\\d+$","", colnames(df)))
    par(mfrow=c(length(categories), 1), mar=rep(2,4))
    
    # Choose datatypes depending on the length of dataframes list
    if (length(dataframes) == 1 & have_prot){
      data_type <- c("Protein Groups")
    } else if (length(dataframes) == 1 & !have_prot) {
      data_type <- c("Phosphosites (STY)")
    } else if (length(dataframes) == 2){
      data_type <- c("Protein Groups", "Phosphosites (STY)")
    } else {
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
      
      # Save the metrics of those histograms in the summary table (mean, variance, skewness)
      summary_table_df[summary_df_index,] <- c(mean(mean_df, na.rm=TRUE), var(mean_df, na.rm=TRUE), skewness(mean_df)[1])
      row.names(summary_table_df)[summary_df_index] <- paste0("Global_", new_category, "_", data_type[df_index])
      
      summary_table_df[summary_df_index + 1,] <- c(mean(MV_df, na.rm=TRUE), var(MV_df, na.rm=TRUE), skewness(MV_df)[1])
      row.names(summary_table_df)[summary_df_index + 1] <- paste0("MV_", new_category, "_", data_type[df_index])
      
      # Add 2 to the index of the summary table
      summary_df_index <- summary_df_index + 2
    }
  }
  
  # Close the pdf device
  dev.off()
  
  return(summary_table_df)
}

if (make_histograms){
  have_prot <- TRUE
  # Create list of dataframes by labels, intensities and conditions
  if (get_proteinGroups & get_Phosphosites & get_normedPhosphosites){
    df_by_labels_list <- list(labels_proteinGroups, labels_Phosphosites, labels_normedPhosphosites)
    df_by_intensities_list <- list(intensities_proteinGroups, intensities_Phosphosites)
    df_by_conditions_labels_list <- list(conditions_by_label_proteinGroups, conditions_by_label_Phosphosites, conditions_by_label_normedPhosphosites)
    df_by_conditions_intensities_list <- list(conditions_by_intensity_proteinGroups, conditions_by_intensity_Phosphosites)
  } else if (get_proteinGroups & get_Phosphosites){
    df_by_labels_list <- list(labels_proteinGroups, labels_Phosphosites)
    df_by_intensities_list <- list(intensities_proteinGroups, intensities_Phosphosites)
    df_by_conditions_labels_list <- list(conditions_by_label_proteinGroups, conditions_by_label_Phosphosites)
    df_by_conditions_intensities_list <- list(conditions_by_intensity_proteinGroups, conditions_by_intensity_Phosphosites)
  } else if (get_proteinGroups){
    df_by_labels_list <- list(labels_proteinGroups)
    df_by_intensities_list <- list(intensities_proteinGroups)
    df_by_conditions_labels_list <- list(conditions_by_label_proteinGroups)
    df_by_conditions_intensities_list <- list(conditions_by_intensity_proteinGroups)
  } else if (get_Phosphosites){
    have_prot <- FALSE
    df_by_labels_list <- list(labels_Phosphosites)
    df_by_intensities_list <- list(intensities_Phosphosites)
    df_by_conditions_labels_list <- list(conditions_by_label_Phosphosites)
    df_by_conditions_intensities_list <- list(conditions_by_intensity_Phosphosites)
  }
  
  # Plot histograms for every possible dataset that is applicable and save the metrics in summary tables
  summary_hist_labels <- plot_histogram(df_by_labels_list, have_prot, "Ratio", "Histogram_by_labels_MV.pdf")
  summary_hist_intensities <- plot_histogram(df_by_intensities_list, have_prot, "", "Histogram_by_intensities_MV.pdf")
  summary_hist_cond_labels <- plot_histogram(df_by_conditions_labels_list, have_prot, "Ratio", "Histogram_by_conditions_labels_MV.pdf")
  summary_hist_cond_intensities <- plot_histogram(df_by_conditions_intensities_list, have_prot, "", "Histogram_by_conditions_intensities_MV.pdf")
  
  # Bind all those summary distributions per measure type and save them in a txt file
  summary_overall_hist_labels <- rbind(summary_hist_labels, summary_hist_cond_labels)
  summary_overall_hist_intensities <- rbind(summary_hist_intensities, summary_hist_cond_intensities)
  
  write.table(summary_overall_hist_labels, file="Summary_histogram_ratios.txt", sep="\t", col.names=NA)
  write.table(summary_overall_hist_intensities, file="Summary_histogram_intensities.txt", sep="\t", col.names=NA)
}

# Plot histograms per ratio dataframe by intensities to ratios presence state (Both, None, Num, Denom)

# Filtering function for classification assignment
filter_df_by_class <- function(values_df, classification_df, class_to_keep){
  # Check if both dataframes have the same dimension
  if (!all(dim(values_df) == dim(classification_df))){
    stop("Error: The dimensions of the value and classified dataframes must match")
  }
  
  # Initiate the filtered dataframe with the original values
  filtered_df <- values_df
  
  # Apply filter to each row
  for (row_data in seq_len(nrow(values_df))){
    # Check if any label in the row matches the label to keep
    if(!any(classification_df[row_data, ] == class_to_keep)){
      # If no match set the entire row to NaN
      filtered_df[row_data,] <- NaN
    }
  }
  
  return(filtered_df)
}

# Histogram plotting function (all together) with summary metrics
plot_histogram_state <- function(dataframes, classification_dataframes, pdf_name){
  # Start the saving of pdf files
  pdf(file=pdf_name, width=11, height=8)
  
  # Set up the dataframes index for later iterations and classification dataframes indexing
  df_index <- 0
  
  # Initiate a dataframe table that summarizes basic metrics of histograms
  summary_table_df <- data.frame(Mean=numeric(), Variance=numeric(), Skewness=numeric(), stringsAsFactors=FALSE)
  
  # Set up the summary index as 1 initially
  summary_df_index <- 1
  
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
      # Subset dataframe and classification dataframe for the current category
      cat_cols <- grep(paste0("^", category), colnames(df))
      cat_df <- df[, cat_cols, drop = FALSE]
      class_cat_df <- classification_dataframes[[df_index]][, cat_cols, drop=FALSE]
      
      # Rename the category without the ".Exp"
      new_category <- gsub("(.*)\\..*","\\1", category)
      
      # Log2 transform the data before anything
      log2_df <- sapply(cat_df, log2)
      
      # From the classification dataframe, build artificial dataframes that only contain values from mean_df if the given state is present
      Both_df <- filter_df_by_class(log2_df, class_cat_df, "Both")
      None_df <- filter_df_by_class(log2_df, class_cat_df, "None")
      Num_df <- filter_df_by_class(log2_df, class_cat_df, "Num")
      Denom_df <- filter_df_by_class(log2_df, class_cat_df, "Denom")
      
      # Create the means per row data from those filtered dataframes
      Both_means_df <- apply(Both_df, 1, function(x) mean(na.omit(x)))
      None_means_df <- apply(None_df, 1, function(x) mean(na.omit(x)))
      Num_means_df <- apply(Num_df, 1, function(x) mean(na.omit(x)))
      Denom_means_df <- apply(Denom_df, 1, function(x) mean(na.omit(x)))
      
      # Make dummy histograms that are not plotted to get the maximum frequencies of numerators and denominators
      hist_both_df <- hist(Both_means_df, breaks="Scott", plot=FALSE)
      hist_none_df <- hist(None_means_df, breaks=hist_both_df$breaks, plot=FALSE)
      hist_num_df <- hist(Num_means_df, breaks=hist_both_df$breaks, plot=FALSE)
      hist_denom_df <- hist(Denom_means_df, breaks=hist_both_df$breaks, plot=FALSE)
      
      max_frequency_num <- max(hist_num_df$counts)
      max_frequency_denom <- max(hist_denom_df$counts)
      maximum_limit <- max(c(10, max_frequency_num, max_frequency_denom))
      
      # Build histogram for the current category's global distribution (y range set for "zoom")
      c_both <- rgb(255,255,255, max=255, alpha=25, names="lt.white")
      if(data_type[df_index] == "Protein Groups"){
        hist_df <- hist(Both_means_df, breaks="Scott", col=c_both, ylim=c(0,10),
                        main=paste("Histogram for intensity presence states of", new_category, "for", data_type[df_index]), xlab=paste("Ratio", new_category, "measure"))
      } else if (data_type[df_index] == "Phosphosites (STY)"){
        hist_df <- hist(Both_means_df, breaks="Scott", col=c_both, ylim=c(0,maximum_limit),
                        main=paste("Histogram for intensity presence states of", new_category, "for", data_type[df_index]), xlab=paste("Ratio", new_category, "measure"))
      }
      
      # Build histogram in order for numerator and denominator presence with similar bins from the global distribution (transparent colors included)
      c_none <- rgb(255,0,0, max=255, alpha=25, names="lt.red")
      c_num <- rgb(0,0,255, max=255, alpha=125, names="lt.blue" )
      c_denom <- rgb(0,255,0, max=255, alpha=125, names="lt.green")
      hist(None_means_df, breaks=hist_df$breaks, col=c_none, add=TRUE)
      hist(Num_means_df, breaks=hist_df$breaks, col=c_num, add=TRUE)
      hist(Denom_means_df, breaks=hist_df$breaks, col=c_denom, add=TRUE)
      
      # Add a legend to indicate distribution difference
      legend("topright", legend=c("Global", "None", "Numerator only", "Denominator only"), fill=c(c_both, c_none, c_num, c_denom))
      
      # Save the metrics of those histograms in the summary table (mean, variance, skewness)
      summary_table_df[summary_df_index,] <- c(mean(Both_means_df, na.rm=TRUE), var(Both_means_df, na.rm=TRUE), skewness(Both_means_df)[1])
      row.names(summary_table_df)[summary_df_index] <- paste0("Global_", new_category, "_", data_type[df_index])
      
      summary_table_df[summary_df_index + 1,] <- c(mean(None_means_df, na.rm=TRUE), var(None_means_df, na.rm=TRUE), skewness(None_means_df)[1])
      row.names(summary_table_df)[summary_df_index + 1] <- paste0("None_", new_category, "_", data_type[df_index])
      
      summary_table_df[summary_df_index + 2,] <- c(mean(Num_means_df, na.rm=TRUE), var(Num_means_df, na.rm=TRUE), skewness(Num_means_df)[1])
      row.names(summary_table_df)[summary_df_index + 2] <- paste0("Num_", new_category, "_", data_type[df_index])
      
      summary_table_df[summary_df_index + 3,] <- c(mean(Denom_means_df, na.rm=TRUE), var(Denom_means_df, na.rm=TRUE), skewness(Denom_means_df)[1])
      row.names(summary_table_df)[summary_df_index + 3] <- paste0("Denom_", new_category, "_", data_type[df_index])
      
      # Add 4 to the index of the summary table
      summary_df_index <- summary_df_index + 4
      
      # Print the loop state in the console for the progress information
      print(paste("Histogram for", new_category, "in", data_type[df_index], "done"))
    }
  }
  
  # Close the pdf device
  dev.off()
  
  return(summary_table_df)
}

if (make_combined_state_histograms){
  # Create necessary variables for calling the histogram by state function
  df_by_labels_for_state_hist <- list(labels_proteinGroups, labels_Phosphosites)
  df_by_conditions_for_state_hist <- list(conditions_by_label_proteinGroups, conditions_by_label_Phosphosites)
  class_df_by_labels <- list(classification_int_to_lab_proteinGroups, classification_int_to_lab_Phosphosites)
  class_df_by_conditions <- list(classification_cond_int_to_lab_proteinGroups, classification_cond_int_to_lab_Phosphosites)
  
  # Plot histograms by state for every possible dataset that is applicable (5 minutes per call)
  summary_hist_state_labels <- plot_histogram_state(df_by_labels_for_state_hist, class_df_by_labels, "Histogram_labels_intensities_to_ratios_presence.pdf")
  summary_hist_state_conditions <- plot_histogram_state(df_by_conditions_for_state_hist, class_df_by_conditions, "Histogram_conditions_intensities_to_ratios_presence.pdf")
  
  # Bind all those summary distributions per and save them in a txt file
  summary_overall_hist_state_labels <- rbind(summary_hist_state_labels, summary_hist_state_conditions)
  write.table(summary_overall_hist_state_labels, file="Summary_histogram_state_intensities_to_ratios.txt", sep="\t", col.names=NA)
}

# Histogram plotting function (separated by classification type with global comparison)
plot_histogram_state_separated <- function(dataframes, classification_dataframes, pdf_name){
  # Start the saving of pdf files
  pdf(file=pdf_name, width=14, height=8)
  
  # Set up the dataframes index for later iterations and classification dataframes indexing
  df_index <- 0
  
  for (df in dataframes){
    # Add 1 to the dataframes index
    df_index <- df_index + 1
    
    # Replace zeroes with NA
    df[df==0] <- NA
    
    # Extract categories from column names and adjust mfrow as such (multiple columns this time)
    categories <- unique(gsub("\\d+$","", colnames(df)))
    par(mfrow=c(length(categories), 3), mar=rep(2,4))
    
    # Choose datatypes depending on the length of dataframes list
    if (length(dataframes) == 2){
      data_type <- c("Protein Groups", "Phosphosites (STY)")
    }
    else{
      data_type <- c("Protein Groups", "Phosphosites (STY)", "Normed ph - ft")
    }
    
    for (category in categories){
      # Subset dataframe and classification dataframe for the current category
      cat_cols <- grep(paste0("^", category), colnames(df))
      cat_df <- df[, cat_cols, drop = FALSE]
      class_cat_df <- classification_dataframes[[df_index]][, cat_cols, drop=FALSE]
      
      # Rename the category without the ".Exp"
      new_category <- gsub("(.*)\\..*","\\1", category)
      
      # Log2 transform the data before anything
      log2_df <- sapply(cat_df, log2)
      
      # From the classification dataframe, build artificial dataframes that only contain values from mean_df if the given state is present
      Both_df <- filter_df_by_class(log2_df, class_cat_df, "Both")
      None_df <- filter_df_by_class(log2_df, class_cat_df, "None")
      Num_df <- filter_df_by_class(log2_df, class_cat_df, "Num")
      Denom_df <- filter_df_by_class(log2_df, class_cat_df, "Denom")
      
      # Create the means per row data from those filtered dataframes
      Both_means_df <- apply(Both_df, 1, function(x) mean(na.omit(x)))
      None_means_df <- apply(None_df, 1, function(x) mean(na.omit(x)))
      Num_means_df <- apply(Num_df, 1, function(x) mean(na.omit(x)))
      Denom_means_df <- apply(Denom_df, 1, function(x) mean(na.omit(x)))
      
      # Make dummy histograms that are not plotted to get the maximum frequencies of numerators and denominators
      hist_both_df <- hist(Both_means_df, breaks="Scott", plot=FALSE)
      hist_none_df <- hist(None_means_df, breaks=hist_both_df$breaks, plot=FALSE)
      hist_num_df <- hist(Num_means_df, breaks=hist_both_df$breaks, plot=FALSE)
      hist_denom_df <- hist(Denom_means_df, breaks=hist_both_df$breaks, plot=FALSE)
      
      max_frequency_num <- max(hist_num_df$counts)
      max_frequency_denom <- max(hist_denom_df$counts)
      
      # Build seperated histograms for the current category's presence global distribution (ylim considered based on maximum frequency)
      hist_both_none_df <- hist(Both_means_df, breaks="Scott", col="yellow", border=FALSE,
                                main=paste("None of the intensities present for", new_category))
      hist(None_means_df, breaks=hist_both_none_df$breaks, col="red", add=TRUE)
      
      hist_both_num_df <- hist(Both_means_df, breaks="Scott", col="yellow", border=FALSE, ylim = c(0, max(c(10,max_frequency_num))),
                               main=paste("Only numerator intensities present for", new_category))
      hist(Num_means_df, breaks=hist_both_num_df$breaks, col="blue", add=TRUE)
      
      hist_both_denom_df <- hist(Both_means_df, breaks="Scott", col="yellow", border=FALSE, ylim = c(0,max(c(10,max_frequency_denom))),
                                 main=paste("Only denominator intensities present for", new_category))
      hist(Denom_means_df, breaks=hist_both_denom_df$breaks, col="green", add=TRUE)
      
      # Print the loop state in the console for the progress information
      print(paste("Histograms for", new_category, "in", data_type[df_index], "done"))
    }
    
    # Add a legend in the end to indicate the type of data we're working with and the labels
    legend("topright", legend = (c(data_type[df_index],"Global", "None", "Numerator", "Denominator")), fill=c("white", "yellow", "red", "blue", "green"))
  }
  
  # Close the pdf device
  dev.off()
}

if (make_separate_state_histograms){
  # Create necessary variables for calling the histogram by state function
  df_by_labels_for_state_hist <- list(labels_proteinGroups, labels_Phosphosites)
  df_by_conditions_for_state_hist <- list(conditions_by_label_proteinGroups, conditions_by_label_Phosphosites)
  class_df_by_labels <- list(classification_int_to_lab_proteinGroups, classification_int_to_lab_Phosphosites)
  class_df_by_conditions <- list(classification_cond_int_to_lab_proteinGroups, classification_cond_int_to_lab_Phosphosites)
  
  plot_histogram_state_separated(df_by_labels_for_state_hist, class_df_by_labels, "Histogram_separated_labels_intensities_to_ratios_presence.pdf")
  plot_histogram_state_separated(df_by_conditions_for_state_hist, class_df_by_conditions, "Histogram_separated_conditions_intensities_to_ratios_presence.pdf")
}