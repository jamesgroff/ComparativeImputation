# Notes: Make sure you generalize as much as possible all the different steps, remember you generate MVs for this script, do not impute!

# Notes for the 2nd option MV generation: the process is quite long, you may want to reduce or smart up your code for it to be faster
# For example, are there too many for loops? Do we really need to process row by row? Are there some steps that are unnecessary?
# It was faster when treating with binary dataframes than with full dataframes
# The generated intensities (not ratios) are not left distributed when looking at their histograms, but they have the right distributions
# (Definitely not for protein groups, but maybe for phosphosites?), how can we generate MNAR values?
# From what is observed in the original data of intensities, all the categories have the same MV attribution, meaning it is in the low values of intensities and MV placement is the same
# The thing we have in generation of MVs is that it respects the MV distribution, but it does not behave like a biased left gaussian distribution in intensities
# Just ignore it in the end and move on?

# Goal: Generate, based on a full data matrix (by category type without MV and split in two with no DE and only DE),
# several artificial matrices (20 times) that reflect the MV distribution found in the original data.
# Another way of generating these artificial matrices is to be based on missing intensities that were also generated according to its missing distribution.

# Required packages to install: stringr, dplyr, tidyr, ggplot2, gridExtra, reshape2

# Disclaimer: You should run DescriptionMV.R and BenchmarkStatistics.R before running this code script
# so that you get the distributions of MV in txt formats and the significant DE proteins/phosphosites as a dataframe basis

# First, set session or go to the directory where all the MaxQuant results are available (normally in combined/txt folder)

# General setup

# Set a general random seed for replicative purposes, or ignore it and continue for true randomness
set.seed(123)

# Get all the datasets of proteins, phosphosites (STY and normed) and their respective distribution matrices and significant dataframes
proteinGroups <- read.table("proteinGroups.txt", header=TRUE, sep="\t", fill=TRUE, quote="")
Phosphosites <- read.table("Phospho (STY)Sites.txt", header=TRUE, sep="\t", fill=TRUE, quote="")
normedPhosphosites <- read.table("Phospho (STY)Sites_normed.txt", header=TRUE, sep="\t", fill=TRUE, quote="")

distribMV_prot <- read.table("proteinsMVdistribution.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)
distribMV_phos <- read.table("phosphositesMVdistribution.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)
distribMV_norm <- read.table("phosphosites_normedMVdistribution.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)

significant_proteinGroups <- read.table("proteinGroups_significant.txt", header=TRUE, sep="\t")
significant_Phosphosites <- read.table("Phospho (STY)Sites_significant.txt", header=TRUE, sep="\t")
significant_normedPhosphosites <- read.table("Phospho (STY)Sites_normed_significant.txt", header=TRUE, sep="\t")

# Combine the significant proteins/phosphosites with the existing original dataframe
proteinGroups <- cbind(proteinGroups, significant_proteinGroups)
Phosphosites <- cbind(Phosphosites, significant_Phosphosites)
normedPhosphosites <- cbind(normedPhosphosites, significant_normedPhosphosites)

# We are only going to be interested into the condition label ratios and intensities for this project
# so pick only those ones as distribution and significance basis
desired_distribMV_columns <- names(distribMV_prot)[grep("^\\w{2,}/\\w{2,}.[A-z0-9]+$", names(distribMV_prot))]
condition_distribMV_prot <- distribMV_prot[,desired_distribMV_columns]
condition_distribMV_phos <- distribMV_phos[,desired_distribMV_columns]
condition_distribMV_norm <- distribMV_norm[,desired_distribMV_columns]

desired_distribMV_intensities <- names(distribMV_prot)[grep("^\\w{2,}.[A-z0-9]+$", names(distribMV_prot))]
intensities_distribMV_prot <- distribMV_prot[,desired_distribMV_intensities]
intensities_distribMV_phos <- distribMV_phos[,desired_distribMV_intensities]

desired_significant_columns <-  names(significant_proteinGroups)[grep("^Significant.\\w{2,}.\\w{2,}.[A-z0-9]+$", names(significant_proteinGroups))]
condition_significant_proteinGroups <- significant_proteinGroups[,desired_significant_columns]
condition_significant_Phosphosites <- significant_Phosphosites[,desired_significant_columns]
condition_significant_normedPhosphosites <- significant_normedPhosphosites[,desired_significant_columns]

# Remove rows that have potential contaminants and reversed sequences for full data and for significant data
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

condition_significant_proteinGroups <- proteinGroups[,desired_significant_columns]
condition_significant_Phosphosites <- Phosphosites[,desired_significant_columns]

# Substitute the significant columns so that they resemble to the actual names of distribution and full dataset
substitute_colnames <- function(columns_significant_df){
  substitute_columns <- sub("Significant.", "", columns_significant_df)
  substitute_columns <- sub("\\.", "/", substitute_columns)
  
  return(substitute_columns)
}
colnames(condition_significant_proteinGroups) <- substitute_colnames(desired_significant_columns)
colnames(condition_significant_Phosphosites) <- substitute_colnames(desired_significant_columns)
colnames(condition_significant_normedPhosphosites) <- substitute_colnames(desired_significant_columns)

# Filter only normalized ratios as columns and only as separate experiments for labels grouping
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

# Filter only intensities as columns and only as separate experiments
desired_columns_int <- names(proteinGroups)[grep("^Intensity.[H|M|L].[A-z]+[0-9]+$", names(proteinGroups))]
intensities_proteinGroups <- proteinGroups[,desired_columns_int]
intensities_Phosphosites <- Phosphosites[,desired_columns_int]

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

# Rename the columns of the conditions table to just initials depending if its 2-plex or 3-plex (L,M,H)
if (length(colnames(conditions_table))==2){
  colnames(conditions_table) <- c("L","H")
}
if (length(colnames(conditions_table))==3){
  colnames(conditions_table) <- c("L","M","H")
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

conditions_by_label_proteinGroups <- rearrange_matrix_label(conditions_table, labels_proteinGroups)
conditions_by_label_Phosphosites <- rearrange_matrix_label(conditions_table, labels_Phosphosites)
conditions_by_label_normedPhosphosites <- rearrange_matrix_label(conditions_table, labels_normedPhosphosites)

# 1st OPTION: Generate MVs according to separation of significance data and its original MV distribution

# Separate the data by conditions while removing MVs and if it is either DE or not DE by their significance, store them in a list

separate_matrices <- function(full_data, significant_data){
  # Extract column names indicating types and significance
  type_columns <- unique(gsub("\\d+$", "", names(full_data)))
  significance_columns <- colnames(significant_data)
    
  # Initialize a list to store matrices
  matrices_list <- list()
  
  # Iterate over type and significance combinations
  for (type_col in type_columns){
    for (significance_col in significance_columns){
      # Check if the significance column matches the type
      if (type_col == significance_col){
        # Extract indices of significant rows
        significant_rows <- which(significant_data[[type_col]]=="+")
        
        # Extract the full dataset for the current type
        type_data <- full_data[, grep(type_col, names(full_data))]
        
        # Filter significant rows from current type data
        significant_type_data <- type_data[significant_rows,]
        
        # Filter non-significant rows from current type data
        non_significant_type_data <- type_data[-significant_rows,]
        
        # Remove NA values into the two filtered data
        significant_type_data <- na.omit(significant_type_data)
        non_significant_type_data <- na.omit(non_significant_type_data)
        
        # Add the two filtered datasets into the matrices list
        matrices_list[[paste("DE", type_col, sep="_")]] <- significant_type_data
        matrices_list[[paste("Non-DE", type_col, sep="_")]] <- non_significant_type_data
      }
    }
  }
  return(matrices_list)
}

# Apply the function to all possible condition datasets
split_proteinGroups_list <- separate_matrices(conditions_by_label_proteinGroups, condition_significant_proteinGroups)
split_Phosphosites_list <- separate_matrices(conditions_by_label_Phosphosites, condition_significant_Phosphosites)
split_normedPhoshposites_list <- separate_matrices(conditions_by_label_normedPhosphosites, condition_significant_normedPhosphosites)

# Generate artificial missing values with similar MV distribution as original data

# Function to generate a certain number of matrices with full data and specified missing distribution in mind (for only one type)
generate_missing_matrices <- function(full_data, MV_distribution, num_matrices){
  # Create a list of artificial matrices that has both full and missing data
  matrices_with_missing <- list()
  for (i in 1:num_matrices){
    # Make a temporary copy of the full data matrix (to not write over the full data!)
    copy_full_data <- full_data
    
    # Get the number of columns and rows from the full data
    num_cols <- ncol(full_data)
    num_rows <- nrow(full_data)
    
    # For every generated matrix, create a missing binary dataset (1=MV, 0=no MV)
    missing_data <- matrix(0, nrow=num_rows, ncol=num_cols)
    
    for (row in 1:num_rows){
      # For each row of the matrix, decide the number of MVs to add depending on the given MV distribution
      num_missing <- sample(0:num_cols, size=1, prob=MV_distribution)
      # If there are missing values to assign, replace them accordingly and randomly
      if (num_missing > 0) {
        missing_indices <- sample(1:num_cols, num_missing)
        copy_full_data[row, missing_indices] <- NA
        missing_data[row, missing_indices] <- 1
      }
    }
    # Compile everything into the list of matrices
    matrices_with_missing[[i]] <- copy_full_data
  }
  
  return (matrices_with_missing)
}

# Function to pass over every matrices to generate artificial missing values from the given split list
apply_missing_matrices_to_list <- function(split_list, MV_distrib, num_matrices){
  # Create a list with the same structure as split_list
  matrices_list <- list()
  
  # Apply generate_missing_matrices function to every element in split_list
  for (i in 1:length(split_list)){
    # Get the element name
    element_name <- names(split_list)[i]
    # Replace the element name to have it similarly named to the MV distribution columns
    new_name <- gsub("DE_", "", gsub("Non-DE_", "", element_name))
    
    # Check if this new element name is in the MV distribution columns
    if (new_name %in% colnames(MV_distrib)){
      # Extract the MV_distribution type for the current element name
      MV_data <- MV_distrib[, grep(new_name, colnames(MV_distrib))]
      # Apply the generating of artificial MV matrices
      matrices_list[[element_name]] <- generate_missing_matrices(split_list[[i]], MV_data, num_matrices)
    }else{
      matrices_list[[element_name]] <- "Error: No Matches for this element"
    }
  }
  
  return(matrices_list)
}

# Specify the number of matrices you want to generate
N_matrices <- 20

# Call the main generator of MV matrices function as a test case for one DE/Non-DE matrix only
missing_matrices_example <- generate_missing_matrices(split_Phosphosites_list[[1]], condition_distribMV_phos[[1]], N_matrices)

# Function to verify the MV distribution from a list of generated matrices
check_distribution <- function(matrix_list){
  # Create a list to store the results for each matrix
  result_list <- list()
  
  for (i in 1:length(matrix_list)){
    # Access the i-th generated matrix
    matrix_with_missing <- matrix_list[[i]]
    # Calculate the distribution of MVs for verification on the first generated matrix
    missing_counts <- rowSums(is.na(matrix_with_missing))
    # Determine the maximum number of missing values
    max_missing <- max(missing_counts)
    # Create a matrix to store the counts
    bins <- matrix(0, nrow = max_missing + 1, ncol=1)
    # Count the frequencies of rows with different numbers of missing values
    for (j in 0:max_missing) {
      bins[j+1,1] <- sum(missing_counts==j)
    }
    rownames(bins) <- paste(0:max_missing, "MV", sep=" ")
    result_list[[i]] <- bins / nrow(matrix_with_missing)
  }
  # Convert the list into a proper matrix
  final_result <- do.call(cbind,result_list)
  
  return(final_result)
}
MVdistrib_missing_example <- check_distribution(missing_matrices_example)

# Get the average of those MV distributions
rowSums(MVdistrib_missing_example) / ncol(MVdistrib_missing_example)

# Functions to verify the integrity of those generated missing values
check_missing_rows <- function(matrix_data){
  apply(matrix_data, 1, function(row) all(is.na(row)))
}

check_integrity <- function(matrix_list){
  # Apply it to the example list across all matrices
  result <- sapply(matrix_list, check_missing_rows)
  # Check if any row has all missing values in all matrices
  if (any(rowSums(result) == length(matrix_list))){
    # If it does, then display those row indexes into a list
    matching_rows_indices <- list()
    for (i in 1:length(result)){
      matching_rows_indices[[i]] <- which(rowSums(result[i]) == length(matrix_list))
    }
    # Display matching rows
    for (i in 1:length(result)){
      cat("Matrix", i, "matching rows:", matching_rows_indices[[i]], "\n")
    }
  }else{
    cat("No rows match the case where all values in a row are missing across all matrices", "\n")
  }
  
}

# Optional: Repeat this integrity experiment example 10 times to see if there are chances to have one or multiple rows to have all missing values in all the generated matrices
#for (experiment_number in 1:10){
  #missing_matrices_example <- generate_missing_matrices(split_normedPhoshposites_list[[1]], condition_distribMV_norm[[1]], N_matrices)
  #check_integrity(missing_matrices_example)
#}

# Overall call if we do it by list of split matrices
missing_matrices_prot <- apply_missing_matrices_to_list(split_proteinGroups_list, condition_distribMV_prot, N_matrices)
missing_matrices_phos <- apply_missing_matrices_to_list(split_Phosphosites_list, condition_distribMV_phos, N_matrices)
missing_matrices_norm <- apply_missing_matrices_to_list(split_normedPhoshposites_list, condition_distribMV_norm, N_matrices)

# Function to plot histogram for a list with missing matrices
plot_histogram <- function(MV_matrices_list, data_type){
  # Get the length of the list
  list_length <- length(MV_matrices_list)
  
  # Adjust mfrow depending on the length of the list
  col_mfrow <- ceiling(sqrt(list_length))
  row_mfrow <- ceiling(list_length / col_mfrow)
  par(mfrow=c(row_mfrow,col_mfrow), mar=rep(2,4))
  
  for (i in 1:list_length){
    df <- MV_matrices_list[[i]]
    
    # Extract categories from column names
    categories <- unique(gsub("\\d+$","", colnames(df)))
    
    for (category in categories){
      # Subset dataframe for the current category
      cat_cols <- grep(paste0("^", category), colnames(df))
      cat_df <- df[, cat_cols, drop = FALSE]
      
      # Rename the category without the ".Exp"
      new_category <- gsub("(.*)\\..*","\\1", category)
      
      # Log2 transform the data before anything (except for normed data as it is already done)
      if (data_type != "Normed ph - ft"){
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
                      main=paste("Histogram n°", i,"of", new_category, "for", data_type), xlab=paste(new_category, "measure"))
      
      # Build histogram for missing values with the same bins from the global distribution
      hist(MV_df, breaks=hist_df$breaks, col="red", add=TRUE)
      
      # Add a legend to indicate distribution difference (if needed)
      #legend("topright", legend=c("Global", "MV only"), fill=c("blue", "red"))
    }
  }
}

# Check if it works for all your matrices with different datatypes as examples
plot_histogram(missing_matrices_prot[[2]], "Non-DE proteinGroups")
plot_histogram(missing_matrices_phos[[4]], "Non-DE Phosphosites")
plot_histogram(missing_matrices_norm[[6]], "Normed ph - ft")

# Combine both DE and non-DE matrices that were artificially generated with MV for later imputations
recombine_matrices <- function(missing_matrices_list){
  # Get the names of matrices
  matrix_names <- names(missing_matrices_list)
  
  # Extract the number of generated matrices per list element
  length_generated_matrices <- length(missing_matrices_list[[matrix_names[1]]])
  
  # Extract the unique types and significance levels
  unique_types <- unique(sub(".*?_","", matrix_names))
  unique_significances <- unique(gsub("_.*","", matrix_names))
  
  # Initialize a list to store those combined matrices
  recombined_matrices <- list()
  
  # Iterate over the number of generated matrices
  for (i in 1:length_generated_matrices){
    # Iterate over unique types
    for (type in unique_types){
      # Initialize an empty matrix to store recombined data for the current type
      combined_data <- matrix(nrow = 0, ncol = 0)
      # Iterate over unique significances
      for (significance in unique_significances){
        # Get the matrix corresponding to the current type and significance
        matrix_name <- paste0(significance,"_",type)
        matrix_data <- missing_matrices_list[[matrix_name]][[i]]
        
        # Add the matrix data to the combined data matrix
        combined_data <- rbind(combined_data, matrix_data)
      }
      # Add the combined data matrix to the list of recombined matrices
      recombined_matrices[[type]][[i]] <- combined_data
      
      # Reorder the current combined data by row names for clarity
      copy_df <- recombined_matrices[[type]][[i]]
      recombined_matrices[[type]][[i]] <- copy_df[order(as.numeric(row.names(copy_df))),]
    }
  }
  return(recombined_matrices)
}

recombined_missing_matrices_prot <- recombine_matrices(missing_matrices_prot)
recombined_missing_matrices_phos <- recombine_matrices(missing_matrices_phos)
recombined_missing_matrices_norm <- recombine_matrices(missing_matrices_norm)

# Check the MV distributions on each of those lists
MVdistrib_missing_prot <- check_distribution(recombined_missing_matrices_prot[[1]])
MVdistrib_missing_phos <- check_distribution(recombined_missing_matrices_phos[[1]])
MVdistrib_missing_norm <- check_distribution(recombined_missing_matrices_norm[[1]])

# Check once again their histograms
plot_histogram(recombined_missing_matrices_prot[[1]], "proteinGroups")
plot_histogram(recombined_missing_matrices_phos[[1]], "Phosphosites")
plot_histogram(recombined_missing_matrices_norm[[1]], "Normed ph - ft")

# Save those matrices in txt formats in different folders per datatype (all categories are in one same folder)
save_matrices <- function(folder_name, matrix_list){
  # Delete the folder if it already exists
  if(file.exists(folder_name)){
    unlink(folder_name, recursive=TRUE)
  }
  
  # Create a new folder
  if(!file.exists(folder_name)){
    dir.create(folder_name)
  }
  
  # Get the names of matrices
  matrix_names <- names(matrix_list)
  
  # Iterate through the nested list and save each matrix as a separate file
  for (name in matrix_names){
    current_list <- matrix_list[[name]]
    # Sanitize this name for file naming (because of "/")
    sanitized_name <- gsub("/","_over_", name)
    for (i in seq_along(current_list)){
      file_name <- paste0(folder_name, "/", sanitized_name, "_matrix_n°", i, ".txt")
      write.table(matrix_list[[name]][[i]], file=file_name, sep="\t", col.names=NA)
    }
  }
}

save_matrices("Generated_MV_matrices_proteinGroups", recombined_missing_matrices_prot)
save_matrices("Generated_MV_matrices_Phosphosites", recombined_missing_matrices_phos)
save_matrices("Generated_MV_matrices_normedPhoshposites", recombined_missing_matrices_norm)

# 2nd OPTION: Generate MVs according to the generated missing intensities into the ratios

# Make the condition intensities into full data by removing all rows that have at least one MV
make_intensities_full_data <- function(intensities_data){
  # Replace zeroes with NA
  intensities_data[intensities_data==0] <- NA
  
  # Remove those NA values to get full data
  full_data <- na.omit(intensities_data)
  
  return(full_data)
}

full_conditions_by_intensity_proteinGroups <- make_intensities_full_data(conditions_by_intensity_proteinGroups)
full_conditions_by_intensity_Phosphosites <- make_intensities_full_data(conditions_by_intensity_Phosphosites)

# Load the intensities histogram summary metrics
summary_hist <- read.table("Summary_histogram_intensities.txt", header=TRUE, sep="\t", row.names=1)

# Create artificial missing condition intensities based on its full data and its MV distribution
library(dplyr)
library(tidyr)

# Function to generate artificial missing intensities to every matrices from the given dataframe
apply_missing_matrices_to_dataframe <- function(full_dataframe, MV_distrib, num_matrices){
  # Create a list for the generated missing matrices storage
  matrices_list <- list()
  
  # Extract categories from column names of the dataframe
  categories <- unique(gsub("\\d+$","", colnames(full_dataframe)))
  
  for (category in categories){
    # Subset dataframe for the current category
    cat_cols <- grep(paste0("^", category), colnames(full_dataframe))
    cat_df <- full_dataframe[, cat_cols, drop = FALSE]
    
    # Check if this new category name is in the MV distribution columns
    if (category %in% colnames(MV_distrib)){
      # Extract the MV_distribution type for the current element name
      MV_data <- MV_distrib[, grep(category, colnames(MV_distrib))]
      # Apply the generating of artificial MV matrices
      matrices_list[[category]] <- generate_missing_matrices(cat_df, MV_data, num_matrices)
    }else{
      matrices_list[[category]] <- "Error: No Matches for this element"
    }
  }
  
  return(matrices_list)
}

apply_missing_matrices_to_dataframe <- function(full_dataframe, MV_distrib, summary_histogram, data_type, num_matrices){
  # Create a list for the generated missing matrices storage
  matrices_list <- list()
  
  # Extract categories from column names of the dataframe
  categories <- unique(gsub("\\d+$","", colnames(full_dataframe)))
  
  for (category in categories){
    # Subset dataframe for the current category
    cat_cols <- grep(paste0("^", category), colnames(full_dataframe))
    cat_df <- full_dataframe[, cat_cols, drop = FALSE]
    
    # Get the number of rows and columns from the dataset
    num_rows <- nrow(cat_df)
    num_cols <- ncol(cat_df)
    
    # Check if this new category name is in the MV distribution columns
    if (category %in% colnames(MV_distrib)){
      # log2 transform the dataset
      log2_df <- log2(cat_df)
      
      # Calculate the mean of each row
      row_means <- rowMeans(log2_df)
      
      # Rename the category without the ".Exp"
      new_category <- gsub("(.*)\\..*","\\1", category)
      
      # Extract the MV distribution type for the current element name
      MV_data <- MV_distrib[, grep(category, colnames(MV_distrib))]
      
      # Generate the number of MV per row
      num_missing_per_row <- sample(0:(length(MV_data)-1), num_rows, replace=TRUE, prob=MV_data)
      
      # Calculate the expected number of MVs on the data based on this MV distribution
      total_missing_values <- sum(num_missing_per_row)
      
      # Get the summary that are based on data type, conditions and MV summary histograms (not global)
      summary_hist_cat <- summary_histogram[grep(paste0("MV_",new_category,"_",data_type), rownames(summary_histogram)),]
      
      # Generate a normal distribution based on the data given from the summary histogram
      biased_mean <- summary_hist_cat$Mean
      biased_sd <- sqrt(summary_hist_cat$Variance)
      MV_normal_distrib <- rnorm(total_missing_values, mean=biased_mean, sd=biased_sd)
      
      # Ensure the generated biased values are within the range of the dataset
      biased_values <- pmax(pmin(MV_normal_distrib, max(log2_df)), min(log2_df))
      
      for (i in 1:num_matrices){
        # Make a temporary copy of the full data matrix (to not write over the full data!)
        copy_full_data <- log2_df
        
        # Initiate a count to create missing values
        current_missing_count <- 1
        
        # Recreate MVs with the desired distribution
        for (r in 1:num_rows){
          if (num_missing_per_row[r] > 0){
            for(j in 1:num_missing_per_row[r]){
              # Ensure we do not exceed the number of generated MVs
              if(current_missing_count <= total_missing_values){
                # Find the column index in the current row that is closest to the biased missing value
                closest_col <- which.min(abs(copy_full_data[r,] - biased_values[current_missing_count]))
                copy_full_data[r, closest_col] <- NA
                current_missing_count <- current_missing_count + 1
              }
            }
          }
        }
        
        # Get the generated matrix into the list while converting it back to the original scale
        matrices_list[[category]][[i]] <- 2^copy_full_data
      }
    }else{
      matrices_list[[category]][[i]] <- "Error: No Matches for this element"
    }
  }
  
  return(matrices_list)
}

N_matrices <- 20
missing_intensities_prot <- apply_missing_matrices_to_dataframe(full_conditions_by_intensity_proteinGroups, intensities_distribMV_prot, summary_hist, "Protein Groups", N_matrices)
missing_intensities_phos <- apply_missing_matrices_to_dataframe(full_conditions_by_intensity_Phosphosites, intensities_distribMV_phos, summary_hist, "Phosphosites (STY)", N_matrices)

# Verify that those matrices are distributed more on the left compared to the global distribution
plot_histogram(missing_intensities_prot[[1]], "proteinGroups")

# Join those separated categories into full intensities data summary just like in the original data format
join_generated_matrices <- function(separated_matrices_list){
  # Create a list for the joined matrices
  joined_matrices <- list()
  
  for (j in 1:length(separated_matrices_list[[1]])){
    for (i in 1:length(separated_matrices_list)){
      if(i==1){
        # If we are starting with a new matrix, initiate a matrix with the current data
        current_matrix <- separated_matrices_list[[i]][[j]]
      } else {
        # We bind the columns of the current matrix and the current data if it's not the first call
        current_matrix <- cbind(current_matrix, separated_matrices_list[[i]][[j]])
      }
    }
    joined_matrices[[j]] <- current_matrix
  }
  return(joined_matrices)
}

missing_joined_intensities_prot <- join_generated_matrices(missing_intensities_prot)
missing_joined_intensities_phos <- join_generated_matrices(missing_intensities_phos)

# Get the full data for condition ratios
full_conditions_by_label_proteinGroups <- na.omit(conditions_by_label_proteinGroups)
full_conditions_by_label_Phosphosites <- na.omit(conditions_by_label_Phosphosites)

# Only 1 generated MV matrix from intensities to ratios version (testing purposes)
generate_missing_ratios_test <- function(full_ratio_data, missing_intensities_list){
  # Get the current missing matrix from the intensities list
  missing_intensity_matrix <- missing_intensities_list[[1]]
  
  # Get the common list of rows that are the same between the full ratio data and the current missing intensities matrix
  common_list <- base::intersect(rownames(full_ratio_data), rownames(missing_intensity_matrix))
  
  # Make a temporary copy of the full data matrix (to not write over the full data!)
  copy_full_data <- full_ratio_data
  
  # Check for missing values in the intensities dataset (Binary df of TRUE/FALSE)
  missing_intensities <- is.na(missing_intensity_matrix)
  
  # Initiate a MV dataframe of ratios that has the same dimension as the original
  missing_ratios <- is.na(full_ratio_data)
  
  # Extract categories from column names of the dataframe
  categories <- unique(gsub("\\d+$","", colnames(full_ratio_data)))
  
  # Iterate by commmon id list
  for (row in common_list){
    # Check if the current row of the intensity matrix has a missing value, skip the row if not
    if (any(is.na(missing_intensity_matrix[row,]))){
      # If there are missing values to assign, add TRUE to the binary dataframe accordingly by ratio category cases
      
      for (category in categories){
        # Subset dataframe for the current category
        cat_cols <- grep(paste0("^", category), colnames(full_ratio_data))
        
        # Initiate a binary row that will be transplanted to the final binary dataframe
        binary_row <- missing_ratios[row, cat_cols, drop = FALSE]
        
        # Get the numerator and denominator name of this category ratio
        numerator_col <- str_extract(category, "[^/]+") # Extract numerator
        denominator_col <- unlist(str_extract_all(category, "(?<=/).*?(?=\\.)")) # Extract denominator
        experiment_name <- unlist(str_extract_all(category, "(?<=\\.).*")) # Extract experiment name without the number
        
        # Get the equivalent names from the intensities data column names
        numerator_name <- paste0(numerator_col,".",experiment_name)
        denominator_name <- paste0(denominator_col,".",experiment_name)
        
        # Get the indexes of those column names for either the numerator or denominator
        col_index_num <- grep(numerator_name, colnames(missing_intensities))
        col_index_denom <- grep(denominator_name, colnames(missing_intensities))
        
        # Check if this new category name's numerator or denominator is in the columns of the missing intensities matrix
        if (length(col_index_num) > 0 || length(col_index_denom) > 0){
          # Extract the missing row data for the current common row in either numerator or denominator
          MV_numerator_data <- missing_intensities[row, col_index_num]
          MV_denominator_data <- missing_intensities[row, col_index_denom]
          
          # Find columns with NAs in each row data
          na_cols_num <- which(MV_numerator_data > 0)
          na_cols_denom <- which(MV_denominator_data > 0)
          
          # Combine and get unique column indices that have MVs (only when both, denominator or numerator values are absent)
          missing_indices <- unique(c(na_cols_num, na_cols_denom))
          
          # Replace with TRUE to the current binary row based on presence of values
          binary_row[,missing_indices] <- TRUE
          missing_ratios[row, grep(category, colnames(missing_ratios))] <- binary_row
        }
      }
    }
  }
  # Transfer NA values to the full data using the binary dataframe
  copy_full_data[missing_ratios > 0] <- NA
  
  return (copy_full_data)
}

example_missing_int_to_ratios_prot <- generate_missing_ratios_test(full_conditions_by_label_proteinGroups, missing_joined_intensities_prot)
example_missing_int_to_ratios_phos <- generate_missing_ratios_test(full_conditions_by_label_Phosphosites, missing_joined_intensities_phos)

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

# Test it on the example missing intensities matrix
MV_bins_example_prot <- get_missing_value_bins_by_category(example_missing_int_to_ratios_prot)
MV_bins_example_phos <- get_missing_value_bins_by_category(example_missing_int_to_ratios_phos)

# The result is different from the original distribution of ratios!
# Even when the intensities distribution of MVs are the same!

# Heatmap visualization for comparison of distributions between original and one of the generated MVs
library(ggplot2)
library(gridExtra)
library(reshape2)

plot_heatmap <- function(original_distrib, generated_distrib, main_title){
  # Convert the original distribution into a matrix for melt function to work
  original_distrib <- as.matrix(original_distrib)
  
  # Rename the columns of the two distribution so that they're different enough in text format
  colnames(original_distrib) <- paste0("original_", colnames(original_distrib))
  colnames(generated_distrib) <- paste0("generated_", colnames(generated_distrib))
  
  # Bind those two distributions data into one
  matrix_data <- cbind(original_distrib, generated_distrib)
  
  # Reshape the matrix into a long format
  missing_data <- melt(matrix_data)
  
  # Plot the heatmap of the matrix
  ggplot(missing_data, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile() +
    geom_text(aes(label=format(round(value, 4), nsmall=4)), color="white", size=5) +
    scale_fill_gradient(low="lightblue", high="darkblue") +
    guides(fill=guide_colorbar(title="Proportion", barwidth=0.5, barheight=20)) +
    labs(x="Data type", y="Number of missing values", title=paste0("Heatmap of MV Distribution from intensities to ratios for ", main_title, ", MaxQuant vs Manual comparison")) +
    theme_minimal()
}

# Save those heatmaps in pdf format
pdf(file="Heatmap_comparison_distribution_original_vs_generated_MV.pdf", width=11, height=8)
plot_heatmap(condition_distribMV_prot, MV_bins_example_prot, "Protein Groups")
plot_heatmap(condition_distribMV_phos, MV_bins_example_phos, "Phosphosites (STY)")
dev.off()

# Function to generate MVs in ratios based on the missing intensities for list format
generate_missing_ratios_from_intensities <- function(full_ratio_data, missing_intensities_list){
  # Create a list of artificial matrices that will have missing data
  matrices_with_missing <- list()

  for (i in 1:length(missing_intensities_list)){
    # Time the process per matrix
    start_time <- Sys.time()
    
    # Get the current missing matrix from the intensities list
    missing_intensity_matrix <- missing_intensities_list[[i]]
    
    # Get the common list of rows that are the same between the full ratio data and the current missing intensities matrix
    common_list <- base::intersect(rownames(full_ratio_data), rownames(missing_intensity_matrix))
    
    # Make a temporary copy of the full data matrix (to not write over the full data!)
    copy_full_data <- full_ratio_data
    
    # Check for missing values in the intensities dataset (Binary df of TRUE/FALSE)
    missing_intensities <- is.na(missing_intensity_matrix)
    
    # Initiate a MV dataframe of ratios that has the same dimension as the original
    missing_ratios <- is.na(full_ratio_data)
    
    # Extract categories from column names of the dataframe
    categories <- unique(gsub("\\d+$","", colnames(full_ratio_data)))
    
    # Iterate by commmon id list
    for (row in common_list){
      # Check if the current row of the intensity matrix has a missing value, skip the row if not
      if (any(is.na(missing_intensity_matrix[row,]))){
        # If there are missing values to assign, add TRUE to the binary dataframe accordingly by ratio category cases
        
        for (category in categories){
          # Subset dataframe for the current category
          cat_cols <- grep(paste0("^", category), colnames(full_ratio_data))
          
          # Initiate a binary row that will be transplanted to the final binary dataframe
          binary_row <- missing_ratios[row, cat_cols, drop = FALSE]
          
          # Get the numerator and denominator name of this category ratio
          numerator_col <- str_extract(category, "[^/]+") # Extract numerator
          denominator_col <- unlist(str_extract_all(category, "(?<=/).*?(?=\\.)")) # Extract denominator
          experiment_name <- unlist(str_extract_all(category, "(?<=\\.).*")) # Extract experiment name without the number
          
          # Get the equivalent names from the intensities data column names
          numerator_name <- paste0(numerator_col,".",experiment_name)
          denominator_name <- paste0(denominator_col,".",experiment_name)
          
          # Get the indexes of those column names for either the numerator or denominator
          col_index_num <- grep(numerator_name, colnames(missing_intensities))
          col_index_denom <- grep(denominator_name, colnames(missing_intensities))
          
          # Check if this new category name's numerator or denominator is in the columns of the missing intensities matrix
          if (length(col_index_num) > 0 || length(col_index_denom) > 0){
            # Extract the missing row data for the current common row in either numerator or denominator
            MV_numerator_data <- missing_intensities[row, col_index_num]
            MV_denominator_data <- missing_intensities[row, col_index_denom]
            
            # Find columns with NAs in each row data
            na_cols_num <- which(MV_numerator_data > 0)
            na_cols_denom <- which(MV_denominator_data > 0)
            
            # Combine and get unique column indices that have MVs (only when both, denominator or numerator values are absent)
            missing_indices <- unique(c(na_cols_num, na_cols_denom))
            
            # Replace with TRUE to the current binary row based on presence of values
            binary_row[,missing_indices] <- TRUE
            missing_ratios[row, grep(category, colnames(missing_ratios))] <- binary_row
          }
        }
      }
    }
    # Transfer NA values to the full data using the binary dataframe
    copy_full_data[missing_ratios > 0] <- NA
    
    # Compile everything into the list of matrices
    matrices_with_missing[[i]] <- copy_full_data
    
    # Stop the time process here
    end_time <- Sys.time()
    time_taken <- end_time - start_time
    
    # Print the progress of the process into the console matrix by matrix
    print(paste0("Matrix n°", i, " done in ", time_taken, " minutes"))
  }
  
  return (matrices_with_missing)
}

# Long process (1h30 in total), so do it either on the server or when you have time to kill
missing_int_to_ratios_prot <- generate_missing_ratios_from_intensities(full_conditions_by_label_proteinGroups, missing_joined_intensities_prot)
missing_int_to_ratios_phos <- generate_missing_ratios_from_intensities(full_conditions_by_label_Phosphosites, missing_joined_intensities_phos)

# Plot those missing ratios in histogram format to compare global and MV distribution
plot_histogram(missing_int_to_ratios_prot, "proteinGroups")
plot_histogram(missing_int_to_ratios_phos, "Phosphosites")

# Save those matrices in other folders than the first option
save_matrices_individual <- function(folder_name, matrix_list, data_name){
  # Delete the folder if it already exists
  if(file.exists(folder_name)){
    unlink(folder_name, recursive=TRUE)
  }
  
  # Create a new folder
  if(!file.exists(folder_name)){
    dir.create(folder_name)
  }
  
  # Iterate through the nested list and save each matrix as a separate file
  for (i in seq_along(matrix_list)){
    current_matrix <- matrix_list[[i]]
    file_name <- paste0(folder_name, "/", "Generated_matrix_n°", i, "_", data_name, ".txt")
    write.table(current_matrix, file=file_name, sep="\t", col.names=NA)
  }
}

# Do it for the generated missing intensities and the missing ratios from intensities for correction and later steps
save_matrices_individual("Generated_MV_intensities_proteinGroups", missing_joined_intensities_prot, "intensities_proteinGroups")
save_matrices_individual("Generated_MV_intensities_Phosphosites", missing_joined_intensities_phos, "intensities_Phosphosites")

save_matrices_individual("Generated_MV_intensities_to_ratios_proteinGroups", missing_int_to_ratios_prot, "int_to_ratio_proteinGroups")
save_matrices_individual("Generated_MV_intensities_to_ratios_Phosphosites", missing_int_to_ratios_phos, "int_to_ratio_Phosphosites")

# Next step would be to impute those generated data based on intensities MV with either existing strategies or new ones (???)