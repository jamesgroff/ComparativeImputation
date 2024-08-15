# Notes: Check if imputations are random values or not. (minmax, gaussian and mle will be random unless seed is added, knn and rf will have the same values)
# There still could be boundary failures when doing the identification of DE proteins/phosphosites process (fix around the warnings with nlme package, see mail)
# There are lapply warnings during the identification of DE proteins/phosphosites (in gaussian and mle DE identification mostly)
# REML=FALSE when calling the linear model to fix problems (Rohr's advice in mail)

# There is a problem with IB_01 generated missing data in phosphosites with Exp1 having abnormal high values (origin: value 13095000 at row 5590 from original phos data)

# The minmax imputed distribution not looking like the uniform distribution in lower ratio region is because of the log2 transformation when plotting the means histogram (and minmax range by column)

# Goal: Impute the missing values from the generated matrices stored in folders with different methods:
# Negative control with random min-max range and Gaussian distribution sampling, then impute with kNN, MLE and RandomForest
# Then make confusion matrices and average them out for method comparisons with metrics such as FDR, FNR, accuracy and Matthews correlation coefficient (NRMSE too as an option)

# Required packages to install: stringr, mice, impute, imp4p, missForest, lme4, pracma, caret

# Disclaimer: You should run GenerateMatrixMV.R before running this code script so that you have all the necessary folders

# First, set session or go to the directory where all the MaxQuant results are available (normally in combined/txt folder)

# General setup

# Set a general random seed for replicative purposes, or ignore it and continue for true randomness
set.seed(123)

# Get all the original datasets of proteins and phosphosites (STY and normed) and their respective significant dataframes
proteinGroups <- read.table("proteinGroups.txt", header=TRUE, sep="\t", fill=TRUE, quote="")
Phosphosites <- read.table("Phospho (STY)Sites.txt", header=TRUE, sep="\t", fill=TRUE, quote="")
normedPhosphosites <- read.table("Phospho (STY)Sites_normed.txt", header=TRUE, sep="\t", fill=TRUE, quote="")

significant_proteinGroups <- read.table("proteinGroups_significant.txt", header=TRUE, sep="\t")
significant_Phosphosites <- read.table("Phospho (STY)Sites_significant.txt", header=TRUE, sep="\t")
significant_normedPhosphosites <- read.table("Phospho (STY)Sites_normed_significant.txt", header=TRUE, sep="\t")

# Combine the significant proteins/phosphosites with the existing original dataframe
proteinGroups <- cbind(proteinGroups, significant_proteinGroups)
Phosphosites <- cbind(Phosphosites, significant_Phosphosites)
normedPhosphosites <- cbind(normedPhosphosites, significant_normedPhosphosites)

# We are only going to be interested into the condition label ratios for this project, so pick only those ones as significance basis

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

# Getting the original matrix list

# Labels bias

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

# Conditions bias

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

# Careful, the rearrangements between generated matrices and original matrices are different in their categories

get_original_matrix_list <- function(complete_df) {
  # Get the unique categories from the complete dataset
  categories <- unique(gsub("\\d+$","", colnames(complete_df)))
  
  # Create a list to store the results for each category
  result_list <- list()
  
  # Iterate through each category
  for (category in categories){
    # Subset for the current category
    cat_cols <- grep(paste0("^", category), colnames(complete_df))
    cat_df <- complete_df[, cat_cols, drop = FALSE]
    
    # Replace every "NaN", "Inf" and "-Inf" with NA
    cat_df[is.na(cat_df) | cat_df=="Inf" | cat_df=="-Inf"] <- NA
    
    # Remove all rows that have at least one missing value to only get complete matrices
    cat_df <- na.omit(cat_df)
    
    # Determine the indexes by extracting the complete data's row names
    row_indexes <- as.integer(rownames(cat_df))
    
    # Rename the categories inside and out
    sanitized_cols <- gsub("/","_over_", colnames(cat_df))
    colnames(cat_df) <- sanitized_cols
    
    sanitized_category <- gsub("/","_over_", category)
    
    # Store the original data in the current category
    result_list[[sanitized_category]] <- as.data.frame(cat_df, row.names=row_indexes)
  }
  # Reorder the list in alphabetical order
  reordered_list <- result_list[order(names(setNames(result_list, result_list)))]
  
  # Rename those lists back to the original form
  names(reordered_list) <- gsub("_over_","/", names(reordered_list))
  
  return(reordered_list)
}

original_matrices_prot <- get_original_matrix_list(conditions_by_label_proteinGroups)
original_matrices_phos <- get_original_matrix_list(conditions_by_label_Phosphosites)
original_matrices_norm <- get_original_matrix_list(conditions_by_label_normedPhosphosites)

# Get the specific rows from the significant dataframes
filter_rows <- function(original_matrices, significant_dataframe){
  # Create a list to store those significant hits
  result_list <- list()
  
  # Iterate by category
  for (i in (1:length(original_matrices))){
    # Get the row indexes from the original matrix 
    row_names <- rownames(original_matrices[[i]])
    
    # Get the category name of the current element from the matrices list
    category_name <- names(original_matrices)[i]
    
    # Get the column index where the significant dataframe column and the category name matches
    column_index <- grep(category_name, colnames(significant_dataframe))
    
    # Store the original data in the current category
    result_list[[category_name]] <- as.data.frame(significant_dataframe[row_names,column_index], row.names=row_names)
  }
  
  return(result_list)
}

filtered_significant_prot <- filter_rows(original_matrices_prot, condition_significant_proteinGroups)
filtered_significant_phos <- filter_rows(original_matrices_phos, condition_significant_Phosphosites)
filtered_significant_norm <- filter_rows(original_matrices_norm, condition_significant_normedPhosphosites)

# Get all the generated artificial missing values from the different folders (not intensities!) and store them into a list
retrieve_generated_matrices <- function(folder_path){
  # Get a list of all txt files in the folder
  txt_files <- list.files(folder_path, full.names=TRUE)
  # Initialize an empty list to store those matrices
  restored_list <- list()
  
  # Iterate through each txt files and read it into a matrix
  for (single_file in txt_files){
    # Read the txt file into a dataframe
    df_data <- as.data.frame(read.table(single_file, header=TRUE, sep="\t", check.names=FALSE, row.names=1))
    
    # Extract the name of the file
    file_name <- tools::file_path_sans_ext(basename(single_file))
    # Make a grep find to only extract the first part of the file name (before the "." character)
    grep_name <- sub("\\..*","",file_name)
    grep_name <- sub("_over_", "/", grep_name)
    # Extract the index from the file name to identify the matrix position for the list
    matrix_index <- as.numeric(gsub(".*\\D+(\\d+)$", "\\1", file_name))
    
    # Assign the matrix  to the corresponding position in the list
    restored_list[[grep_name]][[matrix_index]] <- df_data
  }
  
  return(restored_list)
}

missing_matrices_prot <- retrieve_generated_matrices("Generated_MV_matrices_proteinGroups")
missing_matrices_phos <- retrieve_generated_matrices("Generated_MV_matrices_Phosphosites")
missing_matrices_norm <- retrieve_generated_matrices("Generated_MV_matrices_normedPhoshposites")

# You will need to replace every column names in every matrix into a readable one without the "/" character for mice package
change_column_names <- function(matrix_list){
  # Iterate through each element of the list that are matrices
  for (i in 1:length(matrix_list)){
    for (j in 1:length(matrix_list[[i]])){
      current_matrix <- matrix_list[[i]][[j]]
      column_names <- colnames(current_matrix)
      # Make a grep find to only extract the "/" character
      grep_names <- sub("/", "_over_", column_names)
      colnames(current_matrix) <- grep_names
      # Assign the change to the matrix list
      matrix_list[[i]][[j]] <- current_matrix
    }
  }
  return(matrix_list)
}

missing_matrices_prot <- change_column_names(missing_matrices_prot)
missing_matrices_phos <- change_column_names(missing_matrices_phos)
missing_matrices_norm <- change_column_names(missing_matrices_norm)

# Make example matrices to test for one matrix imputation per generation type
test_missing_matrix <- missing_matrices_prot[[1]][[1]]

# IMPUTATIONS PART

# Negative control imputation methods
library(mice)

# Random number in min-max range
min_max_range <- function(matrix_column){
  minimum <- min(matrix_column, na.rm = TRUE)
  maximum <- max(matrix_column, na.rm = TRUE)
  return(c(minimum, maximum))
}

# Apply it to the nested lists
random_min_max <- function(matrix_list){
  for (i in 1:length(matrix_list)){
    for (j in 1:length(matrix_list[[i]])){
      current_matrix <- matrix_list[[i]][[j]]
      for (k in 1:ncol(current_matrix)){
        min_max_vector <- min_max_range(current_matrix[,k])
        current_matrix[,k] <- ifelse(is.na(current_matrix[,k]), 
                                     runif(sum(is.na(current_matrix[,k])), min=min_max_vector[1], max=min_max_vector[2]), current_matrix[,k])
      }
      matrix_list[[i]][[j]] <- current_matrix
    }
  }
  return(matrix_list)
}

min_max_imputed_matrices_prot <- random_min_max(missing_matrices_prot)
min_max_imputed_matrices_phos <- random_min_max(missing_matrices_phos)
min_max_imputed_matrices_norm <- random_min_max(missing_matrices_norm)

# Random Gaussian distribution
gaussian_sampling <- function(matrix_list){
  for (i in 1:length(matrix_list)){
    for (j in 1:length(matrix_list[[i]])){
      current_matrix <- matrix_list[[i]][[j]]
      # Store the row indexes of the matrix for easy identification on later processes
      row_indexes <- as.integer(rownames(current_matrix))
      
      # Make the imputation process
      imputed_matrix <- mice(current_matrix, method="norm", m=1, maxit=1, print=FALSE)
      
      # Get the completed matrix as it was imputed with proper row indexes
      complete_matrix <- complete(imputed_matrix)
      matrix_list[[i]][[j]] <- as.data.frame(complete_matrix, row.names=row_indexes)
    }
  }
  return(matrix_list)
}

gaussian_imputed_matrices_prot <- gaussian_sampling(missing_matrices_prot)
gaussian_imputed_matrices_phos <- gaussian_sampling(missing_matrices_phos)
gaussian_imputed_matrices_norm <- gaussian_sampling(missing_matrices_norm)

# Actual imputation methods (libraries required)

# kNN algorithm (for reproducibility, the random seed is by default 362436069)
library(impute)

# There are two versions of this function imputation depending on the data_type
knn_imputation_prot <- function(matrix_list){
  for (i in 1:length(matrix_list)){
    for (j in 1:length(matrix_list[[i]])){
      current_matrix <- matrix_list[[i]][[j]]
      
      # Store the row indexes of the matrix for easy identification on later processes
      row_indexes <- as.integer(rownames(current_matrix))
      
      # Convert non-numeric columns to numeric values
      data_numeric <- sapply(current_matrix, as.numeric)
      
      # Calculate the ideal k as the square root of the total number of proteins samples
      ideal_k <- round(sqrt(nrow(data_numeric)))
      
      # Impute the knn method with default max block size
      imputed_matrix <- impute.knn(data_numeric, k=ideal_k, rowmax=1)
      complete_matrix <- imputed_matrix$data
      # If needed, convert back to dataframe
      matrix_list[[i]][[j]] <- as.data.frame(complete_matrix, row.names=row_indexes)
    }
  }
  return(matrix_list)
}

knn_imputation_phos <- function(matrix_list){
  for (i in 1:length(matrix_list)){
    for (j in 1:length(matrix_list[[i]])){
      current_matrix <- matrix_list[[i]][[j]]
      
      # Store the row indexes of the matrix for easy identification on later processes
      row_indexes <- as.integer(rownames(current_matrix))
      
      # Convert non-numeric columns to numeric values
      data_numeric <- sapply(current_matrix, as.numeric)
      
      # Calculate the ideal k as the square root of the total number of phoshposites samples
      ideal_k <- round(sqrt(nrow(data_numeric)))
      
      # Impute the knn method with an increased block
      imputed_matrix <- impute.knn(data_numeric, k=ideal_k, maxp=5000, rowmax=1)
      complete_matrix <- imputed_matrix$data
      # If needed, convert back to dataframe
      matrix_list[[i]][[j]] <- as.data.frame(complete_matrix, row.names=row_indexes)
    }
  }
  return(matrix_list)
}

knn_imputed_matrices_prot <- knn_imputation_prot(missing_matrices_prot)
knn_imputed_matrices_phos <- knn_imputation_phos(missing_matrices_phos)
knn_imputed_matrices_norm <- knn_imputation_phos(missing_matrices_norm)

# Maximum Likelihood Estimation (MLE)
library(imp4p)

mle_imputation <- function(matrix_list){
  for (i in 1:length(matrix_list)){
    for (j in 1:length(matrix_list[[i]])){
      current_matrix <- matrix_list[[i]][[j]]
      
      # Get the column factors for replicates considered as 1 condition and impute with mle method using EM algorithm
      column_factors <- factor(rep(1, ncol(current_matrix)))
      imputed_matrix <- impute.mle(tab=as.matrix(current_matrix), conditions=column_factors)
      
      # Get the completed matrix as it was imputed with proper row indexes
      matrix_list[[i]][[j]] <- as.data.frame(imputed_matrix)
    }
  }
  return(matrix_list)
}

mle_imputed_matrices_prot <- mle_imputation(missing_matrices_prot)
mle_imputed_matrices_phos <- mle_imputation(missing_matrices_phos)
mle_imputed_matrices_norm <- mle_imputation(missing_matrices_norm)

# Random forest
library(missForest)

rf_imputation <- function(matrix_list){
  for (i in 1:length(matrix_list)){
    for (j in 1:length(matrix_list[[i]])){
      # Time the process per matrix
      start_time <- Sys.time()
      current_matrix <- matrix_list[[i]][[j]]
      
      # Impute the random forest method imputation to the current matrix
      imputed_matrix <- missForest(current_matrix)
      matrix_list[[i]][[j]] <- imputed_matrix$ximp
      
      # Stop the time process here
      end_time <- Sys.time()
      time_taken <- end_time - start_time
      
      # Print the current state and time of which imputed matrix is done for the console information
      print(paste0("RF imputation for (",i,",",j,") done in ", time_taken, "minutes"))
    }
  }
  return(matrix_list)
}

# Only run this in a remote computer for a while as the computation time is very lengthy (estimated time of 1 day for a whole dataset folder)
rf_imputed_matrices_prot <- rf_imputation(missing_matrices_prot)
rf_imputed_matrices_phos <- rf_imputation(missing_matrices_phos)
rf_imputed_matrices_norm <- rf_imputation(missing_matrices_norm)

# Identify DE proteins/phosphosites from those imputed matrices
# Essentially a copy from BenchmarkStatistics.R, from the mlm implementation to the addition of significant ids

# Function to find the significant proteins/phosphosites from a CI table
find_significant_lme <- function(CI_table){
  if(!is.nan(CI_table[1]) & !is.nan(CI_table[2])){
    if((CI_table[1]<0) & (CI_table[2]>0)){
      significant<-c("")
    }else{
      significant<-c("+")
    }
  }else{
    significant<-c("")
  }
  
  return(significant)
}

# Let's apply the mixed linear model with lme4 package (might yield to boundary problems)
library(lme4)

# Application part of the mixed linear model

# Function that applies the mixed linear model (if the test is complete and with no problems...)
mixed_linear_model_significant <- function(complete_df, norm_test = FALSE) {
  # Transform the data in log2 base to apply log2 fold change except if it already is (normed phoshposites cases)
  if (norm_test) {
    log2_complete_df <- complete_df
  }else{
    log2_complete_df <- log2(complete_df)
  }
  
  # Get the unique categories from the complete dataset
  categories <- unique(gsub("\\d+$","", colnames(complete_df)))
  
  # Create a list to store the results for each category
  result_list <- list()
  
  # Iterate through each category
  for (category in categories){
    # Subset the log2 dataframe for the current category
    cat_cols <- grep(paste0("^", category), colnames(log2_complete_df))
    cat_df <- log2_complete_df[, cat_cols, drop = FALSE]
    
    # Replace every "NaN", "Inf" and "-Inf" with NA
    cat_df[is.na(cat_df) | cat_df=="Inf" | cat_df=="-Inf"] <- NA
    
    # Remove all rows that have at least one missing value to only get complete matrices
    cat_df <- na.omit(cat_df)
    
    # Determine the indexes by extracting the complete data's row names
    row_indexes <- as.integer(rownames(cat_df))
    
    # Determine the number of protein groups/phosphosites from the rows of the complete matrix
    nb_complete <- nrow(cat_df)
    
    # Determine the number of replicates from the number of columns from cat_df
    nb_rep <- ncol(cat_df)
    
    # Linear modelling stuff
    positions <- rep(1:nb_complete, times=nb_rep) # construct random effect factor column: protein group/phosphosites
    replicates <- rep(letters[1:nb_rep], each=nb_complete) # construct random effect factor column: replicate
    
    # You would need to apply the log2fold change yourself (especially for intensities)
    # Construct target value column: log2fold ratios 
    log2fold_vector <- c(cat_df[,1])
    for(i in c(2:nb_rep)){
      log2fold_vector <- c(log2fold_vector,cat_df[,i])
    }
    
    # Create a dataframe that mimic the expected values for linear modelling with all the constructed response values and random effects
    linear_df <- data.frame(log2fold_vector, positions, replicates)
    
    # Create the mixed linear model and extract its random effects
    linear_model <- lmer(log2fold_vector ~ 1 + (1|replicates) + (1|positions), na.action=na.omit, data=linear_df)
    fit <- ranef(linear_model, condVar=TRUE)
    
    # Create empty matrix that summarizes the fit model with column data
    column_names <- c("fit_log2fold", "var_log2fold", "sd_log2fold", "nb_replicate", "CI_low", "CI_high", "significant")
    summary_df <- as.data.frame(matrix(NA, nrow=nb_complete, ncol=length(column_names)), row.names=row_indexes)
    colnames(summary_df) <- column_names
    
    # Fill the columns with linear model fit data
    summary_df$fit_log2fold <- fit$positions$`(Intercept)`
    
    summary_df$var_log2fold <- as.vector(attr(fit$positions,'postVar'))
    summary_df$sd_log2fold <- sqrt(summary_df$var_log2fold)
    
    summary_df$nb_replicate <- nb_rep
    
    summary_df$CI_low <- summary_df$fit_log2fold - qt(0.975, df = summary_df$nb_replicate) * summary_df$sd_log2fold
    summary_df$CI_high <- summary_df$fit_log2fold + qt(0.975, df = summary_df$nb_replicate) * summary_df$sd_log2fold
    
    # Find all the significant proteins/phosphosites from the CI table
    CI <- base::cbind(summary_df$CI_low, summary_df$CI_high)
    summary_df$significant <- apply(CI, 1, find_significant_lme)
    
    # Identify the number of significant proteins/phosphosites and store them as row ids
    row_ids_significant_lmer <- as.integer(row.names(summary_df[which(summary_df$significant=="+"),]))
    
    # Store the results for the current category
    result_list[[category]] <- row_ids_significant_lmer
  }
  
  return (result_list)
}

# Significance A test on distribution of mean values (or outlier ratio significance)
# Do it two-sided rather than one sided for significance test
library(pracma)

# Function to find the significant proteins/phosphosites from a p-value table
find_significant_sigA <- function(sigA_vector){
  significant <- c()
  for (s in sigA_vector){
    if(s <= 0.05){
      significant <- c(significant, "+")
    }else{
      significant <- c(significant, "")
    }
  }
  return(significant)
}

# Function that applies the outlier or significance A test per complete dataframe (without MVs)
outlier_test_significant <- function(complete_df, have_norm = FALSE){
  # Transform the data in log2 base to apply log2 fold change except if it already is (normed phoshposites cases)
  if (have_norm){
    log2_complete_df <- complete_df
  }else{
    log2_complete_df <- log2(complete_df)
  }
  
  # Get the unique categories from the complete dataset
  categories <- unique(gsub("\\d+$","", colnames(complete_df)))
  
  # Create a list to store the results for each category
  result_list <- list()
  
  # Iterate through the different categories
  for (category in categories){
    # Subset dataframe for the current category
    cat_cols <- grep(paste0("^", category), colnames(log2_complete_df))
    cat_df <- log2_complete_df[, cat_cols, drop = FALSE]
    
    # Replace every "NaN", "Inf" and "-Inf" with NA
    cat_df[is.na(cat_df) | cat_df=="Inf" | cat_df=="-Inf"] <- NA
    
    # Remove all rows that have at least one missing value to only get complete matrices
    cat_df <- na.omit(cat_df)
    
    # Determine the indexes by extracting the complete data's row names
    row_indexes <- as.integer(rownames(cat_df))
    
    # Determine the number of protein groups/phosphosites from the rows of the complete matrix
    nb_complete <- nrow(cat_df)
    
    # Create an empty dataframe to store all the calculated values
    column_names <- c("mean", "z_distance", "p_value", "p_adjusted", "significant")
    summary_df <- as.data.frame(matrix(NA, nrow=nb_complete, ncol=length(column_names)), row.names=row_indexes)
    colnames(summary_df) <- column_names
    
    # Compute the mean for each row that will be the crux of the distribution
    row_means <- as.vector(rowMeans(cat_df))
    summary_df$mean <- row_means
    
    # Create the distribution of mean values through a histogram for visualization
    hist(row_means, breaks="Scott", main=paste("Distribution of means of", category), xlab="Mean value")
    
    # Get the percentiles from the main distribution as explained in the paper
    r_neg <- as.numeric(quantile(row_means, probs=0.1587))
    r_0 <- as.numeric(quantile(row_means, probs=0.5))
    r_pos <- as.numeric(quantile(row_means, probs=0.8413))
    
    # Introduce a vector z that stores every distance measures per row id
    z <- c()
    
    # Calculate the distance measure according to the percentiles per ratio value
    for(r in row_means){
      if (r > r_0){
        z <- c(z, (r-r_0) / (r_pos-r_0))
      }else{
        z <- c(z, (r_0-r) / (r_0-r_neg))
      }
    }
    summary_df$z_distance <- z
    
    # Calculate the p-value for detection of significant outlier ratios
    significance_A <- erfc(z/sqrt(2))
    summary_df$p_value <- significance_A
    
    # Apply the Bonferroni-Hochberg correction for multiple testing to the set of p-values
    summary_df$p_adjusted <- p.adjust(summary_df$p_value, method="BH")
    
    # Only consider those who have an adjusted p-value of less than 0.05 as significant outliers after BH correction
    summary_df$significant <- find_significant_sigA(summary_df$p_adjusted)
    
    # Identify the number of significant proteins/phosphosites and store them as row ids
    row_ids_significant_sigA <- as.integer(row.names(summary_df[which(summary_df$significant=="+"),]))
    
    # Store the results for the current category
    result_list[[category]] <- row_ids_significant_sigA
  }
  
  return(result_list)
}

# Combine the two tests significant results to get DE proteins/phosphosites ids as vector or as matrix
# The idea is to get the two lists of mlm and sigA and get significant ids per conditions
compare_lists_ids <- function(combined_list){
  # Find similarly named experiment variables that stores ids
  common_names <- intersect(names(combined_list$mlm), names(combined_list$sigA))
  
  # Initialize a list to store common variables
  common_list <- list()
  
  # Iterate over similarly named id variables
  for (name in common_names){
    # Convert id variables to sets
    set1 <- as.list(combined_list$mlm[[name]])
    set2 <- as.list(combined_list$sigA[[name]])
    
    # Find common elements
    common_ids <- intersect(set1, set2)
    
    # Store common id variable in common_list
    common_list[[name]] <- unlist(common_ids)
  }
  
  return(common_list)
}

# Add columns to the original data to signify if the proteins/peptides are significantly DE depending on the conditions with a "+"
add_id_columns <- function(original_df, id_list){
  # Create a matrix to store the values with original row names from the original dataframe
  significant_matrix <- matrix(NA, nrow=nrow(original_df), ncol=length(id_list))
  row.names(significant_matrix) <- row.names(original_df)
  
  # Iterate through the list of IDs
  for (i in seq_along(id_list)){
    id <- names(id_list)[i]
    significant_matrix[,i] <- ifelse(rownames(original_df) %in% id_list[[id]], "+", "")
  }
  # Convert the matrix to a dataframe
  result_df<- as.data.frame(significant_matrix)
  colnames(result_df) <- sub("_over_", "/", names(id_list))
  
  return(result_df)
}

# Function that calls the mlm, sigA and list comparisons with matrix list and original matrices as input arguments
call_significant <- function(matrix_list, original_list){
  # Test if the matrix list contains normed phosphosites or not as a logical value
  have_norm <- grepl("norm", deparse(substitute(matrix_list)))
  
  # Make some space for displaying histogram figures
  list_length <- length(matrix_list[[1]])
  col_mfrow <- ceiling(sqrt(list_length))
  row_mfrow <- ceiling(list_length / col_mfrow)
  par(mfrow=c(row_mfrow,col_mfrow), mar=rep(2,4))
  
  # Iterate over the matrix list
  for (i in 1:length(matrix_list)){
    # Get the original matrix from the category of choice
    original_matrix <- original_list[[i]]
    for (j in 1:length(matrix_list[[i]])){
      current_matrix <- matrix_list[[i]][[j]]
      
      # Call the mlm model to the current matrix
      significant_mlm_matrix <- mixed_linear_model_significant(current_matrix, have_norm)
      
      # Call the sigA model to the current_matrix
      significant_sigA_matrix <- outlier_test_significant(current_matrix, have_norm)
      
      # Store those two results in a combined list placement for the function call to compare them
      combined_results_list <- list(mlm=significant_mlm_matrix, sigA=significant_sigA_matrix)
      compared_list <- compare_lists_ids(combined_results_list)
      
      # Add "+" to those significant ids
      matrix_list[[i]][[j]] <- add_id_columns(original_matrix, compared_list)
      
      # Print the current state of which matrix is done for the information
      print(paste0("(",i,",",j,") done"))
    }
  }
  return(matrix_list)
}

# We can call once without any imputation methods for comparison as blank control (especially for false negatives)
significant_without_prot <- call_significant(missing_matrices_prot, original_matrices_prot)
significant_without_phos <- call_significant(missing_matrices_phos, original_matrices_phos)
significant_without_norm <- call_significant(missing_matrices_norm, original_matrices_norm)

# Eventual calls to identify DE proteins/phosphosites for min_max imputation 
# (remark: so few or no significant phosphosites due to significance A test finding nothing, expected for pure randomness?)
significant_min_max_prot <- call_significant(min_max_imputed_matrices_prot, original_matrices_prot)
significant_min_max_phos <- call_significant(min_max_imputed_matrices_phos, original_matrices_phos)
significant_min_max_norm <- call_significant(min_max_imputed_matrices_norm, original_matrices_norm)

# Eventual calls to identify DE proteins/phosphosites for gaussian imputation (better than min_max for negative control)
significant_gaussian_prot <- call_significant(gaussian_imputed_matrices_prot, original_matrices_prot)
significant_gaussian_phos <- call_significant(gaussian_imputed_matrices_phos, original_matrices_phos)
significant_gaussian_norm <- call_significant(gaussian_imputed_matrices_norm, original_matrices_norm)

# Eventual calls to identify DE proteins/phosphosites for knn imputation
significant_knn_prot <- call_significant(knn_imputed_matrices_prot, original_matrices_prot)
significant_knn_phos <- call_significant(knn_imputed_matrices_phos, original_matrices_phos)
significant_knn_norm <- call_significant(knn_imputed_matrices_norm, original_matrices_norm)

# Eventual calls to identify DE proteins/phosphosites for mle imputation
significant_mle_prot <- call_significant(mle_imputed_matrices_prot, original_matrices_prot)
significant_mle_phos <- call_significant(mle_imputed_matrices_phos, original_matrices_phos)
significant_mle_norm <- call_significant(mle_imputed_matrices_norm, original_matrices_norm)

# Eventual calls to identify DE proteins/phosphosites for rf imputation
significant_rf_prot <- call_significant(rf_imputed_matrices_prot, original_matrices_prot)
significant_rf_phos <- call_significant(rf_imputed_matrices_phos, original_matrices_phos)
significant_rf_norm <- call_significant(rf_imputed_matrices_norm, original_matrices_norm)

# Benchmark the different strategies with the original significant matrix (confusion matrices per method)

# Make confusion matrices with the filtered significant list and the imputed significant list
library(caret)
make_confusion_matrices <- function(original_significant_list, imputed_significant_list){
  # Create a list to store all the confusion matrices per category
  result_list <- list()
  
  for (i in 1:length(imputed_significant_list)){
    # Get the column name of the original significant matrix at position i
    original_name <- names(original_significant_list)[i]
    # Get the data vector of the original significant matrix
    original_data <- original_significant_list[[i]]
    # Replace every "" characters into 0 and every "+" characters into 1
    original_data[original_data==""] <- 0
    original_data[original_data=="+"] <- 1
    # Make the original dataframe a factor data (using lapply which is counter-intuitive)
    factored_original_data <- unlist(lapply(original_data, factor))
    
    for (j in 1:length(imputed_significant_list[[i]])){
      # Get the data vector of the imputed significant list
      imputed_data <- imputed_significant_list[[i]][[j]]
      
      # Make the confusion matrix calculation if the imputed data has no missing column
      if (ncol(imputed_data)!=0){
        # Replace every "" characters into 0 and every "+" characters into 1
        imputed_data[imputed_data==""] <- 0
        imputed_data[imputed_data=="+"] <- 1
        # Make the imputed dataframe a factor data
        factored_imputed_data <- unlist(lapply(imputed_data, factor))
        
        confusion_matrix <- confusionMatrix(data=factored_imputed_data, reference=factored_original_data)
        result_list[[original_name]][[j]] <- confusion_matrix$table
      }else{
        # If the column is missing, assign it to "NA" rather than "NULL"
        result_list[[original_name]][[j]] <- NA
      }
    }
  }
  
  return(result_list)
}

# Apply the confusion matrices calculation to every applicable data (even with blank control)
confusion_matrices_without_prot <- make_confusion_matrices(filtered_significant_prot, significant_without_prot)
confusion_matrices_without_phos <- make_confusion_matrices(filtered_significant_phos, significant_without_phos)
confusion_matrices_without_norm <- make_confusion_matrices(filtered_significant_norm, significant_without_norm)

confusion_matrices_min_max_prot <- make_confusion_matrices(filtered_significant_prot, significant_min_max_prot)
confusion_matrices_min_max_phos <- make_confusion_matrices(filtered_significant_phos, significant_min_max_phos)
confusion_matrices_min_max_norm <- make_confusion_matrices(filtered_significant_norm, significant_min_max_norm)

confusion_matrices_gaussian_prot <- make_confusion_matrices(filtered_significant_prot, significant_gaussian_prot)
confusion_matrices_gaussian_phos <- make_confusion_matrices(filtered_significant_phos, significant_gaussian_phos)
confusion_matrices_gaussian_norm <- make_confusion_matrices(filtered_significant_norm, significant_gaussian_norm)

confusion_matrices_knn_prot <- make_confusion_matrices(filtered_significant_prot, significant_knn_prot)
confusion_matrices_knn_phos <- make_confusion_matrices(filtered_significant_phos, significant_knn_phos)
confusion_matrices_knn_norm <- make_confusion_matrices(filtered_significant_norm, significant_knn_norm)

confusion_matrices_mle_prot <- make_confusion_matrices(filtered_significant_prot, significant_mle_prot)
confusion_matrices_mle_phos <- make_confusion_matrices(filtered_significant_phos, significant_mle_phos)
confusion_matrices_mle_norm <- make_confusion_matrices(filtered_significant_norm, significant_mle_norm)

confusion_matrices_rf_prot <- make_confusion_matrices(filtered_significant_prot, significant_rf_prot)
confusion_matrices_rf_phos <- make_confusion_matrices(filtered_significant_phos, significant_rf_phos)
confusion_matrices_rf_norm <- make_confusion_matrices(filtered_significant_norm, significant_rf_norm)

# Average all confusion matrices (CM) from a list
average_CM <- function(confusion_matrices_list){
  # Create a list to store all the average CMs per category
  result_list <- list()
  
  for (i in 1:length(confusion_matrices_list)){
    # Get the element name of the list at position i
    original_name <- names(confusion_matrices_list)[i]
    
    # If all of the elements on the list are "NA", label the average matrix as "NA" and ignore the next commands
    if(all(is.na(confusion_matrices_list[[i]]))){
      result_list[[original_name]] <- NA; next
    }
    
    # Ignore the confusion matrices that are not there in the list (labelled as "NA")
    refined_list <- confusion_matrices_list[[i]][!sapply(confusion_matrices_list[[i]], anyNA)]
    
    # Calculate the total number of matrices from the category
    num_matrices <- length(refined_list)
    # Initialize a matrix that contains the same class typo as the original confusion matrix
    sum_matrix <- refined_list[[1]]
    # Sum up all the matrices from this category except the first one
    for (m in refined_list[-1]){
      sum_matrix <- sum_matrix + m
    }
    # Average out the sum
    average_matrix <- sum_matrix / num_matrices
    
    # Store this average CM in the final list
    result_list[[original_name]] <- average_matrix
  }
  
  return(result_list)
}

average_CM_without_prot <- average_CM(confusion_matrices_without_prot)
average_CM_without_phos <- average_CM(confusion_matrices_without_phos)
average_CM_without_norm <- average_CM(confusion_matrices_without_norm)

average_CM_min_max_prot <- average_CM(confusion_matrices_min_max_prot)
average_CM_min_max_phos <- average_CM(confusion_matrices_min_max_phos)
average_CM_min_max_norm <- average_CM(confusion_matrices_min_max_norm)

average_CM_gaussian_prot <- average_CM(confusion_matrices_gaussian_prot)
average_CM_gaussian_phos <- average_CM(confusion_matrices_gaussian_phos)
average_CM_gaussian_norm <- average_CM(confusion_matrices_gaussian_norm)

average_CM_knn_prot <- average_CM(confusion_matrices_knn_prot)
average_CM_knn_phos <- average_CM(confusion_matrices_knn_phos)
average_CM_knn_norm <- average_CM(confusion_matrices_knn_norm)

average_CM_mle_prot <- average_CM(confusion_matrices_mle_prot)
average_CM_mle_phos <- average_CM(confusion_matrices_mle_phos)
average_CM_mle_norm <- average_CM(confusion_matrices_mle_norm)

average_CM_rf_prot <- average_CM(confusion_matrices_rf_prot)
average_CM_rf_phos <- average_CM(confusion_matrices_rf_phos)
average_CM_rf_norm <- average_CM(confusion_matrices_rf_norm)

# Make sure that those averages are grouped up by datatype as variables
average_CM_overall_proteinGroups <- list(no=average_CM_without_prot, minmax=average_CM_min_max_prot, gaussian=average_CM_gaussian_prot, knn=average_CM_knn_prot, mle=average_CM_mle_prot, rf=average_CM_rf_prot)
average_CM_overall_Phophosites <- list(no=average_CM_without_phos, minmax=average_CM_min_max_phos, gaussian=average_CM_gaussian_phos, knn=average_CM_knn_phos, mle=average_CM_mle_phos, rf=average_CM_rf_phos)
average_CM_overall_normedPhosphosites <- list(no=average_CM_without_norm, minmax=average_CM_min_max_norm, gaussian=average_CM_gaussian_norm, knn=average_CM_knn_norm, mle=average_CM_mle_norm, rf=average_CM_rf_norm)

# Save those average confusion matrices in RData format for visualization later
save(average_CM_overall_proteinGroups, file="average_CM_proteinGroups.RData")
save(average_CM_overall_Phophosites, file="average_CM_Phospho(STY)Sites.RData")
save(average_CM_overall_normedPhosphosites, file="average_CM_Phospho(STY)Sites_normed.RData")

# Metrics calculations for each imputation type with average and standard deviation in mind
calculate_metrics <- function(average_cm_list, cm_list, imputation_type){
  # Create a list to store all the metrics of CMs per category
  result_list <- list()
  for (i in 1:length(average_cm_list)){
    # Get the element name of the list at position i
    original_name <- names(average_cm_list)[i]
    
    # Get the average CM from the list
    conf_matrix <- average_cm_list[[i]]
    
    # If the current average CM is NA, the metrics will be NA too and we ignore the next commands
    if(all(is.na(conf_matrix))){
      metrics <- data.frame(FDR = NA, FDR_sd = NA, FNR = NA, FNR_sd = NA , Accuracy = NA, Accuracy_sd = NA, MCC = NA, MCC_sd = NA)
      result_list[[paste0(imputation_type,"_",original_name)]] <- metrics
      next
    }
    
    # Get the values themselves from the confusion matrix individually
    TP <- conf_matrix[2,2] # True Positives
    TN <- conf_matrix[1,1] # True Negatives
    FP <- conf_matrix[2,1] # False Positives
    FN <- conf_matrix[1,2] # False Negatives
    
    # Calculate precision (or False Discovery Rate)
    FDR_value <- FP / (FP + TP)
    
    # Calculate recall (or False Negative Rate)
    FNR_value <- FN / (FN + TP)
    
    # Calculate accuracy (TP + TN / total)
    accuracy_value <- (TP + TN) / (TP + TN + FP + FN)
    
    # Calculate Matthews correlation coefficient (see formula below)
    numerator  <- TP*TN - FP*FN
    denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    MCC_value <- numerator / denominator
    
    # To get the variations or standard deviations, we'll need to look at the confusion matrices list
    # Initialize two empty dataframes where all the confusion values are stored by columns (in this order: TN, FP, FN, TP)
    # and where all the metrics are stored by columns (in this order: FDR, FNR, Accuracy, MCC)
    CM_table <- data.frame(matrix(ncol=4, nrow=20, dimnames=list(NULL, c("TN","FP","FN","TP"))))
    metrics_table <- data.frame(matrix(ncol=4, nrow=20, dimnames=list(NULL, c("FDR","FNR","Accuracy","MCC"))))
    
    for (j in 1:length(cm_list[[i]])){
      CM <- cm_list[[i]][[j]]
      
      if(all(is.na(CM))){next}
      
      TP_value <- as.numeric(CM[2,2])
      TN_value <- as.numeric(CM[1,1])
      FP_value <- as.numeric(CM[2,1])
      FN_value <- as.numeric(CM[1,2])
      
      # Store those confusion values in the current row (it is the random matrix number)
      CM_table[j,] <- c(TN_value, FP_value, FN_value, TP_value)
      
      # Calculate precision (or False Discovery Rate)
      FDR_val <- FP_value / (FP_value + TP_value)
      
      # Calculate recall (or False Negative Rate)
      FNR_val <- FN_value / (FN_value + TP_value)
      
      # Calculate accuracy (TP + TN / total)
      accuracy_val <- (TP_value + TN_value) / (TP_value + TN_value + FP_value + FN_value)
      
      # Calculate Matthews correlation coefficient (see formula below)
      num  <- TP_value*TN_value - FP_value*FN_value
      denom <- sqrt((TP_value + FP_value) * (TP_value + FN_value) * (TN_value + FP_value) * (TN_value + FN_value))
      MCC_val <- num / denom
      
      # Store those metrics values in the current row
      metrics_table[j,] <- c(FDR_val, FNR_val, accuracy_val, MCC_val)
    }
    
    # Combine all metrics average and sd into a dataframe
    metrics <- data.frame(FDR = FDR_value, FDR_sd = sd(metrics_table$FDR, na.rm=TRUE), 
                          FNR = FNR_value, FNR_sd = sd(metrics_table$FNR, na.rm=TRUE),
                          Accuracy = accuracy_value, Accuracy_sd = sd(metrics_table$Accuracy, na.rm=TRUE), 
                          MCC = MCC_value, MCC_sd = sd(metrics_table$MCC, na.rm=TRUE))
    
    # Store those metrics by category type
    result_list[[paste0(imputation_type,"_",original_name)]] <- metrics
  }
  # Convert the list into a dataframe and return it
  list_to_df <- do.call(rbind.data.frame, result_list)
  
  return(list_to_df)
}

without_metrics_prot <- calculate_metrics(average_CM_without_prot, confusion_matrices_without_prot, "without")
without_metrics_phos <- calculate_metrics(average_CM_without_phos, confusion_matrices_without_phos, "without")
without_metrics_norm <- calculate_metrics(average_CM_without_norm, confusion_matrices_without_norm, "without")

min_max_metrics_prot <- calculate_metrics(average_CM_min_max_prot, confusion_matrices_min_max_prot, "minmax")
min_max_metrics_phos <- calculate_metrics(average_CM_min_max_phos, confusion_matrices_min_max_phos, "minmax")
min_max_metrics_norm <- calculate_metrics(average_CM_min_max_norm, confusion_matrices_min_max_norm, "minmax")

gaussian_metrics_prot <- calculate_metrics(average_CM_gaussian_prot, confusion_matrices_gaussian_prot, "gaussian")
gaussian_metrics_phos <- calculate_metrics(average_CM_gaussian_phos, confusion_matrices_gaussian_phos, "gaussian")
gaussian_metrics_norm <- calculate_metrics(average_CM_gaussian_norm, confusion_matrices_gaussian_norm, "gaussian")

knn_metrics_prot <- calculate_metrics(average_CM_knn_prot, confusion_matrices_knn_prot, "knn")
knn_metrics_phos <- calculate_metrics(average_CM_knn_phos, confusion_matrices_knn_phos, "knn")
knn_metrics_norm <- calculate_metrics(average_CM_knn_norm, confusion_matrices_knn_norm, "knn")

mle_metrics_prot <- calculate_metrics(average_CM_mle_prot, confusion_matrices_mle_prot, "mle")
mle_metrics_phos <- calculate_metrics(average_CM_mle_phos, confusion_matrices_mle_phos, "mle")
mle_metrics_norm <- calculate_metrics(average_CM_mle_norm, confusion_matrices_mle_norm, "mle")

rf_metrics_prot <- calculate_metrics(average_CM_rf_prot, confusion_matrices_rf_prot, "rf")
rf_metrics_phos <- calculate_metrics(average_CM_rf_phos, confusion_matrices_rf_phos, "rf")
rf_metrics_norm <- calculate_metrics(average_CM_rf_norm, confusion_matrices_rf_norm, "rf")

# Store those results in row-binded tables for each datatype
metrics_CM_proteinGroups <- base::rbind(without_metrics_prot, min_max_metrics_prot, gaussian_metrics_prot, knn_metrics_prot, mle_metrics_prot, rf_metrics_prot)
metrics_CM_Phosphosites <- base::rbind(without_metrics_phos, min_max_metrics_phos, gaussian_metrics_phos, knn_metrics_phos, mle_metrics_phos, rf_metrics_phos)
metrics_CM_normedPhosphosites <- base::rbind(without_metrics_norm, min_max_metrics_norm, gaussian_metrics_norm, knn_metrics_norm, mle_metrics_norm, rf_metrics_norm)

# Other metrics

# Normalized Root Mean Squared Error
# For example you could use the nrmse function from missForest package if you have imputed, generated MV and original data

calculate_nrmse <- function(imputed_list, missing_list, original_list){
  # Create a list to store all the nrmse values per matrix
  nrmse_list <- list()
  
  # Iterate over the matrix lists
  for (i in 1:length(imputed_list)){
    # Get the element name of the list at position i
    original_name <- names(original_list[i])
    
    # Get the original matrix that need to be compared with the other ones
    original_matrix <- original_list[[i]]
    
    for (j in 1:length(imputed_list[[i]])){
      # Get all the necessary imputed and missing matrices that are in each list
      imputed_matrix <- imputed_list[[i]][[j]]
      missing_matrix <- missing_list[[i]][[j]]
      
      # Compute its NRMSE
      nrmse_value <- nrmse(imputed_matrix, missing_matrix, original_matrix)
      
      # Store this value in the nrmse list
      nrmse_list[[original_name]][[j]] <- nrmse_value
    }
  }
  
  return(nrmse_list)
}

nrmse_minmax_prot <- calculate_nrmse(min_max_imputed_matrices_prot, missing_matrices_prot, original_matrices_prot)
nrmse_minmax_phos <- calculate_nrmse(min_max_imputed_matrices_phos, missing_matrices_phos, original_matrices_phos)
nrmse_minmax_norm <- calculate_nrmse(min_max_imputed_matrices_norm, missing_matrices_norm, original_matrices_norm)

nrmse_gaussian_prot <- calculate_nrmse(gaussian_imputed_matrices_prot, missing_matrices_prot, original_matrices_prot)
nrmse_gaussian_phos <- calculate_nrmse(gaussian_imputed_matrices_phos, missing_matrices_phos, original_matrices_phos)
nrmse_gaussian_norm <- calculate_nrmse(gaussian_imputed_matrices_norm, missing_matrices_norm, original_matrices_norm)

nrmse_knn_prot <- calculate_nrmse(knn_imputed_matrices_prot, missing_matrices_prot, original_matrices_prot)
nrmse_knn_phos <- calculate_nrmse(knn_imputed_matrices_phos, missing_matrices_phos, original_matrices_phos)
nrmse_knn_norm <- calculate_nrmse(knn_imputed_matrices_norm, missing_matrices_norm, original_matrices_norm)

nrmse_mle_prot <- calculate_nrmse(mle_imputed_matrices_prot, missing_matrices_prot, original_matrices_prot)
nrmse_mle_phos <- calculate_nrmse(mle_imputed_matrices_phos, missing_matrices_phos, original_matrices_phos)
nrmse_mle_norm <- calculate_nrmse(mle_imputed_matrices_norm, missing_matrices_norm, original_matrices_norm)

nrmse_rf_prot <- calculate_nrmse(rf_imputed_matrices_prot, missing_matrices_prot, original_matrices_prot)
nrmse_rf_phos <- calculate_nrmse(rf_imputed_matrices_phos, missing_matrices_phos, original_matrices_phos)
nrmse_rf_norm <- calculate_nrmse(rf_imputed_matrices_norm, missing_matrices_norm, original_matrices_norm)

#Average out those different values and make it a proper dataframe
average_nrmse <- function(nrmse_list, imputation_type){
  # Create a list to store all the average NRMSE values per category
  result_list <- list()
  
  for (i in 1:length(nrmse_list)){
    # Get the element name of the list at position i
    original_name <- names(nrmse_list)[i]
    
    # Calculate the total number of values from the category
    num_values <- length(nrmse_list[[i]])
    
    # Initialize the sum with the first value given
    sum_values <- nrmse_list[[i]][[1]]
    
    # Sum up all the values from this category except the first one
    for (v in nrmse_list[[i]][-1]){
      sum_values <- sum_values + v
    }
    
    # Average out the sum
    average_value <- sum_values / num_values
    
    # To get the variations or standard deviations, we'll need to look at the nrmse list
    # Initialize an empty table with all the values from the NRMSE calculations
    NRMSE_table <- data.frame(matrix(ncol=1, nrow=20, dimnames=list(NULL, c("NRMSE"))))
    
    for (j in 1:length(nrmse_list[[i]])){
      NRMSE_table[j,] <- nrmse_list[[i]][[j]]
    }
    
    # Store this average NRMSE and sd in the final list
    result_list[[paste0(imputation_type, "_", original_name)]] <- data.frame(NRMSE=average_value, NRMSE_sd=sd(NRMSE_table$NRMSE, na.rm=TRUE))
  }
  # Convert the list into a dataframe
  list_to_df <- do.call(rbind.data.frame, result_list)
  
  return(list_to_df)
}

avg_nrmse_minmax_prot <- average_nrmse(nrmse_minmax_prot, "minmax")
avg_nrmse_minmax_phos <- average_nrmse(nrmse_minmax_phos, "minmax")
avg_nrmse_minmax_norm <- average_nrmse(nrmse_minmax_norm, "minmax")

avg_nrmse_gaussian_prot <- average_nrmse(nrmse_gaussian_prot, "gaussian")
avg_nrmse_gaussian_phos <- average_nrmse(nrmse_gaussian_phos, "gaussian")
avg_nrmse_gaussian_norm <- average_nrmse(nrmse_gaussian_norm, "gaussian")

avg_nrmse_knn_prot <- average_nrmse(nrmse_knn_prot, "knn")
avg_nrmse_knn_phos <- average_nrmse(nrmse_knn_phos, "knn")
avg_nrmse_knn_norm <- average_nrmse(nrmse_knn_norm, "knn")

avg_nrmse_mle_prot <- average_nrmse(nrmse_mle_prot, "mle")
avg_nrmse_mle_phos <- average_nrmse(nrmse_mle_phos, "mle")
avg_nrmse_mle_norm <- average_nrmse(nrmse_mle_norm, "mle")

avg_nrmse_rf_prot <- average_nrmse(nrmse_rf_prot, "rf")
avg_nrmse_rf_phos <- average_nrmse(nrmse_rf_phos, "rf")
avg_nrmse_rf_norm <- average_nrmse(nrmse_rf_norm, "rf")

# Add a blank dataframe of NRMSE values for the without imputation method (we can't calculate NRMSE if there are NA values in the data)
avg_nrmse_without_prot <- data.frame(NRMSE=rep(NA,3), NRMSE_sd=rep(NA,3), row.names = rownames(without_metrics_prot))
avg_nrmse_without_phos <- data.frame(NRMSE=rep(NA,3), NRMSE_sd=rep(NA,3), row.names = rownames(without_metrics_phos))
avg_nrmse_without_norm <- data.frame(NRMSE=rep(NA,3), NRMSE_sd=rep(NA,3), row.names = rownames(without_metrics_norm))

# Bind all the NRMSE averages by datatype
avg_nrmse_prot <- base::rbind(avg_nrmse_without_prot, avg_nrmse_minmax_prot, avg_nrmse_gaussian_prot, avg_nrmse_knn_prot, avg_nrmse_mle_prot, avg_nrmse_rf_prot)
avg_nrmse_phos <- base::rbind(avg_nrmse_without_phos, avg_nrmse_minmax_phos, avg_nrmse_gaussian_phos, avg_nrmse_knn_phos, avg_nrmse_mle_phos, avg_nrmse_rf_phos)
avg_nrmse_norm <- base::rbind(avg_nrmse_without_norm, avg_nrmse_minmax_norm, avg_nrmse_gaussian_norm, avg_nrmse_knn_norm, avg_nrmse_mle_norm, avg_nrmse_rf_norm)

# Add those NRMSE values to the metrics table
metrics_overall_proteinGroups <- base::cbind(metrics_CM_proteinGroups, avg_nrmse_prot)
metrics_overall_Phosphosites <- base::cbind(metrics_CM_Phosphosites, avg_nrmse_phos)
metrics_overall_normedPhosphosites <- base::cbind(metrics_CM_normedPhosphosites, avg_nrmse_norm)

# Save those metrics in a txt table format
write.table(metrics_overall_proteinGroups, file="proteinGroups_avg_metrics.txt", sep="\t", col.names=NA)
write.table(metrics_overall_Phosphosites, file="Phospho (STY)Sites_avg_metrics.txt", sep="\t", col.names=NA)
write.table(metrics_overall_normedPhosphosites, file="Phospho (STY)Sites_normed_avg_metrics.txt", sep="\t", col.names=NA)

# You can get the barplots visualization of these tables in VisualizeMetrics.R