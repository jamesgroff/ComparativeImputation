# Notes: Do it for each category again with the available data, use lme4 package in the end...
# Very much a WIP for mixed linear model, remind yourself with the mixed linear models on applied statistics course.
# A fix would be to ignore the mixed linear model entirely and only focus on the significance A method which yields a good amount of DE proteins/phosphosites (5-10% range).

# Other note: Make sure that the names of the condition labels are correct from the original label (more mish-mash in the overviewConditions table)
# Make it for condition intensities too and not just condition labels (optional)

# Goal: Pass the complete matrix to two statistical tests (mixed linear model and significance A test) for finding DE proteins/phosphosites

# First, set session or go to the directory where all the MaxQuant results are available (normally in combined/txt folder)

# General Setup of datasets

# Get all the datasets of proteins and phosphosites
proteinGroups <- read.table("proteinGroups.txt", header=TRUE, sep="\t", fill=TRUE, quote="")
Phosphosites <- read.table("Phospho (STY)Sites.txt", header=TRUE, sep="\t", fill=TRUE, quote="")

# If you have the Perseus combined dataset (after normalization by log2), you can also get the dataset for combined data
# It should already be filtered and ridden of undesirable rows.
normedPhosphosites <- read.table("Phospho (STY)Sites_normed.txt", header=TRUE, sep="\t", fill=TRUE, quote="")

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

removed_proteinGroups <- remove_rows_protein(proteinGroups)
removed_Phosphosites <- remove_rows_phos(Phosphosites)

# Labels bias

# Filter only normalized ratios as columns and only as separate experiments for labels grouping
desired_columns_lab <- names(proteinGroups)[grep("^Ratio.[H|M].[M|L].normalized.[A-z]+[0-9]+$", names(proteinGroups))]
labels_proteinGroups <- removed_proteinGroups[,desired_columns_lab]
labels_Phosphosites <- removed_Phosphosites[,desired_columns_lab]

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
    pair_combinations_matrix[i,] <- get_pair_combinations(i, conditions_matrix[i,])
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
    for (j in 1:num_conditions){
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
  
  # If it is already normed, apply another transformation where it's just a sign change
  if (deparse(substitute(numeric_matrix))=="labels_normedPhosphosites"){
    rearranged_matrix[,invert_indexes] <- (-1) * rearranged_matrix[,invert_indexes]
  }else{
    # Apply an inversion transformation to the columns if the invert vector is TRUE
    rearranged_matrix[,invert_indexes] <- 1 / rearranged_matrix[,invert_indexes]
  }
  # Reorder columns to have ordered conditions and replicates
  rearranged_matrix <- rearranged_matrix[, order(colnames(rearranged_matrix))]
  
  return(rearranged_matrix)
}

conditions_by_label_proteinGroups <- rearrange_matrix_label(conditions_table, labels_proteinGroups)
conditions_by_label_Phosphosites <- rearrange_matrix_label(conditions_table, labels_Phosphosites)
conditions_by_label_normedPhosphosites <- rearrange_matrix_label(conditions_table, labels_normedPhosphosites)

# Mixed linear model (inspired by Rudolf's code)

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
mixed_linear_model_significant <- function(complete_df) {
  # Transform the data in log2 base to apply log2 fold change except if it already is (normed phoshposites cases)
  if (deparse(substitute(complete_df))=="labels_normedPhosphosites" | deparse(substitute(complete_df))=="conditions_by_label_normedPhosphosites"){
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
    CI<-cbind(summary_df$CI_low, summary_df$CI_high)
    summary_df$significant <- apply(CI, 1, find_significant_lme)
    
    # Identify the number of significant proteins/phosphosites and store them as row ids
    row_ids_significant_lmer <- as.integer(row.names(summary_df[which(summary_df$significant=="+"),]))
    
    # Store the results for the current category
    result_list[[category]] <- row_ids_significant_lmer
  }
  
  return (result_list)
}

# Apply the mixed linear model to all the complete data available
list_significant_mlm_conditions_proteinGroups <- mixed_linear_model_significant(conditions_by_label_proteinGroups)
list_significant_mlm_conditions_Phosphosites <- mixed_linear_model_significant(conditions_by_label_Phosphosites)
list_significant_mlm_conditions_normedPhosphosites <- mixed_linear_model_significant(conditions_by_label_normedPhosphosites)

# Print one of the lists for verification
print(list_significant_mlm_conditions_proteinGroups)

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
outlier_test_significant <- function(complete_df){
  # Transform the data in log2 base to apply log2 fold change except if it already is (normed phoshposites cases)
  if (deparse(substitute(complete_df))=="labels_normedPhosphosites" | deparse(substitute(complete_df))=="conditions_by_label_normedPhosphosites"){
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

# Apply the significance A test to all the complete data available
list_significant_sigA_conditions_proteinGroups <- outlier_test_significant(conditions_by_label_proteinGroups)
list_significant_sigA_conditions_Phosphosites <- outlier_test_significant(conditions_by_label_Phosphosites)
list_significant_sigA_conditions_normedPhosphosites <- outlier_test_significant(conditions_by_label_normedPhosphosites)

# Print one of the lists for verification
print(list_significant_sigA_conditions_proteinGroups)

# Combine the two tests significant results to get DE proteins/phosphosites ids as vector or as matrix
# The idea is to get the two lists of mlm and sigA and get significant ids per conditions
compare_lists_ids <- function(list1, list2){
  # Find similarly named experiment variables that stores ids
  common_names <- intersect(names(list1), names(list2))
  
  # Initialize a list to store common variables
  common_list <- list()
  
  # Iterate over similarly named id variables
  for (name in common_names){
    # Convert id variables to sets
    set1 <- as.list(list1[[name]])
    set2 <- as.list(list2[[name]])
    
    # Find common elements
    common_ids <- intersect(set1, set2)
    
    # Store common id variable in common_list
    common_list[[name]] <- unlist(common_ids)
  }
  
  return(common_list)
}

# Apply the comparison to all significant lists that are comparable
full_significant_list_conditions_proteinGroups <- compare_lists_ids(list_significant_mlm_conditions_proteinGroups, list_significant_sigA_conditions_proteinGroups)
full_significant_list_conditions_Phosphosites <- compare_lists_ids(list_significant_mlm_conditions_Phosphosites, list_significant_sigA_conditions_Phosphosites)
full_significant_list_conditions_normedPhosphosites <- compare_lists_ids(list_significant_mlm_conditions_normedPhosphosites, list_significant_sigA_conditions_normedPhosphosites)

# Add columns to the original data to signify if the proteins/peptides are significantly DE depending on the conditions/label with a "+"

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
  colnames(result_df) <- paste("Significant", names(id_list))
  
  return(result_df)
}

# For instance, apply the resulting significant dataframes with all applicable data
significant_conditions_proteinGroups <- add_id_columns(proteinGroups, full_significant_list_conditions_proteinGroups)
significant_conditions_Phosphosites <- add_id_columns(Phosphosites, full_significant_list_conditions_Phosphosites)
significant_conditions_normedPhosphosites <- add_id_columns(normedPhosphosites, full_significant_list_conditions_normedPhosphosites)

# Save the significant DE proteins/phosphosites dataframes in txt format
write.table(significant_conditions_proteinGroups, file="proteinGroups_significant.txt", sep="\t", row.names=FALSE)
write.table(significant_conditions_Phosphosites, file="Phospho (STY)Sites_significant.txt", sep="\t", row.names=FALSE)
write.table(significant_conditions_normedPhosphosites, file="Phospho (STY)Sites_normed_significant.txt", sep="\t", row.names=FALSE)

# Those significant tables will be very useful for splitting the matrices in DE and Non-DE proteins/phosphosites
# They will also be the basis for the calculation of confusion matrices after imputing a certain strategy