# This is an unofficial script for verification purposes, do not use for general application

# Goal: Compare full matrix (w/o MVs) and generated matrices histograms distributions with their metrics highlighted for comparisons

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

# Get all the respective distribution matrices
distribMV_prot <- read.table("proteinsMVdistribution.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)
distribMV_phos <- read.table("phosphositesMVdistribution.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)
distribMV_norm <- read.table("phosphosites_normedMVdistribution.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)

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

conditions_by_label_proteinGroups <- rearrange_matrix_label(conditions_table, labels_proteinGroups)
conditions_by_label_Phosphosites <- rearrange_matrix_label(conditions_table, labels_Phosphosites)
conditions_by_label_normedPhosphosites <-rearrange_matrix_label(conditions_table, labels_normedPhosphosites)

# Remove all rows that have at least one missing value to only get complete matrices
labels_proteinGroups <- na.omit(labels_proteinGroups)
conditions_by_label_proteinGroups <- na.omit(conditions_by_label_proteinGroups)

labels_Phosphosites <- na.omit(labels_Phosphosites)
conditions_by_label_Phosphosites <- na.omit(conditions_by_label_Phosphosites)

labels_normedPhosphosites <-na.omit(labels_normedPhosphosites)
conditions_by_label_normedPhosphosites <- na.omit(conditions_by_label_normedPhosphosites)

# Plot histograms per dataframe (once for global in blue, then with averaged missing values in red)
# Make a table of means, variance and skewness of those histograms too
library(parameters)

plot_histogram <- function(dataframes, measure_type, pdf_name){
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
      
      # Build histogram for the current category global distribution
      hist_df <- hist(mean_df, breaks="Scott", col="blue",
                      main=paste("Histogram of complete", new_category, "for", data_type[df_index]), xlab=paste(measure_type, new_category, "measure"))
      
      # Save the metrics of those histograms in the summary table (mean, variance, skewness)
      summary_table_df[summary_df_index,] <- c(mean(mean_df, na.rm=TRUE), var(mean_df, na.rm=TRUE), skewness(mean_df)[1])
      row.names(summary_table_df)[summary_df_index] <- paste0("Global_", new_category, "_", data_type[df_index])
      
      # Add 1 to the index of the summary table
      summary_df_index <- summary_df_index + 1
    }
  }
  
  return(summary_table_df)
}

# Create a list of dataframes by labels and conditions
df_by_labels_list <- list(labels_proteinGroups, labels_Phosphosites, labels_normedPhosphosites)
df_by_conditions_labels_list <- list(conditions_by_label_proteinGroups, conditions_by_label_Phosphosites, conditions_by_label_normedPhosphosites)

# Plot histograms for every possible dataset that is applicable and save the metrics in summary tables
summary_hist_labels <- plot_histogram(df_by_labels_list, "Ratio", "Histogram_by_labels_MV.pdf")
summary_hist_cond_labels <- plot_histogram(df_by_conditions_labels_list, "Ratio", "Histogram_by_conditions_labels_MV.pdf")

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

# Function to generate artificial missing values to every matrices from the given dataframe
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

missing_matrices_labels_prot <- apply_missing_matrices_to_dataframe(labels_proteinGroups, distribMV_prot, 20)
missing_matrices_labels_phos <- apply_missing_matrices_to_dataframe(labels_Phosphosites, distribMV_phos, 20)
missing_matrices_labels_norm <- apply_missing_matrices_to_dataframe(labels_normedPhosphosites, distribMV_norm, 20)

missing_matrices_conditions_ratios_prot <- apply_missing_matrices_to_dataframe(conditions_by_label_proteinGroups, distribMV_prot, 20)
missing_matrices_conditions_ratios_phos <- apply_missing_matrices_to_dataframe(conditions_by_label_Phosphosites, distribMV_phos, 20)
missing_matrices_conditions_ratios_norm <- apply_missing_matrices_to_dataframe(conditions_by_label_normedPhosphosites, distribMV_norm, 20)

# Function to plot histogram for a list with missing matrices

plot_histogram_generated_MV <- function(MV_matrices_list, data_type){
  # Initiate a dataframe table that summarizes basic metrics of histograms
  summary_table_df <- data.frame(Mean=numeric(), Variance=numeric(), Skewness=numeric(), stringsAsFactors=FALSE)
  
  # Set up the summary index as 1 initially (row indexing)
  summary_df_index <- 1
  
  # Iterate over the conditions
  for (i in 1:length(MV_matrices_list)){
    # Get the length of the list
    list_length <- length(MV_matrices_list[[i]])
    
    # Adjust mfrow depending on the length of the list
    col_mfrow <- ceiling(sqrt(list_length))
    row_mfrow <- ceiling(list_length / col_mfrow)
    par(mfrow=c(row_mfrow,col_mfrow), mar=rep(2,4))
    
    for (j in 1:list_length){
      df <- MV_matrices_list[[i]][[j]]
      
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
                        main=paste("Histogram n°", j,"of", new_category, "for", data_type), xlab=paste(new_category, "measure"))
        
        # Build histogram for missing values with the same bins from the global distribution
        hist(MV_df, breaks=hist_df$breaks, col="red", add=TRUE)
        
        # Save the metrics of those histograms in the summary table (mean, variance, skewness)
        summary_table_df[summary_df_index,] <- c(mean(mean_df, na.rm=TRUE), var(mean_df, na.rm=TRUE), skewness(mean_df)[1])
        row.names(summary_table_df)[summary_df_index] <- paste0("Global_", new_category, "_n°", j, "_", data_type)
        
        summary_table_df[summary_df_index + 1,] <- c(mean(MV_df, na.rm=TRUE), var(MV_df, na.rm=TRUE), skewness(MV_df)[1])
        row.names(summary_table_df)[summary_df_index + 1] <- paste0("MV_", new_category, "_n°", j, "_", data_type)
        
        # Add 1 to the index of the summary table
        summary_df_index <- summary_df_index + 2
      }
    }
  }
  return(summary_table_df)
}

summary_hist_generated_labels_prot <- plot_histogram_generated_MV(missing_matrices_labels_prot, "Protein Groups")
summary_hist_generated_labels_phos <- plot_histogram_generated_MV(missing_matrices_labels_phos, "Phosphosites")
summary_hist_generated_labels_norm <- plot_histogram_generated_MV(missing_matrices_labels_norm, "Normed ph - ft")

summary_hist_generated_conditions_prot <- plot_histogram_generated_MV(missing_matrices_conditions_ratios_prot, "Protein Groups")
summary_hist_generated_conditions_phos <- plot_histogram_generated_MV(missing_matrices_conditions_ratios_phos, "Phosphosites")
summary_hist_generated_conditions_norm <- plot_histogram_generated_MV(missing_matrices_conditions_ratios_norm, "Normed ph - ft")