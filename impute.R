#!/usr/bin/Rscript

#LOAD THE REQUIRED LIBRARIES
library(mice, warn.conflicts = FALSE)
library(readxl)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(tibble)
library(readr)
library(tools)
library(optparse)


####################################################################################
###PARSE COMMAND-LINE ARGUMENTS
####################################################################################
option_list <- list(
  make_option(c("-d", "--data"), action = "store",
              help = "Input file containing the metabolomics data in {tsv,csv,xlsx,xls} format", ),
  make_option(
    c("-c", "--chemical_annotation"),
    action = "store",
    help = "Input file containing the metabolite chemical annotation in {tsv,csv,xlsx,xls} format"
  ),
  make_option(c("-o", "--output"), action = "store",
              help = "Path to the file where the imputed datasets will be saved, in .tsv format"),
  make_option(
    c("-t", "--imputed_only"),
    action = "store",
    default = NULL,
    help = "For testing purposes only: path to the file where only the imputed values will be saved, in .tsv format"
  ),
  make_option(
    c("-s", "--seed"),
    type = "integer",
    default = 42,
    help = "Random seed for MICE imputation",
    metavar = "integer number"
  ),
  make_option(
    c("-m", "--method"),
    type = "character",
    default = "pmm",
    help = "Imputation method for MICE; the default and recommended for most use-cases is 'pmm', predictive mean matching. Other examples are rf | logreg | polyreg | polr",
  ),
  make_option(
    c("-n", "--number_imputations"),
    type = "integer",
    default = 5,
    help = "Number of imputed datasets to generate",
  ),
  make_option(
    c("-a", "--additional_predictors"),
    type = "character",
    default = NULL,
    help = "Comma-separated list of additional non-metabolite columns of the dataset to use as predictors, such as sex,weight etc",
    metavar = "col1,col2"
  ),
  make_option(
    c("-r", "--remove_missing_over_threshold"),
    type = "numeric",
    default = 0.25,
    help = "Remove metabolites from the dataset if their percentage of missing values is above this threshold. Ranges from 0 to 1; default and recommended value is 0.25",
    metavar = "Range from 0 to 1"
  ),
  make_option(
    c("-u", "--remove_outliers"),
    type = "numeric",
    default = 1.5,
    help = "Remove outlier metabolite values from the dataset if their value is outside the interquartile range by X times, with X being the threshold set. Default and recommended value for the threshold is 1.5",
    metavar = "numeric"
  ),
  make_option(
    c("-l", "--missing_outlier_stats"),
    action = "store",
    default = NULL,
    help = "Optional: Path to save a table with statistics about missing and outlier values in the metabolites",
  )
)



opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$data)) {
  stop("Input file path is mandatory. Use -d or --data to specify it.")
}
if (is.null(opt$chemical_annotation)) {
  stop("Input file path is mandatory. Use -c or --chemical-annotion to specify it.")
}
if (is.null(opt$output)) {
  stop("Output file path is mandatory. Use -o or --output to specify it.")
}


#Define arguments to use in the code

n_imputations = opt$number_imputations
mice_random_seed = opt$seed
set.seed(mice_random_seed)
method = opt$method
additional_predictors = opt$additional_predictors
output_path = opt$output
imputed_only_path = opt$imputed_only

if (length(additional_predictors) > 0) {
  additional_predictors = strsplit(additional_predictors, ",")[[1]]
  additional_predictors <-
    unlist(lapply(additional_predictors, trimws))
}


####################################################################################
###READ INPUT FILES
####################################################################################

read_file <- function(file_path) {
  # Check file extension and read accordingly
  file_ext <- file_ext(file_path)
  
  if (file_ext == "csv") {
    df <- read_csv(file_path)
  } else if (file_ext == "tsv") {
    df <- read_tsv(file_path)
  } else if (file_ext %in% c("xls", "xlsx")) {
    df <- read_excel(file_path)
  } else {
    stop("Unsupported file type")
  }
  
  return(df)
}

#read main dataset
cat("Reading data...\n\n")
df = read_file(opt$data)
cat("\n\n\n")
#read metabolite annotation
cat("Reading chemical annotation data...\n\n")
metabolites = read_file(opt$chemical_annotation)
cat("\n\n\n")



####################################################################################
###REMOVE METABOLITES WITH TOO MANY MISSING VALUES
####################################################################################


threshold <- opt$remove_missing_over_threshold

cat("Removing metabolites with too many missing data...\n\n")

df <- df %>%
  select(where( ~ mean(is.na(.)) <= threshold))

####################################################################################
###IDENTIFY UNNAMED AND XENOBIOTIC METABOLITES
####################################################################################

assign_annotation <- function(df) {
  df %>%
    mutate(
      ORIGIN = case_when(
        is.na(SUPER_PATHWAY) ~ "UNNAMED",
        SUPER_PATHWAY == "Xenobiotics" ~ "EXOGENOUS",
        TRUE ~ "ENDOGENOUS"
      )
    )
}

metabolites = assign_annotation(metabolites)

#drop exogenous metabolites

metabolites <- metabolites %>%
  filter(CHEM_ID %in% colnames(df))

sample_metadata <- df %>%
  select(-one_of(metabolites$CHEM_ID))

exogenous_metabolites = metabolites %>%
  filter(ORIGIN == "EXOGENOUS") %>%
  pull("CHEM_ID")

endogenous_metabolite_ids <- metabolites %>%
  filter(ORIGIN == "ENDOGENOUS") %>%
  pull("CHEM_ID")

metabolites <- metabolites %>%
  filter(ORIGIN != "EXOGENOUS")

metabolite_ids = c(metabolites$CHEM_ID)


####################################################################################
###REMOVE OUTLIER VALUES FROM METABOLITES
####################################################################################

process_outliers <- function(data, column_names, threshold = 1.5) {
  # Initialize dataframe for storing counts
  count_summary <- data.frame(
    CHEM_ID = character(),
    missing_count = integer(),
    outlier_count = integer(),
    stringsAsFactors = FALSE
  )
  
  processed_data <- data %>%
    mutate(across(all_of(column_names), function(col) {
      # Count missing values
      missing_count <- sum(is.na(col))
      
      # Calculate statistics for non-missing values
      non_missing <- col[!is.na(col)]
      median_val <- median(non_missing)
      q1 <- quantile(non_missing, 0.25)
      q3 <- quantile(non_missing, 0.75)
      iqr <- q3 - q1
      cutoff_lower <- median_val - (threshold * iqr)
      cutoff_upper <- median_val + (threshold * iqr)
      
      # Identify outliers
      is_outlier <-
        (non_missing < cutoff_lower) | (non_missing > cutoff_upper)
      outlier_count <- sum(is_outlier)
      
      # Add counts to summary dataframe
      count_summary <<- count_summary %>%
        add_row(
          CHEM_ID = cur_column(),
          missing_count = missing_count,
          outlier_count = outlier_count
        )
      
      # Replace outliers with NaN
      col[!is.na(col)] <- if_else(is_outlier, NaN, non_missing)
      col
    }))
  
  list(processed_data = processed_data,
       count_summary = count_summary)
}

cat("Removing outlier values from the dataset...\n\n")

processed_outliers <-
  process_outliers(df, metabolite_ids, opt$remove_outliers)
df <- processed_outliers$processed_data
missing_outlier_summary <- processed_outliers$count_summary


if (is.null(opt$missing_outlier_stats))
{
  cat("")
} else
{
  write.table(
    missing_outlier_summary,
    file = opt$missing_outlier_stats,
    sep = '\t',
    row.names = FALSE
  )
}


#extract metabolite data only
metabolite_data <- df %>% select(any_of(metabolite_ids))
metabolite_data = as.data.frame(lapply(metabolite_data, as.numeric))


####################################################################################
###COMPUTE CORRELATIONS BETWEEN METABOLITES AND REMOVE UNNAMED FROM TOP 10
####################################################################################

correlations <-
  cor(metabolite_data, method = 'pearson', use = 'pairwise.complete.obs')
correlations <- as_tibble(correlations)
correlations$CHEM_ID <- colnames(correlations)
correlations = column_to_rownames(correlations, var = "CHEM_ID")
correlations_endogenous_only <- correlations %>%
  rownames_to_column(var = "CHEM_ID") %>%
  filter(CHEM_ID %in% endogenous_metabolite_ids) %>%
  column_to_rownames(var = "CHEM_ID")

top_10_correlated <-
  lapply(colnames(correlations_endogenous_only), function(col) {
    # Extract the correlation values for the current column
    cor_vals <- correlations_endogenous_only[[col]]
    
    # Remove self-correlation (correlation of the column with itself)
    cor_vals <-
      cor_vals[colnames(correlations_endogenous_only) != col]
    
    # Get the names of the remaining columns
    remaining_cols <-
      colnames(correlations_endogenous_only)[colnames(correlations_endogenous_only) != col]
    
    # Return the names of the columns with the top 10 correlations
    return(remaining_cols[order(-abs(cor_vals))[1:10]])
  })

names(top_10_correlated) = colnames(correlations_endogenous_only)


####################################################################################
###MICE IMPUTATION
####################################################################################

#build where matrix to identify missing values

where_matrix <- make.where(df)
where_matrix[, c(colnames(sample_metadata))] <- FALSE

#build default prediction matrix and set all predictors to 0

pred_matrix <- make.predictorMatrix(df)
pred_matrix[, colnames(pred_matrix)] <- 0

#update prediction matrix with top 10 correlated metabolites for each metabolite
update_pred_matrix <- function(mat, predictor_list) {
  for (row_name in rownames(mat)) {
    names_to_set <- predictor_list[[row_name]]
    mat[row_name, ] <- ifelse(colnames(mat) %in% names_to_set, 1, 0)
  }
  return(mat)
}

pred_matrix <- update_pred_matrix(pred_matrix, top_10_correlated)

#if additional predictors were specified in the arguments, set their status as 1 in the prediction matrix
if (length(additional_predictors) > 0) {
  pred_matrix[, additional_predictors] <- 1
}

#impute missing values
cat("Starting imputation...\n\n")

IMP <-
  mice(
    df,
    m = n_imputations,
    where = where_matrix ,
    predictorMatrix = pred_matrix,
    method = method,
    seed = mice_random_seed
  )
cat("Imputation finished! \n")

####################################################################################
###SAVE OUTPUT FILES
####################################################################################
cat("Saving output... \n")

#extract all imputed datasets
output <-
  complete(IMP ,  action = 'long' , include = TRUE)

write.table(output, file = output_path, sep = '\t')

#extract only imputed values (useful for testing)

if (is.null(opt$imputed_only))
{
  cat("")
} else
{
  only_imputed_values <- IMP$imp
  all_imputed_values <-
    bind_rows(only_imputed_values, .id = "metabolite")
  write.table(all_imputed_values, file = imputed_only_path, sep = '\t')
  
}
cat("Saved! Quitting... \n")
