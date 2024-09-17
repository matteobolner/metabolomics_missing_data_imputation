#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("tibble")
#install.packages("mice")
#install.packages("readxl")
#install.packages("optparse")

library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(tools)
library(readxl)
library(optparse)
library(mice)

#setwd("work/workspace/pigphenomics_metabolomics/")
data_path="data/raw_data/data_filtered.tsv"
chemical_annotation_path="config/metabolites.tsv"


mice_random_seed=1000
mice_imputation_method='pmm'
metadata_predictors='sex, weight'
n_imputations=5



option_list <- list(
  make_option(c("-d", "--data"), action="store",
              help="Path to the file containing the metabolomics data"),
  make_option(c("-c", "--chemical_annotation"), action="store",
              help="Path to the file containing the metabolite chemical annotation"),
  make_option(c("-s", "--seed"), type="integer", default=42,
              help="Random seed for MICE imputation",
              metavar="integer number"),
  make_option(c("-m", "--method"), type="string", default="pmm",
              help="Imputation method for MICE the default and recommended for most use-cases is pmm",
              metavar="pmm | rf | logreg | polyreg | polr"),
  make_option(c("-n", "--number_imputations"), type="string", default=5,
              help="Number of imputed datasets to generate",
              metavar="pmm | rf | logreg | polyreg | polr"),
    make_option(c("-ap", "--additional-predictors"), type="string", default="aaa",
              help="Comma-separated list of additional non-metabolite columns of the dataset to use as predictors, such as sex,weight etc",
              metavar="sex,weght")
  
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)



n_imputations=as.numeric(n_imputations)
mice_random_seed=as.numeric(mice_random_seed)
set.seed(mice_random_seed)
method=mice_imputation_method
metadata_predictors=strsplit(metadata_predictors, ",")[[1]]
metadata_predictors <- unlist(lapply(metadata_predictors, trimws))

cli_alert_success("Updated database.")

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
df=read_file(data_path)

#read metabolite annotation
metabolites=read_file(chemical_annotation_path)


### IDENTIFY UNNAMED AND XENOBIOTICS
assign_annotation <- function(df) {
  df %>%
    mutate(ORIGIN = case_when(
      is.na(SUPER_PATHWAY) ~ "UNNAMED",
      SUPER_PATHWAY == "Xenobiotics" ~ "EXOGENOUS",
      TRUE ~ "ENDOGENOUS"
    ))
}

metabolites=assign_annotation(metabolites)

#drop exogenous metabolites

metabolites <- metabolites %>%
  filter(CHEM_ID %in% colnames(df))

sample_metadata <- df %>%
  select(-one_of(metabolites$CHEM_ID))

exogenous_metabolites=metabolites %>%
  filter(ORIGIN == "EXOGENOUS") %>%
  pull("CHEM_ID")

endogenous_metabolite_ids <- metabolites %>%
  filter(ORIGIN == "ENDOGENOUS") %>%
  pull("CHEM_ID")

metabolites <- metabolites %>%
  filter(ORIGIN != "EXOGENOUS")

metabolite_ids=c(metabolites$CHEM_ID)

metabolite_data <- df %>% select(any_of(metabolite_ids))
metabolite_data=as.data.frame(lapply(metabolite_data, as.numeric))

###COMPUTE CORRELATIONS AND REMOVE UNNAMED FROM TOP 10
correlations <- cor(metabolite_data, method='pearson', use='pairwise.complete.obs')
correlations <- as_tibble(correlations)
correlations$CHEM_ID <- colnames(correlations)
correlations=column_to_rownames(correlations, var = "CHEM_ID")
correlations_endogenous_only <- correlations %>%
  rownames_to_column(var = "CHEM_ID") %>%
  filter(CHEM_ID %in% endogenous_metabolite_ids) %>%
  column_to_rownames(var = "CHEM_ID")

top_10_correlated <- lapply(colnames(correlations_endogenous_only), function(col) {
  # Extract the correlation values for the current column
  cor_vals <- correlations_endogenous_only[[col]]
  
  # Remove self-correlation (correlation of the column with itself)
  cor_vals <- cor_vals[colnames(correlations_endogenous_only) != col]
  
  # Get the names of the remaining columns
  remaining_cols <- colnames(correlations_endogenous_only)[colnames(correlations_endogenous_only) != col]
  
  # Return the names of the columns with the top 10 correlations
  return(remaining_cols[order(-abs(cor_vals))[1:10]])
})

names(top_10_correlated)=colnames(correlations_endogenous_only)

### MICE IMPUTATION

### BUILD WHERE MATRIX

where_matrix <- make.where(df)
where_matrix[,c(colnames(sample_metadata))] <- FALSE

### BUILD PREDICTION MATRIX

pred_matrix <- make.predictorMatrix(df)
pred_matrix[, colnames(pred_matrix)] <- 0

update_pred_matrix <- function(mat, predictor_list) {
  for (row_name in rownames(mat)) {
    names_to_set <- predictor_list[[row_name]]
    mat[row_name, ] <- ifelse(colnames(mat) %in% names_to_set, 1, 0)
  }
  return(mat)
}

pred_matrix <- update_pred_matrix(pred_matrix, top_10_correlated)

if (length(metadata_predictors)>0) {
pred_matrix[, metadata_predictors] <- 1
}

###MICE IMPUTATION

IMP <- mice(df, m=1, where = where_matrix , predictorMatrix = pred_matrix, method=mice_imputation_method, seed=mice_random_seed)
output <- complete(IMP ,  action = 'long' , include = TRUE)#[c(".imp",metabolite)]

###SAVE OUTPUTS

only_imputed_values <- IMP$imp
all_imputed_values <- bind_rows(only_imputed_values, .id = "metabolite")

write.table(output, file=output_path, sep='\t')
write.table(all_imputed_values, file=only_imputed_values_path, sep='\t')


