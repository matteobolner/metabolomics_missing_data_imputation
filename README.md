# Metabolomics Data Imputation

## Overview

This R script performs imputation of missing values in metabolomics data using the MICE (Multiple Imputation by Chained Equations) method. It's designed to handle metabolomics datasets in the format provided by Metabolon.  
For any question and issue contact the author Matteo Bolner.
## Summary of the process

1. The script reads the input files (in various formats: CSV, TSV, XLS, XLSX)
2. Exogenous (Xenobiotic) metabolites are identified and excluded from the imputation process
3. Metabolites with too many missing values are removed from the dataset.
4. Individual metabolite values which are outliers (based on the interquartile range) are replaced with missing values
4. Correlations between all metabolites are computed.
5. For each metabolite, the top 10 correlated **endogenous** metabolites are identified as predictors.
6. A prediction matrix is created based on these correlations and any additional predictors specified.
7. MICE imputation is performed using the specified method, generating the specified number of imputed datasets.

## Prerequisites

The script was tested and works on the latest R version (4.3.3); the minimum version compatilble should be 3.6.0, due to package requirements; previous versions may not work.

The script requires the following R packages:

- mice
- readxl
- dplyr
- tidyr
- tibble
- readr
- tools
- optparse

You can install these packages individually using the `install.packages()` function in R or simply running the `install_prerequisites.R` script available in the repository.

```
Rscript install_prerequisites.R 
```

## Input File Requirements

**!!** Filtering based on missing value contents and outliers is performed by this script before imputation. We recommend removing metabolites with more than 25% of missing values across samples, and outliers over 5 times the interquartile range.  
**!!** The core of the script involves matching information between metabolite information and metabolite abundance in the samples. To do so, unique metabolite IDs (one per metabolite) are necessary, and they must be consistent between the data and the chemical annotation files (see below).

1. Metabolomics Data File:
   - Should contain metabolite measurements with samples as rows and metabolites as columns
   - Can include additional metadata columns (e.g., sex, weight, sample ID); they will be ignored in the imputation process, unless they are explicitly set as predictors
   - The column names for metabolite abundance columns should be named consistently with the CHEM_ID column of the chemical annotation file: the column order doesn't matter, since the script automatically extracts them from the dataset, but the names must be identical for the same metabolite
2. Chemical Annotation File:
   - Must contain a column named "CHEM_ID" matching the metabolite IDs in the data file
   - Using purely numeric IDs for metabolites is not recommended, due to R's handling of numeric columns, which automatically adds an "X" before the number. You can handle this by prefixing them yourself with "X" in both the column names of the data file and in the CHEM_ID column.
   - Must include "SUPER_PATHWAY" column for identifying xenobiotic, endogenous and unknown metabolites



## Usage

The script is designed to be run from the command line. Here's the basic usage:

```
Rscript impute.R [options]
```

You can also set it as executable and run it directly:
```
chmod +x impute.R 
./impute.R [options]
```


### Command-line Options

- `-d, --data`: Input file containing the metabolomics data (required)
- `-c, --chemical_annotation`: Input file containing the metabolite chemical annotation (required)
- `-o, --output`: Path to the file where the imputed datasets will be saved in TSV format (required)
- `-r, --remove_missing_over_threshold`: Percentage of missing values under which metabolites will be discarded. For example, `-r 0.25` will discard metabolites with more than 25% of values missing. Default value is 0.25. 
- `-u, --remove_outliers`: Threshold for the interquartile range, over which the individual metabolite values will be considered an outlier, and replaced with NA. The outlier is computed with regards to the whole column. Default and recommended value is 5
- `-l, --missing_outlier_stats`: (Optional) Path to a file where statistics about missing values and outliers per metabolite will be saved.
- `-t, --imputed_only`: (Optional) Path to save only the imputed values (for testing purposes)
- `-s, --seed`: Random seed for MICE imputation (default: 42)
- `-m, --method`: Imputation method for MICE (default: "pmm")
- `-n, --number_imputations`: Number of imputed datasets to generate (default: 5)
- `-a, --additional_predictors`: Comma-separated list of additional non-metabolite columns to use as predictors



### Example Command

```
Rscript impute.R -d input_data.csv -c chemical_annotation.csv -o imputed_output.tsv -s 42 -m pmm -n 5 -r 0.25 -u 5 -a sex,weight -l stats.tsv
```

## Output

The script generates two types of output:

1. Complete Imputed Datasets:
   - Saved in the file specified by the `-o` option
   - Contains all imputed datasets in long format

2. Imputed Values Only (optional):
   - Saved in the file specified by the `-t` option
   - Contains only the imputed values, useful for testing and verification

Both outputs are in TSV format.

## Notes

- The default imputation method is "pmm" (Predictive Mean Matching), which is recommended for most use cases.
- The script uses a seed for reproducibility. Change the seed value to generate different imputation results.