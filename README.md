# Metabolomics Data Imputation

## Overview

This R script performs imputation of missing values in metabolomics data using the MICE (Multiple Imputation by Chained Equations) method. It's designed to handle metabolomics datasets in the format provided by Metabolon.  
For any question and issue contact the author Matteo Bolner.
## Summary of the process

1. The script reads the input files (in various formats: CSV, TSV, XLS, XLSX) and identifies exogenous (xenobiotic) metabolites.
2. Xenobiotic metabolites are excluded from the imputation process.
3. Correlations between all remaining metabolites are computed.
4. For each metabolite, the top 10 correlated **endogenous** metabolites are identified as predictors.
5. A prediction matrix is created based on these correlations and any additional predictors specified.
6. MICE imputation is performed using the specified method, generating the specified number of imputed datasets.

## Prerequisites

The script requires the following R packages:

- mice
- readxl
- dplyr
- tidyr
- tibble
- readr
- tools
- optparse

You can install these packages using the `install.packages()` function in R or simply running the `install_prerequisites.R` script available in the repository.

```
Rscript install_prerequisites.R 
```

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
- `-t, --imputed_only`: (Optional) Path to save only the imputed values (for testing purposes)
- `-s, --seed`: Random seed for MICE imputation (default: 42)
- `-m, --method`: Imputation method for MICE (default: "pmm")
- `-n, --number_imputations`: Number of imputed datasets to generate (default: 5)
- `-a, --additional_predictors`: Comma-separated list of additional non-metabolite columns to use as predictors

### Example Command

```
Rscript impute.R -d input_data.csv -c chemical_annotation.csv -o imputed_output.tsv -s 123 -m pmm -n 5 -a sex,weight
```

## Input File Requirements

1. Metabolomics Data File:
   - Should contain metabolite measurements with samples as rows and metabolites as columns
   - Can include additional metadata columns (e.g., sex, weight, sample ID)

2. Chemical Annotation File:
   - Must contain a column named "CHEM_ID" matching the metabolite IDs in the data file
   - Using purely numeric IDs for metabolites is not recommended, due to R's handling of numeric columns. You can handle this by prefixing them with "X" in both the column names of the data file and in the CHEM_ID column
   - Must include "SUPER_PATHWAY" column for identifying xenobiotics

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