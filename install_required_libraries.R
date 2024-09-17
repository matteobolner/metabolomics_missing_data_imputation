#!/usr/bin/Rscript

# List of required packages
packages <- c("mice",
              "readxl",
              "dplyr",
              "tidyr",
              "tibble",
              "readr",
              "tools",
              "optparse")

# Function to install packages if they're not already installed
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# Main installation process
main <- function() {
  cat("Starting installation of required packages...\n")
  
  # Install and load each package
  for (package in packages) {
    tryCatch({
      install_if_missing(package)
      cat(sprintf(
        "Package '%s' installed and loaded successfully.\n",
        package
      ))
    }, error = function(e) {
      cat(sprintf("Error installing package '%s': %s\n", package, e$message))
    })
  }
  
  cat("\nInstallation process completed. You may now run impute.R\n")
}

main()