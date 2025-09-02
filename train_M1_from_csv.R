#!/usr/bin/env Rscript

#
# M1 Hepatotoxicity Model Training - Reading from CSV Files
# Modified version that reads compound data from CSV instead of hardcoded values
#

library(reticulate)
# Update this path to your conda environment
use_python("/Users/luqmanawoniyi_1/miniconda3/envs/rdkit-env/bin/python", required = TRUE)

cat("\n=======================================================\n")
cat("M1 Hepatotoxicity Model Training from CSV\n")
cat("Enhanced Descriptor Calculator + GA-SVM Framework\n")
cat("=======================================================\n\n")

# Configuration options
APPLY_QSAR_FILTERS <- TRUE      # Filter problematic molecules
APPLY_STAT_FILTERS <- FALSE     # Filter descriptors statistically
VERBOSE_OUTPUT <- TRUE          # Show detailed progress

# CSV input configuration
INPUT_CSV_FILE <- "datasets/compounds_data.csv"  # Default input file name
OUTPUT_PREFIX <- "M1_hepatotoxicity"              # Output file prefix

# Check command line arguments for custom CSV file
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  # If path doesn't start with datasets/, prepend it
  input_file <- args[1]
  if (!startsWith(input_file, "datasets/")) {
    input_file <- paste0("datasets/", input_file)
  }
  INPUT_CSV_FILE <- input_file
  cat(paste("Using input file:", INPUT_CSV_FILE, "\n"))
}

if (length(args) > 1) {
  OUTPUT_PREFIX <- args[2]
  cat(paste("Using output prefix:", OUTPUT_PREFIX, "\n"))
}

cat("\nConfiguration:\n")
cat(paste("- Input CSV file:", INPUT_CSV_FILE, "\n"))
cat(paste("- Output prefix:", OUTPUT_PREFIX, "\n"))
cat(paste("- QSAR molecular filters:", ifelse(APPLY_QSAR_FILTERS, "ENABLED", "DISABLED"), "\n"))
cat(paste("- Statistical descriptor filters:", ifelse(APPLY_STAT_FILTERS, "ENABLED", "DISABLED"), "\n"))
cat(paste("- Verbose output:", ifelse(VERBOSE_OUTPUT, "ENABLED", "DISABLED"), "\n\n"))

# Step 1: Prepare GA-SVM scripts
cat("Step 1: Preparing GA-SVM scripts...\n")

# Load GA-SVM functions if needed
if (!file.exists("tx5b00465_si_002/GA-SVM_v2.5_main_script.r")) {
  stop("Error: GA-SVM scripts not found in tx5b00465_si_002/")
}

# Step 2: Load compound data from CSV
cat(paste("\nStep 2: Loading compound data from", INPUT_CSV_FILE, "...\n"))

# Check if file exists
# Create datasets directory if it doesn't exist
if (!dir.exists("datasets")) {
  dir.create("datasets")
  cat("Created datasets/ directory\n")
}

if (!file.exists(INPUT_CSV_FILE)) {
  cat(paste("\nError: Input file not found:", INPUT_CSV_FILE, "\n"))
  cat("Please create a CSV file in the datasets/ folder with the following format:\n")
  cat("ID,SMILES,Hepatotoxic\n")
  cat("COMPOUND1,CCO,0\n")
  cat("COMPOUND2,CC(=O)NC1=CC=C(C=C1)O,1\n")
  cat("...\n\n")
  cat("Where:\n")
  cat("- ID: Unique compound identifier\n")
  cat("- SMILES: SMILES string of the compound\n")
  cat("- Hepatotoxic: 1 for hepatotoxic, 0 for non-hepatotoxic\n\n")
  
  # Create example CSV file in datasets folder
  example_file <- "datasets/compounds_data_example.csv"
  cat(paste("Creating example file:", example_file, "\n"))
  example_data <- data.frame(
    ID = c("ACETAMINOPHEN", "CAFFEINE", "TROGLITAZONE", "ASPIRIN"),
    SMILES = c(
      "CC(=O)NC1=CC=C(C=C1)O",
      "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
      "Cc1c(C)c2c(c(C)c1O)CCC(C)(COc1ccc(CC3SC(=O)NC3=O)cc1)O2",
      "CC(=O)OC1=CC=CC=C1C(=O)O"
    ),
    Hepatotoxic = c(1, 0, 1, 0),
    stringsAsFactors = FALSE
  )
  write.csv(example_data, example_file, row.names = FALSE)
  
  stop("Please prepare your data in the datasets/ folder and run again.")
}

# Read the CSV file
compounds <- read.csv(INPUT_CSV_FILE, stringsAsFactors = FALSE)

# Validate required columns
required_cols <- c("ID", "SMILES", "Hepatotoxic")
missing_cols <- setdiff(required_cols, names(compounds))
if (length(missing_cols) > 0) {
  stop(paste("Error: Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# Validate data
if (any(is.na(compounds$ID)) || any(compounds$ID == "")) {
  stop("Error: All compounds must have a valid ID")
}

if (any(is.na(compounds$SMILES)) || any(compounds$SMILES == "")) {
  stop("Error: All compounds must have a valid SMILES string")
}

if (any(!compounds$Hepatotoxic %in% c(0, 1))) {
  stop("Error: Hepatotoxic column must contain only 0 or 1 values")
}

cat(paste("Loaded", nrow(compounds), "compounds from CSV\n"))
cat(paste("- Hepatotoxic:", sum(compounds$Hepatotoxic), "compounds\n"))
cat(paste("- Non-hepatotoxic:", sum(1 - compounds$Hepatotoxic), "compounds\n"))

# Check class balance
toxic_ratio <- sum(compounds$Hepatotoxic) / nrow(compounds)
if (toxic_ratio < 0.2 || toxic_ratio > 0.8) {
  cat("\nWarning: Dataset is highly imbalanced!\n")
  cat(sprintf("Toxic ratio: %.1f%%\n", toxic_ratio * 100))
  cat("Consider balancing your dataset for better model performance.\n\n")
}

# Step 3: Calculate descriptors with optional filtering
cat("\nStep 3: Calculating molecular descriptors with optional filtering...\n")

# Import Python descriptor calculator
source_python("hepatotox_descriptors_improved.py")

# Initialize improved calculator
calculator <- ImprovedHepatotoxicityDescriptors()

# Calculate descriptors with optional molecular filtering
all_descriptors <- calculator$calculate_all_descriptors(
  compounds$SMILES, 
  preprocess = TRUE,  # Always do preprocessing
  apply_qsar_filters = APPLY_QSAR_FILTERS
)

# Add compound IDs
all_descriptors <- cbind(ID = compounds$ID, all_descriptors[, -1])  # Remove SMILES column

original_features <- ncol(all_descriptors) - 1
cat(paste("Calculated", original_features, "descriptors for", nrow(all_descriptors), "compounds\n"))

# Optional: Apply statistical descriptor filtering
if (APPLY_STAT_FILTERS) {
  cat("\nApplying statistical descriptor filters...\n")
  
  filtered_result <- calculator$filter_descriptors(
    all_descriptors, 
    apply_statistical_filters = TRUE
  )
  
  filtered_descriptors <- filtered_result[[1]]
  filter_stats <- filtered_result[[2]]
  
  if (filter_stats$filters_applied) {
    filtered_features <- ncol(filtered_descriptors) - 1
    cat(paste("Statistical filtering: kept", filtered_features, "out of", original_features, "descriptors\n"))
    
    if ("constant_removed" %in% names(filter_stats)) {
      cat(paste("- Removed", length(filter_stats$constant_removed), "constant features\n"))
    }
    if ("variance_removed" %in% names(filter_stats)) {
      cat(paste("- Removed", length(filter_stats$variance_removed), "low-variance features\n"))
    }
    if ("correlation_removed" %in% names(filter_stats)) {
      cat(paste("- Removed", length(filter_stats$correlation_removed), "highly correlated features\n"))
    }
    
    # Use filtered descriptors
    all_descriptors <- filtered_descriptors
  }
}

# Step 4: Prepare GA-SVM input files
cat("\nStep 4: Preparing input files for GA-SVM...\n")

# Remove rows with all NA descriptors (failed molecules)
non_na_rows <- !apply(is.na(all_descriptors[,-1]), 1, all)  # Exclude ID column from check
all_descriptors_clean <- all_descriptors[non_na_rows, ]
compounds_clean <- compounds[non_na_rows, ]

cat(sprintf("\nRemoved %d molecules with all NA descriptors\n", 
            nrow(all_descriptors) - nrow(all_descriptors_clean)))
cat(sprintf("Final dataset: %d compounds\n", nrow(all_descriptors_clean)))

# Create datasets directory if it doesn't exist
if (!dir.exists("datasets")) {
  dir.create("datasets")
}

# Save descriptor file
desc_file <- paste0("datasets/", OUTPUT_PREFIX, "_descriptors.desc")
write.table(all_descriptors_clean, desc_file, 
            sep = " ", row.names = FALSE, quote = FALSE)
cat(paste("Saved descriptors to:", desc_file, "\n"))

# Save toxicity file
tox_data <- data.frame(
  ID = compounds_clean$ID,
  Hepatotoxic_Human = compounds_clean$Hepatotoxic
)
tox_file <- paste0("datasets/", OUTPUT_PREFIX, "_toxicity.txt")
write.table(tox_data, tox_file, 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Saved toxicity data to:", tox_file, "\n"))

# Create stratified train/test/validation split (using cleaned data)
set.seed(42)
toxic_idx <- which(compounds_clean$Hepatotoxic == 1)
safe_idx <- which(compounds_clean$Hepatotoxic == 0)

# Ensure we have enough compounds for splitting
min_per_class <- 3  # Minimum compounds per class for train/test/val split
if (length(toxic_idx) < min_per_class || length(safe_idx) < min_per_class) {
  cat("\nWarning: Not enough compounds for proper train/test/validation split!\n")
  cat(sprintf("Toxic compounds: %d, Safe compounds: %d\n", length(toxic_idx), length(safe_idx)))
  cat("Using simple 70/20/10 split without stratification.\n\n")
  
  # Simple random split
  n <- nrow(compounds_clean)
  train_size <- floor(0.7 * n)
  test_size <- floor(0.2 * n)
  
  idx <- sample(1:n)
  train_idx <- idx[1:train_size]
  test_idx <- idx[(train_size + 1):(train_size + test_size)]
  val_idx <- idx[(train_size + test_size + 1):n]
} else {
  # Stratified sampling
  train_toxic <- sample(toxic_idx, floor(0.7 * length(toxic_idx)))
  test_toxic <- sample(setdiff(toxic_idx, train_toxic), floor(0.2 * length(toxic_idx)))
  val_toxic <- setdiff(toxic_idx, c(train_toxic, test_toxic))

  train_safe <- sample(safe_idx, floor(0.7 * length(safe_idx)))
  test_safe <- sample(setdiff(safe_idx, train_safe), floor(0.2 * length(safe_idx)))
  val_safe <- setdiff(safe_idx, c(train_safe, test_safe))

  train_idx <- c(train_toxic, train_safe)
  test_idx <- c(test_toxic, test_safe)
  val_idx <- c(val_toxic, val_safe)
}

# Create split assignment
split_data <- data.frame(
  ID = compounds_clean$ID,
  SET_hepatotoxicity = character(nrow(compounds_clean)),
  stringsAsFactors = FALSE
)

split_data$SET_hepatotoxicity[train_idx] <- "train"
split_data$SET_hepatotoxicity[test_idx] <- "test"
split_data$SET_hepatotoxicity[val_idx] <- "val"

# Save split file
split_file <- paste0("datasets/", OUTPUT_PREFIX, "_split.txt")
write.table(split_data, split_file, 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Saved train/test/val split to:", split_file, "\n"))

cat("\nDataset split summary:\n")
cat(paste("- Training set:", sum(split_data$SET_hepatotoxicity == "train"), "compounds\n"))
cat(paste("- Test set:", sum(split_data$SET_hepatotoxicity == "test"), "compounds\n"))
cat(paste("- Validation set:", sum(split_data$SET_hepatotoxicity == "val"), "compounds\n"))

# Step 5: Create GA-SVM configuration script
cat("\nStep 5: Creating GA-SVM configuration...\n")

config_script <- paste0('#!/usr/bin/env Rscript
#
# GA-SVM Configuration for ', OUTPUT_PREFIX, '
# Generated from CSV input with filtering options: QSAR=', APPLY_QSAR_FILTERS, ', STAT=', APPLY_STAT_FILTERS, '
#

cat("\\n=== ', OUTPUT_PREFIX, ' Model Training ===\\n")
cat("Using Enhanced Descriptor Calculator + Filtering\\n")
cat("Filtering applied: QSAR=', APPLY_QSAR_FILTERS, ', Statistical=', APPLY_STAT_FILTERS, '\\n\\n")

## General settings
DEBUG = FALSE
FILENAME = "', OUTPUT_PREFIX, '_model.r"

## Descriptor files
DESCRIPTOR.FILES = "datasets/', OUTPUT_PREFIX, '_descriptors"
DESCRIPTOR.postfixes = "desc"

## Toxicity data
TOXICITY.FILES = "datasets/', OUTPUT_PREFIX, '_toxicity.txt"
ENDPOINT = "Hepatotoxic_Human"
NUMERIC.ENDPOINT = FALSE

## Train/Test/Validation split
TRAIN.TEST.SPLIT.FILES = "datasets/', OUTPUT_PREFIX, '_split.txt"
SET.column.header = "SET_hepatotoxicity"

## SVM Settings (from Mulliner et al. 2016)
svm.type = "nu-classification"
kernel.type = "radial"
tolerance = 0.01
epsilon = 0.1
NU = 0.7
GAMMA = 0.01

## Genetic Algorithm settings
pS = 50           # Population size
iter = 20         # Iterations per run
Nruns = 10        # Number of runs
N.desc = 30       # Initial number of descriptors
MUTATION = 0.01   # Mutation rate
PROBABILITY = TRUE # Use probability for AUC-based fitness
testweight = 1.05 # Weight for test set performance

## Process control
LOAD = FALSE
OPT = TRUE
EVAL = TRUE

## Output options
APPLICABILITY = TRUE
ADDITIONAL.MODELS = TRUE
POP.EVAL = TRUE

## Source the main GA-SVM script
source("tx5b00465_si_002/GA-SVM_v2.5_main_script.r")
')

config_file <- paste0("run_", OUTPUT_PREFIX, ".R")
writeLines(config_script, config_file)
Sys.chmod(config_file, "755")
cat(paste("Created GA-SVM configuration script:", config_file, "\n"))

# Summary
cat("\n=======================================================\n")
cat("Data preparation complete!\n")
cat("\nInput summary:\n")
cat(paste("- Input file:", INPUT_CSV_FILE, "\n"))
cat(paste("- Total compounds:", nrow(compounds), "\n"))
cat(paste("- Successfully processed:", nrow(compounds_clean), "\n"))
cat(paste("- Failed processing:", nrow(compounds) - nrow(compounds_clean), "\n"))

cat("\nOutput files created:\n")
cat(paste("- Descriptors:", desc_file, "\n"))
cat(paste("- Toxicity data:", tox_file, "\n"))
cat(paste("- Train/test/val split:", split_file, "\n"))
cat(paste("- GA-SVM config:", config_file, "\n"))

cat("\nTo run GA-SVM training:\n")
cat(paste("Rscript", config_file, "\n"))
cat("\n=======================================================\n")