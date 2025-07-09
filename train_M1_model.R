#!/usr/bin/env Rscript

#
# Training Script for M1 Hepatotoxicity Model
# Using GA-SVM framework from Mulliner et al. (2016)
#
# This script prepares data and creates input files for the GA-SVM framework
#

cat("\n================================================\n")
cat("M1 Hepatotoxicity Model Training Setup\n")
cat("Using GA-SVM Framework from Mulliner et al. 2016\n")
cat("================================================\n\n")

# Step 1: Rename GA-SVM scripts from .txt to .r
cat("Step 1: Preparing GA-SVM scripts...\n")

ga_svm_files <- list.files("tx5b00465_si_002", pattern = "\\.txt$", full.names = TRUE)

for (file in ga_svm_files) {
  new_file <- gsub("\\.txt$", ".r", file)
  if (!file.exists(new_file)) {
    file.copy(file, new_file)
    cat(paste("Renamed:", basename(file), "->", basename(new_file), "\n"))
  }
}

# Step 2: Generate Sample Dataset
cat("\nStep 2: Generating sample hepatotoxicity dataset...\n")

# Create directories
dir.create("datasets", showWarnings = FALSE)

# Generate sample compounds with hepatotoxicity labels
# In practice, replace this with real data
set.seed(42)
n_compounds <- 500

# Sample SMILES strings (mix of known hepatotoxic and safe compounds)
sample_smiles <- c(
  # Known hepatotoxic
  "CC(=O)NC1=CC=C(C=C1)O",  # Acetaminophen
  "Cc1c(C)c2c(c(C)c1O)CCC(C)(COc1ccc(CC3SC(=O)NC3=O)cc1)O2",  # Troglitazone
  "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",  # Diclofenac
  "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",  # Indomethacin
  "Cc1oncc1C(=O)Nc1ccc(C(F)(F)F)cc1",  # Leflunomide
  # Known safe/low risk
  "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
  "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
  "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
  "CCO",  # Ethanol
  "CC(C)NCC(O)COc1ccc(OCC(O)CNC(C)C)cc1",  # Metoprolol
  # Additional drug-like molecules
  "Cc1ccccc1",  # Toluene
  "c1ccc(cc1)O",  # Phenol
  "CC(C)(C)NCC(O)c1ccc(O)c(O)c1",  # Salbutamol
  "CN(C)CCOC(c1ccccc1)c1ccccc1",  # Diphenhydramine
  "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
)

# Expand dataset with variations
all_smiles <- character(n_compounds)
all_labels <- integer(n_compounds)

# Fill first part with known compounds
n_known <- length(sample_smiles)
all_smiles[1:n_known] <- sample_smiles
all_labels[1:n_known] <- c(rep(1, 5), rep(0, 10))  # First 5 are toxic, rest are safe

# Generate additional synthetic SMILES for demo
# In real application, use actual compound database
for (i in (n_known + 1):n_compounds) {
  # Simple SMILES generation for demo
  n_carbons <- sample(1:20, 1)
  if (runif(1) < 0.3) {
    # Aromatic compound
    all_smiles[i] <- paste0("c1ccc", paste(rep("c", min(n_carbons, 6)), collapse = ""), "1")
  } else {
    # Aliphatic compound
    all_smiles[i] <- paste(rep("C", n_carbons), collapse = "")
  }
  # Random hepatotoxicity label (30% toxic)
  all_labels[i] <- ifelse(runif(1) < 0.3, 1, 0)
}

# Create compound IDs
compound_ids <- paste0("CMPD_", sprintf("%04d", 1:n_compounds))

# Step 3: Calculate Descriptors
cat("\nStep 3: Calculating molecular descriptors...\n")

library(reticulate)
source_python("hepatotox_descriptors.py")

# Initialize calculator
calculator <- HepatotoxicityDescriptors()

# Calculate descriptors in batches
batch_size <- 50
all_descriptors <- NULL

for (i in seq(1, n_compounds, by = batch_size)) {
  end_idx <- min(i + batch_size - 1, n_compounds)
  batch_smiles <- all_smiles[i:end_idx]
  
  cat(paste("Processing compounds", i, "to", end_idx, "...\n"))
  
  # Calculate descriptors
  batch_descriptors <- calculator$calculate_all_descriptors(batch_smiles)
  
  if (is.null(all_descriptors)) {
    all_descriptors <- batch_descriptors
  } else {
    all_descriptors <- rbind(all_descriptors, batch_descriptors)
  }
}

# Add compound IDs
all_descriptors <- cbind(ID = compound_ids, all_descriptors[, -1])  # Remove SMILES column

cat(paste("Calculated", ncol(all_descriptors) - 1, "descriptors for", nrow(all_descriptors), "compounds\n"))

# Step 4: Prepare GA-SVM Input Files
cat("\nStep 4: Preparing input files for GA-SVM...\n")

# 4.1: Save descriptor files by category (GA-SVM expects space-separated files)
# Categorize descriptors
desc_names <- names(all_descriptors)[-1]  # Exclude ID column

pharmacophore_desc <- grep("^(PP|PD|PA|PL|ND|NA|NL|DA|DL|AA|AL|SH)", desc_names, value = TRUE)
maccs_desc <- grep("^MKEY", desc_names, value = TRUE)
moe_desc <- grep("^(A_|B_|BCUT|CHIRAL|PEOE|PETITJEAN|RSYNTH|SLOGP|SMR|VADJMA|VDISTMA)", desc_names, value = TRUE)
custom_desc <- setdiff(desc_names, c(pharmacophore_desc, maccs_desc, moe_desc))

# Save descriptor files
if (length(pharmacophore_desc) > 0) {
  write.table(all_descriptors[, c("ID", pharmacophore_desc)], 
              "datasets/M1_descriptors.pharmacophore",
              sep = " ", row.names = FALSE, quote = FALSE)
}

if (length(maccs_desc) > 0) {
  write.table(all_descriptors[, c("ID", maccs_desc)], 
              "datasets/M1_descriptors.maccs",
              sep = " ", row.names = FALSE, quote = FALSE)
}

if (length(moe_desc) > 0) {
  write.table(all_descriptors[, c("ID", moe_desc)], 
              "datasets/M1_descriptors.moe",
              sep = " ", row.names = FALSE, quote = FALSE)
}

if (length(custom_desc) > 0) {
  write.table(all_descriptors[, c("ID", custom_desc)], 
              "datasets/M1_descriptors.custom",
              sep = " ", row.names = FALSE, quote = FALSE)
}

# 4.2: Save toxicity file (tab-separated as expected by GA-SVM)
tox_data <- data.frame(
  ID = compound_ids,
  Hepatotoxic_Human = all_labels
)
write.table(tox_data, "datasets/M1_toxicity.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# 4.3: Create train/test/validation split
# GA-SVM expects: "train", "test", "val" labels
set.seed(123)
n <- n_compounds
train_size <- round(0.6 * n)
test_size <- round(0.2 * n)
val_size <- n - train_size - test_size

# Stratified sampling to maintain class balance
toxic_idx <- which(all_labels == 1)
safe_idx <- which(all_labels == 0)

# Sample proportionally from each class
train_toxic <- sample(toxic_idx, round(train_size * length(toxic_idx) / n))
train_safe <- sample(safe_idx, round(train_size * length(safe_idx) / n))
train_idx <- c(train_toxic, train_safe)

remaining_toxic <- setdiff(toxic_idx, train_toxic)
remaining_safe <- setdiff(safe_idx, train_safe)

test_toxic <- sample(remaining_toxic, round(test_size * length(toxic_idx) / n))
test_safe <- sample(remaining_safe, round(test_size * length(safe_idx) / n))
test_idx <- c(test_toxic, test_safe)

val_idx <- setdiff(1:n, c(train_idx, test_idx))

# Create split assignment
split_assignment <- character(n)
split_assignment[train_idx] <- "train"
split_assignment[test_idx] <- "test"
split_assignment[val_idx] <- "val"

split_data <- data.frame(
  ID = compound_ids,
  SET_M1 = split_assignment
)

write.table(split_data, "datasets/M1_split.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nDataset split:\n")
cat(paste("- Training:", sum(split_assignment == "train"), "compounds\n"))
cat(paste("- Test:", sum(split_assignment == "test"), "compounds\n"))
cat(paste("- Validation:", sum(split_assignment == "val"), "compounds\n"))

# Step 5: Create GA-SVM Configuration Script
cat("\nStep 5: Creating GA-SVM configuration script...\n")

config_script <- '#!/usr/bin/env Rscript
#
# GA-SVM Configuration for M1 Hepatotoxicity Model
# Based on parameters from Mulliner et al. (2016)
#

cat("\\n=== M1 Hepatotoxicity Model Training ===\\n")
cat("Using GA-SVM v2.5\\n\\n")

## General settings
DEBUG = FALSE
FILENAME = "M1_hepatotoxicity_model.r"

## Descriptor files
DESCRIPTOR.FILES = "datasets/M1_descriptors"
DESCRIPTOR.postfixes = c("pharmacophore", "maccs", "moe", "custom")

## Toxicity data
TOXICITY.FILES = "datasets/M1_toxicity.txt"
ENDPOINT = "Hepatotoxic_Human"
NUMERIC.ENDPOINT = FALSE  # Binary classification

## Train/Test/Validation split
TRAIN.TEST.SPLIT.FILES = "datasets/M1_split.txt"
SET.column.header = "SET_M1"

## SVM Settings (from Mulliner et al. 2016)
svm.type = "nu-classification"
kernel.type = "radial"
tolerance = 0.01
epsilon = 0.1
NU = 0.7
GAMMA = 0.01

## Genetic Algorithm settings
pS = 100          # Population size (use 1000 for final model)
iter = 20         # Iterations per run
Nruns = 10        # Number of runs
N.desc = 30       # Initial number of descriptors
MUTATION = 0.01   # Mutation rate
PROBABILITY = TRUE # Use probability for AUC-based fitness
testweight = 1.05 # Weight for test set performance
penalty = NA      # No penalty for number of descriptors

## Cross-validation options
CROSSVAL = FALSE  # Use train/test split instead
Nfold = 10
LOOCV = FALSE

## Process control
LOAD = FALSE      # Start fresh
OPT = TRUE        # Run optimization
EVAL = TRUE       # Evaluate model after optimization

## Output options
APPLICABILITY = TRUE     # Calculate applicability domain
ADDITIONAL.MODELS = TRUE # Create additional models with all data
POP.EVAL = TRUE         # Evaluate population diversity

## Source the main GA-SVM script
source("tx5b00465_si_002/GA-SVM_v2.5_main_script.r")
'

writeLines(config_script, "run_M1_training.R")
Sys.chmod("run_M1_training.R", "755")

# Step 6: Create Quick Start Script
cat("\nStep 6: Creating quick start script...\n")

quick_start <- '#!/usr/bin/env Rscript
#
# Quick M1 Model for Testing
# Creates a simple model for immediate use while GA-SVM runs
#

library(e1071)
library(pROC)

# Load prepared data
descriptors <- read.table("datasets/M1_descriptors.moe", header = TRUE, sep = " ")
toxicity <- read.table("datasets/M1_toxicity.txt", header = TRUE, sep = "\\t")
splits <- read.table("datasets/M1_split.txt", header = TRUE, sep = "\\t")

# Merge data
data <- merge(descriptors, toxicity, by = "ID")
data <- merge(data, splits, by = "ID")

# Prepare training data
train_data <- data[data$SET_M1 == "train", ]
train_x <- as.matrix(train_data[, 2:(ncol(descriptors))])
train_y <- factor(train_data$Hepatotoxic_Human)

# Handle missing values
train_x[is.na(train_x)] <- 0

# Remove zero-variance columns
var_cols <- apply(train_x, 2, var) > 0
train_x <- train_x[, var_cols]

# Train simple SVM
cat("Training quick SVM model...\\n")
quick_svm <- svm(
  x = train_x,
  y = train_y,
  type = "nu-classification",
  kernel = "radial",
  nu = 0.7,
  gamma = 0.01,
  probability = TRUE
)

# Create model object for Shiny app
model <- list(
  svm_model = quick_svm,
  DESCRIPTORS = colnames(train_x),
  type = "M1_quick",
  creation_date = Sys.Date()
)

# Save model
save(model, file = "M1_hepatotoxicity_model.RData")
cat("Quick model saved as M1_hepatotoxicity_model.RData\\n")
'

writeLines(quick_start, "create_quick_model.R")
Sys.chmod("create_quick_model.R", "755")

# Final instructions
cat("\n================================================\n")
cat("Setup Complete!\n")
cat("================================================\n\n")

cat("To train the M1 model:\n\n")

cat("1. For QUICK MODEL (immediate use, ~1 minute):\n")
cat("   Rscript create_quick_model.R\n\n")

cat("2. For FULL GA-SVM MODEL (optimized, several hours):\n")
cat("   Rscript run_M1_training.R\n\n")

cat("Notes:\n")
cat("- The quick model uses all MOE descriptors without GA selection\n")
cat("- The full GA-SVM will select optimal descriptor subset\n")
cat("- For production, use pS=1000 in run_M1_training.R\n")
cat("- Full optimization may take 24+ hours\n\n")

cat("Required R packages:\n")
cat("- genalg, e1071, pROC, ggplot2, reshape, plyr, grid\n")
cat("- Install with: install.packages(c('genalg', 'e1071', 'pROC', 'ggplot2', 'reshape', 'plyr'))\n\n")

cat("The trained model will be saved as M1_hepatotoxicity_model.RData\n")
cat("and can be used directly with hepatotox_app.R\n")