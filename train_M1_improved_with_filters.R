#!/usr/bin/env Rscript

#
# M1 Hepatotoxicity Model Training - With Optional Feature Selection
# Enhanced version with configurable filtering options
#

cat("\n=======================================================\n")
cat("M1 Hepatotoxicity Model Training with Optional Filters\n")
cat("Enhanced Descriptor Calculator + GA-SVM Framework\n")
cat("=======================================================\n\n")

# Configuration options
APPLY_QSAR_FILTERS <- TRUE      # Filter problematic molecules
APPLY_STAT_FILTERS <- FALSE     # Filter descriptors statistically
VERBOSE_OUTPUT <- TRUE          # Show detailed progress

cat("Configuration:\n")
cat(paste("- QSAR molecular filters:", ifelse(APPLY_QSAR_FILTERS, "ENABLED", "DISABLED"), "\n"))
cat(paste("- Statistical descriptor filters:", ifelse(APPLY_STAT_FILTERS, "ENABLED", "DISABLED"), "\n"))
cat(paste("- Verbose output:", ifelse(VERBOSE_OUTPUT, "ENABLED", "DISABLED"), "\n\n"))

# Step 1: Prepare GA-SVM scripts
cat("Step 1: Preparing GA-SVM scripts...\n")

ga_svm_files <- list.files("tx5b00465_si_002", pattern = "\\.txt$", full.names = TRUE)
for (file in ga_svm_files) {
  new_file <- gsub("\\.txt$", ".r", file)
  if (!file.exists(new_file)) {
    file.copy(file, new_file)
    if (VERBOSE_OUTPUT) cat(paste("Renamed:", basename(file), "->", basename(new_file), "\n"))
  }
}

# Step 2: Generate hepatotoxicity dataset
cat("\nStep 2: Generating hepatotoxicity dataset...\n")

library(reticulate)
source_python("hepatotox_descriptors_improved.py")

# Create directories
dir.create("datasets", showWarnings = FALSE)

# Known hepatotoxic and non-hepatotoxic compounds (same as before)
compounds <- data.frame(
  ID = c(
    # Known hepatotoxic compounds
    "ACETAMINOPHEN", "TROGLITAZONE", "DICLOFENAC", "INDOMETHACIN", 
    "KETOCONAZOLE", "VALPROIC_ACID", "PHENYTOIN", "HALOTHANE",
    "ISONIAZID", "RIFAMPIN", "CHLORPROMAZINE", "FELBAMATE",
    "NEFAZODONE", "PEMOLINE", "TROVAFLOXACIN", "BROMFENAC",
    "TOLCAPONE", "LABETALOL", "METHYLDOPA", "NITROFURANTOIN",
    
    # Known safe/low hepatotoxicity compounds  
    "IBUPROFEN", "CAFFEINE", "ASPIRIN", "METFORMIN",
    "LISINOPRIL", "ATENOLOL", "OMEPRAZOLE", "LORATADINE",
    "SIMVASTATIN", "CIMETIDINE", "RANITIDINE", "FUROSEMIDE",
    "HYDROCHLOROTHIAZIDE", "AMLODIPINE", "WARFARIN", "DIGOXIN",
    "PROPRANOLOL", "CAPTOPRIL", "NIFEDIPINE", "VERAPAMIL"
  ),
  SMILES = c(
    # Hepatotoxic compounds
    "CC(=O)NC1=CC=C(C=C1)O",
    "Cc1c(C)c2c(c(C)c1O)CCC(C)(COc1ccc(CC3SC(=O)NC3=O)cc1)O2",
    "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
    "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",
    "CC(C)(C#N)N1C=NC(=CN1)C2=CC=C(C=C2)OC3=CC=CC=C3Cl",
    "CCCC(C(=O)O)CCC",
    "C1=CC=C(C=C1)C2C(=O)NC(=O)N2C3=CC=CC=C3",
    "FC(F)CBr(Cl)C(F)(F)F",
    "C1=CN=CC(=C1)C(=O)NN",
    "CC1=CC(=C(C(=C1C)O)C(=O)NCCN(C)C)C",
    "CN(C)CCCN1C2=CC=CC=C2SC3=C1C=C(C=C3)Cl",
    "NC(=O)OCC1=CC=CC=C1C#N",
    "CCC1=CC2=C(C=C1)N=C(N2CC3=CC=C(C=C3)Cl)C4=CC=CC=N4",
    "CC1=CC=C(C=C1)C2(CCNC2=O)C3=CC=CC=C3",
    "COC1=C(C=C2C(=C1)N(C=C(C2=O)C(=O)O)C3CC3)N4CCN(CC4)C",
    "CC1=CC=CC=C1NC2=C(C=CC(=C2)CC(=O)O)Br",
    "CC(C)(C)C1=CC(=C(C=C1)O)C(=O)NC2=CC=C(C=C2)C(F)(F)F",
    "CC(C)NCC(COC1=CC=C(C=C1)CC(C)C)O",
    "CC(CC1=CC(=C(C=C1)O)O)C(=O)O",
    "CC1=CC=C(C=C1)C(=NNC(=O)N)C2=CC=CC=N2",
    
    # Safe compounds
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "CC(=O)OC1=CC=CC=C1C(=O)O",
    "CN(C)C(=N)NC(=N)N",
    "CCCCC1=CC=C(C=C1)C(CCN2CCCC2C(=O)O)O",
    "CC(C)NCC(COC1=CC=C(C=C1)CC(=O)N)O",
    "COC1=CC2=C(C=C1)C(=CN2)CS(=O)C3=NC4=CC=CC=C4N3",
    "CCC1=CC=CC=C1C(=O)N2CCC(CC2)C(C3=CC=CC=N3)C4=CC=C(C=C4)Cl",
    "CCC(C)(C)C(=O)OC1CC(C=C2C1=CC=C(C2)C)C=CC3CC(CC(C3=O)O)O",
    "C1=CN=C(N1)CNCC2=C(N=CN2)CSC",
    "CNC(=NCCSC1=NC2=CC=CC=C2N1)N[N+](=O)[O-]",
    "C1=CC=C(C=C1)NS(=O)(=O)C2=CC(=C(C=C2)Cl)N",
    "C1=CC=C(C=C1)S(=O)(=O)NC2=NC3=CC=CC=C3S2",
    "CCOC(=O)C1=C(NC(=C(C1C2=CC=CC=C2[N+](=O)[O-])C(=O)OC)C)C",
    "CC1=CC2=C(C=C1)C(=CN2)CC(=O)NC3=CC=C(C=C3)C(=O)O",
    "CCC1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)O)O",
    "CC(C)NCC(COC1=CC=CC2=C1C=CN2)O",
    "CC(CS)C(=O)N1CCCC1C(=O)O",
    "COC(=O)C1=C(NC(=C(C1C2=CC=CC=C2[N+](=O)[O-])C)C)C",
    "CC(C)NCC(COC1=CC=C(C=C1)CC(C)(C)C#N)O"
  ),
  Hepatotoxic = c(
    rep(1, 20),  # Hepatotoxic
    rep(0, 20)   # Non-hepatotoxic
  ),
  stringsAsFactors = FALSE
)

cat(paste("Created dataset with", nrow(compounds), "compounds\n"))
cat(paste("- Hepatotoxic:", sum(compounds$Hepatotoxic), "compounds\n"))
cat(paste("- Non-hepatotoxic:", sum(1 - compounds$Hepatotoxic), "compounds\n"))

# Step 3: Calculate descriptors with optional filtering
cat("\nStep 3: Calculating molecular descriptors with optional filtering...\n")

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

# Save descriptor file
write.table(all_descriptors, "datasets/hepatotoxicity_filtered_descriptors.desc", 
            sep = " ", row.names = FALSE, quote = FALSE)

# Save toxicity file
tox_data <- data.frame(
  ID = compounds$ID,
  Hepatotoxic_Human = compounds$Hepatotoxic
)
write.table(tox_data, "datasets/hepatotoxicity_toxicity.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Create stratified train/test/validation split (same as before)
set.seed(42)
toxic_idx <- which(compounds$Hepatotoxic == 1)
safe_idx <- which(compounds$Hepatotoxic == 0)

train_toxic <- sample(toxic_idx, 12)
train_safe <- sample(safe_idx, 12)
train_idx <- c(train_toxic, train_safe)

remaining_toxic <- setdiff(toxic_idx, train_toxic)
remaining_safe <- setdiff(safe_idx, train_safe)

test_toxic <- sample(remaining_toxic, 4)
test_safe <- sample(remaining_safe, 4)
test_idx <- c(test_toxic, test_safe)

val_idx <- setdiff(1:nrow(compounds), c(train_idx, test_idx))

split_assignment <- character(nrow(compounds))
split_assignment[train_idx] <- "train"
split_assignment[test_idx] <- "test"
split_assignment[val_idx] <- "val"

split_data <- data.frame(
  ID = compounds$ID,
  SET_hepatotoxicity = split_assignment
)

write.table(split_data, "datasets/hepatotoxicity_split.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Print split statistics
train_toxic_count <- sum(compounds$Hepatotoxic[train_idx])
test_toxic_count <- sum(compounds$Hepatotoxic[test_idx])
val_toxic_count <- sum(compounds$Hepatotoxic[val_idx])

cat("\nDataset split (stratified):\n")
cat(paste("- Training:", length(train_idx), "compounds (", train_toxic_count, "toxic,", 
          length(train_idx) - train_toxic_count, "safe)\n"))
cat(paste("- Test:", length(test_idx), "compounds (", test_toxic_count, "toxic,", 
          length(test_idx) - test_toxic_count, "safe)\n"))
cat(paste("- Validation:", length(val_idx), "compounds (", val_toxic_count, "toxic,", 
          length(val_idx) - val_toxic_count, "safe)\n"))

# Step 5: Create GA-SVM configuration script
cat("\nStep 5: Creating GA-SVM configuration...\n")

config_script <- paste0('#!/usr/bin/env Rscript
#
# GA-SVM Configuration for M1 Hepatotoxicity Model
# Generated with filtering options: QSAR=', APPLY_QSAR_FILTERS, ', STAT=', APPLY_STAT_FILTERS, '
#

cat("\\n=== M1 Hepatotoxicity Model Training ===\\n")
cat("Using Enhanced Descriptor Calculator + Filtering\\n")
cat("Filtering applied: QSAR=', APPLY_QSAR_FILTERS, ', Statistical=', APPLY_STAT_FILTERS, '\\n\\n")

## General settings
DEBUG = FALSE
FILENAME = "M1_hepatotoxicity_filtered_model.r"

## Descriptor files
DESCRIPTOR.FILES = "datasets/hepatotoxicity_filtered_descriptors"
DESCRIPTOR.postfixes = "desc"

## Toxicity data
TOXICITY.FILES = "datasets/hepatotoxicity_toxicity.txt"
ENDPOINT = "Hepatotoxic_Human"
NUMERIC.ENDPOINT = FALSE

## Train/Test/Validation split
TRAIN.TEST.SPLIT.FILES = "datasets/hepatotoxicity_split.txt"
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

writeLines(config_script, "run_M1_filtered.R")
Sys.chmod("run_M1_filtered.R", "755")

# Step 6: Create quick validation script
validation_script <- paste0('#!/usr/bin/env Rscript
#
# Quick validation with filtering options
# QSAR=', APPLY_QSAR_FILTERS, ', Statistical=', APPLY_STAT_FILTERS, '
#

library(e1071)
library(pROC)

# Load data
descriptors <- read.table("datasets/hepatotoxicity_filtered_descriptors.desc", header = TRUE, sep = " ")
toxicity <- read.table("datasets/hepatotoxicity_toxicity.txt", header = TRUE, sep = "\\t")
splits <- read.table("datasets/hepatotoxicity_split.txt", header = TRUE, sep = "\\t")

# Merge data
data <- merge(descriptors, toxicity, by = "ID")
data <- merge(data, splits, by = "ID")

# Prepare training data
train_data <- data[data$SET_hepatotoxicity == "train", ]
test_data <- data[data$SET_hepatotoxicity == "test", ]

train_x <- as.matrix(train_data[, 2:(ncol(descriptors))])
train_y <- factor(train_data$Hepatotoxic_Human)
test_x <- as.matrix(test_data[, 2:(ncol(descriptors))])
test_y <- factor(test_data$Hepatotoxic_Human)

# Handle missing values and zero variance
train_x[is.na(train_x)] <- 0
test_x[is.na(test_x)] <- 0

var_cols <- apply(train_x, 2, var) > 1e-6
train_x <- train_x[, var_cols]
test_x <- test_x[, var_cols]

cat("Training with filtering (QSAR=', APPLY_QSAR_FILTERS, ', Stat=', APPLY_STAT_FILTERS, '):\\n")
cat("- Training data:", nrow(train_x), "compounds,", ncol(train_x), "descriptors\\n")
cat("- Test data:", nrow(test_x), "compounds\\n")

# Train SVM
svm_model <- svm(
  x = train_x,
  y = train_y,
  type = "nu-classification",
  kernel = "radial",
  nu = 0.7,
  gamma = 0.01,
  probability = TRUE
)

# Evaluate
test_pred <- predict(svm_model, test_x, probability = TRUE)
test_prob <- attr(test_pred, "probabilities")[, "1"]

test_accuracy <- mean(as.numeric(as.character(test_pred)) == as.numeric(as.character(test_y)))
test_roc <- roc(as.numeric(as.character(test_y)), test_prob, quiet = TRUE)
test_auc <- as.numeric(test_roc$auc)

cat("\\nPerformance with filtering:\\n")
cat("- Accuracy:", round(test_accuracy, 3), "\\n")
cat("- AUC:", round(test_auc, 3), "\\n")

# Create model for Shiny app
model <- list(
  svm_model = svm_model,
  DESCRIPTORS = colnames(train_x),
  type = "M1_filtered",
  qsar_filters = ', APPLY_QSAR_FILTERS, ',
  stat_filters = ', APPLY_STAT_FILTERS, ',
  creation_date = Sys.Date(),
  test_accuracy = test_accuracy,
  test_auc = test_auc
)

save(model, file = "M1_hepatotoxicity_model.RData")
cat("\\nFiltered model saved as M1_hepatotoxicity_model.RData\\n")
')

writeLines(validation_script, "validate_filtered_model.R")
Sys.chmod("validate_filtered_model.R", "755")

# Final summary
cat("\n=======================================================\n")
cat("Enhanced M1 Model Training Setup Complete!\n")
cat("=======================================================\n\n")

cat("Configuration used:\n")
cat(paste("- QSAR molecular filters:", ifelse(APPLY_QSAR_FILTERS, "ENABLED", "DISABLED"), "\n"))
cat(paste("- Statistical descriptor filters:", ifelse(APPLY_STAT_FILTERS, "ENABLED", "DISABLED"), "\n"))

cat("\nTo run:\n")
cat("1. Quick validation: Rscript validate_filtered_model.R\n")
cat("2. Full GA-SVM: Rscript run_M1_filtered.R\n\n")

cat("Key benefits of this version:\n")
cat("- Configurable filtering options at the top of the script\n")
cat("- QSAR filters remove problematic molecules\n")
cat("- Statistical filters reduce descriptor redundancy\n")
cat("- Compatible with all existing GA-SVM functionality\n")
cat("- Easy to compare filtered vs unfiltered results\n")