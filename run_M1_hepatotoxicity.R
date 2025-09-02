#!/usr/bin/env Rscript
#
# GA-SVM Configuration for M1_hepatotoxicity
# Generated from CSV input with filtering options: QSAR=TRUE, STAT=FALSE
#

cat("\n=== M1_hepatotoxicity Model Training ===\n")
cat("Using Enhanced Descriptor Calculator + Filtering\n")
cat("Filtering applied: QSAR=TRUE, Statistical=FALSE\n\n")

## General settings
DEBUG = FALSE
FILENAME = "M1_hepatotoxicity_model.r"

## Descriptor files
DESCRIPTOR.FILES = "datasets/M1_hepatotoxicity_descriptors"
DESCRIPTOR.postfixes = "desc"

## Toxicity data
TOXICITY.FILES = "datasets/M1_hepatotoxicity_toxicity.txt"
ENDPOINT = "Hepatotoxic_Human"
NUMERIC.ENDPOINT = FALSE

## Train/Test/Validation split
TRAIN.TEST.SPLIT.FILES = "datasets/M1_hepatotoxicity_split.txt"
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

