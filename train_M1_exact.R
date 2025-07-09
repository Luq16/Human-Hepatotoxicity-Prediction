#!/usr/bin/env Rscript

#
# M1 Hepatotoxicity Model Training - Using Exact GA-SVM Scripts
# Following the exact methodology from Mulliner et al. (2016)
#

cat("\n=======================================================\n")
cat("M1 Hepatotoxicity Model Training\n")
cat("Using EXACT GA-SVM Scripts from Mulliner et al. (2016)\n")
cat("=======================================================\n\n")

# Step 1: Prepare the exact GA-SVM scripts
cat("Step 1: Preparing GA-SVM scripts...\n")

# First, rename all .txt files to .r as required by the original scripts
ga_svm_files <- list.files("tx5b00465_si_002", pattern = "\\.txt$", full.names = TRUE)
for (file in ga_svm_files) {
  new_file <- gsub("\\.txt$", ".r", file)
  if (!file.exists(new_file)) {
    file.copy(file, new_file)
    cat(paste("Renamed:", basename(file), "->", basename(new_file), "\n"))
  }
}

# Step 2: Generate sample data using the EXACT format from the original scripts
cat("\nStep 2: Generating sample data in exact GA-SVM format...\n")

# Create datasets directory
dir.create("datasets", showWarnings = FALSE)

# Use the EXACT fake data generation code from the original GA-SVM script
# This is copied directly from lines 219-256 in GA-SVM_v2.5_main_script.txt

## Make fake input for GA-SVM script (EXACT COPY from original)
A.size = 200 # size of dataset A
B.size = 50  # size of dataset B

## make fake tox data
A.tox <- data.frame("xIDx"=paste("X",1:A.size,sep=""),"X_class"=round(runif(A.size),digits=3),"Y_class"=as.integer(round(runif(A.size),digits=0)))
write.table(A.tox,file="datasets/test_tox_A.txt",row.names=F,sep="\t")
B.tox <- data.frame("xIDx"=paste("Y",1:B.size,sep=""),"X_class"=round(runif(B.size),digits=3),"Y_class"=as.integer(round(runif(B.size),digits=0)))
write.table(B.tox,file="datasets/test_tox_B.txt",row.names=F,sep="\t")

# make fake descriptors
A.xx <- data.frame("xIDx"=paste("X",1:A.size,sep=""),"X1"=runif(A.size))
for (i in 2:10) {A.xx[,ncol(A.xx)+1] <- runif(A.size); colnames(A.xx)[ncol(A.xx)] <- paste("X",i,sep="")}
write.table(A.xx,file="datasets/test_A.xx",row.names=F,sep=" ")
A.xx <- data.frame("xIDx"=paste("X",1:A.size,sep=""),"Y1"=runif(A.size))
for (i in 2:10) {A.xx[,ncol(A.xx)+1] <- runif(A.size); colnames(A.xx)[ncol(A.xx)] <- paste("Y",i,sep="")}
write.table(A.xx,file="datasets/test_A.yy",row.names=F,sep=" ")
A.xx <- data.frame("xIDx"=paste("X",1:A.size,sep=""),"Z1"=runif(A.size))
for (i in 2:10) {A.xx[,ncol(A.xx)+1] <- runif(A.size); colnames(A.xx)[ncol(A.xx)] <- paste("Z",i,sep="")}
write.table(A.xx,file="datasets/test_A.zz",row.names=F,sep=" ")

B.xx <- data.frame("xIDx"=paste("Y",1:B.size,sep=""),"X1"=runif(B.size))
for (i in 2:10) {B.xx[,ncol(B.xx)+1] <- runif(B.size); colnames(B.xx)[ncol(B.xx)] <- paste("X",i,sep="")}
write.table(B.xx,file="datasets/test_B.xx",row.names=F,sep=" ")
B.xx <- data.frame("xIDx"=paste("Y",1:B.size,sep=""),"Y1"=runif(B.size))
for (i in 2:10) {B.xx[,ncol(B.xx)+1] <- runif(B.size); colnames(B.xx)[ncol(B.xx)] <- paste("Y",i,sep="")}
write.table(B.xx,file="datasets/test_B.yy",row.names=F,sep=" ")
B.xx <- data.frame("xIDx"=paste("Y",1:B.size,sep=""),"Z1"=runif(B.size))
for (i in 2:10) {B.xx[,ncol(B.xx)+1] <- runif(B.size); colnames(B.xx)[ncol(B.xx)] <- paste("Z",i,sep="")}
write.table(B.xx,file="datasets/test_B.zz",row.names=F,sep=" ")

# make fake test training set split data
B.set <- data.frame("xIDx"=paste("Y",1:B.size,sep=""),"Y_SET_1"=c(rep("train",B.size-20),rep("test",10),rep("val",10))
                                                     ,"Y_SET_2"=c(rep("train",B.size-10),rep("test",5),rep("val",5)))
write.table(B.set,file="datasets/test_B_SET.txt",row.names=F,sep="\t")
A.set <- data.frame("xIDx"=paste("X",1:A.size,sep=""),"Y_SET_1"=c(rep("train",A.size-100),rep("test",50),rep("val",50))
                                                     ,"Y_SET_2"=c(rep("train",A.size-50),rep("test",25),rep("val",25)))
write.table(A.set,file="datasets/test_A_SET.txt",row.names=F,sep="\t")

cat("Sample data generated using exact GA-SVM format\n")

# Step 3: Create the EXACT input script as described in the original documentation
cat("\nStep 3: Creating exact GA-SVM input script...\n")

# This is the EXACT input file format from lines 69-213 in GA-SVM_v2.5_main_script.txt
input_script <- '
## input file for GA-SVM_V2.5
## by Denis Mulliner

## general settings
DEBUG=F # default DEBUG=F; DEBUG=T for a lot more text info on what is done to check for errors

# FILENAME: This is the master file name. It will be used to construct
# a base name (cutting away the .r suffix) which will be the base for all
# output file names and for the loaded file if continuing a run!
FILENAME="M1_hepatotoxicity_model.r"

##
## Descriptor files must have the following format:
## space seperated text file with header for every column
## First column: unique compound ID, header should be
##               defined but name does not matter because it
##               will be set to ID
## ALL Other Columns: the descriptors with a uniq header each
DESCRIPTOR.FILES=c("datasets/test_A","datasets/test_B") ## here using 2 sets of compounds
DESCRIPTOR.postfixes=c("xx","yy","zz")

##
## The files containing the Toxicity data mut have the following format
## Tab seperated text file with a header for every column
## First column: unique compound ID, header should be
##               defined but name does not matter because it
##               will be set to ID anyway
## ALL Other Columns: numeric or binary endpoint data.
##                    if column is binary values must be 0 and 1.
TOXICITY.FILES=c("datasets/test_tox_A.txt","datasets/test_tox_B.txt") ## here using two sets of compounds
ENDPOINT="Y_class"
NUMERIC.ENDPOINT=F # use: NUMERIC.ENDPOINT=T if a a regression should be done
                   # use: NUMERIC.ENDPOINT=F (default) if a classification model should be built

##
## The file(s) containing the data set split information
## must have the following format:
## First column: unique compound ID, header should be
##               defined but name does not matter because it
##               will be set to ID anyway
## ALL Other Columns: one of the values:
##                    "train" (use compound for model training)
##                    "test"  (use compound for model testing during feature selection)
##                    "val"   (use compound for model validation after each run of the genetic algorithm)
TRAIN.TEST.SPLIT.FILES=c("datasets/test_A_SET.txt","datasets/test_B_SET.txt")  ## here using two set of compound
SET.column.header="Y_SET_1" # defines the header of the column to use

##
## SVM Settings
##
svm.type="nu-classification"  ## other possible options:
kernel.type="radial"          ## other possible options: polynomial
                              ##                         linear
tolerance=0.01                ## a numeric value. for more info see manual of the
                              ## e1071 package containing the svm function
epsilon=0.1                   ## a numeric value. for more info see manual of the
                              ## e1071 package containing the svm function
NU=0.7                        ## a numeric value. for more info see manual of the
                              ## e1071 package containing the svm function
GAMMA=0.01                    ## a numeric value. for more info see manual of the
                              ## e1071 package containing the svm function

##
## genetic algorithm settings
pS = 20           ## population size (minimum 5); I used 1000 in most runs
iter = 10         ## number of iterations befor output is generated
Nruns = 5         ## number of runs repiting the specified number of iterations (iter) each time
N.desc = 10       ## number of descriptors choosen in first run
MUTATION=0.01     ## mutation rate used in the genetic algorithm
PROBABILITY=T     ## only for classification models
                  ## default: PROBABILITY=F  using the SVM model classification output
                  ##                         and fitness will be measured with ACC (Accuracy)
                  ##                         of classification with AUC=#N(correct classifications)/#N(total number of classifications)
                  ##       if PROBABILITY=T  use probability value
                  ##                         calculated by the SVM model to evaluate
                  ##                         predictive performance. So for fitness
                  ##                         the AUC (area under ROC curve)
                  ##                         will be used!
testweight=1.05   ## defining the weight of train and test set performance for
                  ## calculating the fitness
                  ## fitness = a*fitness(TEST) + (2-a)*fitness(train)
                  ## with a=testweight!
penalty=NA        ## default=NA; this can be used to penelize a choice of too many
                  ## descriptors and force the optimization towards less descriptors
                  ## penalty should be a numeric vector of length 2 e.g. c(40,0.5)
                  ## the penalty is introduced by reducing the fitness
                  ## if #N(number of descriptors > penalty[1]
                  ## by [ (#N(number of descriptors)-penalty[1])*penalty[2] ]^2
                  ## WARNING: this option was not well tested and not used regularly!



## cross-validation options
CROSSVAL=F ## default; if CROSSVAL=T fitness of an individual will not be calculated from
           ##          training and test set but using train+test compounds to build
           ##          the model and a cross-validation to test the model
Nfold=10   ## da a N-fold cross-validation
LOOCV=F    ## default; LOOCV=T: only effective when CROSSVAL=T; will do a leave one out
           ##          cross-validation. This is very time consuming and should not
           ##          be the choice for large optimizations

## PROCESS:
LOAD=F    ## default: LOAD=F; LOAD=T will try to load last run of optimization
          ##                  based on the FILENAME base name
OPT=T     ## default: OPT=F; OPT=T do an optimization run
EVAL=T    ## default: EVAL=F; EVAL=T evaluate the model output.
          ## if EVAL=T and OPT=F a model must be loaded that than can be evaluated

## using a default set of descriptors
#use.these.descriptors="test_descriptor_set.txt" ## this option runs a pseudo optimization
                                                ## using only the descriptor subset
                                                ## specified in the given text file
                                                ## with one descriptor name (header in descriptor files)
                                                ## per line and no headers!
                                                ## This will automatically set
                                                ## POP.EVAL=F and LOAD=F, they are not used!

APPLICABILITY=F  ## default; if APPLICABILITY=T an applicability domain
                 ##          analysis is run based on an euclidean distance calculation
                 ##          with the selected descriptors.
                 ##          so far this has not lead to any improvement
                 ##          so it is not default!

ADDITIONAL.MODELS=F ## default is FALSE; if ADDITIOnAL.MODELS=T than produce a lot of additional models
                    ##                   using all train, test and validation compounds for model training
                    ##                   and a pruning of the chosen descriptor set
POP.EVAL=T  ## default: POP.EVAL=F; if POP.EVAL=T the population is evaluated at
            ##                      the end of each run checking individuals for
            ##                      overlap and giving a feel for the optimization status
            ##                      of the model

## now loading the GA-SVM script and start running
source("tx5b00465_si_002/GA-SVM_v2.5_main_script.r")
q(save="no")
'

writeLines(input_script, "run_M1_exact.R")
cat("Created exact GA-SVM input script: run_M1_exact.R\n")

# Step 4: Create version with real descriptors for actual hepatotoxicity prediction
cat("\nStep 4: Creating version with real molecular descriptors...\n")

# Generate sample compound data with real SMILES and calculated descriptors
library(reticulate)
source_python("hepatotox_descriptors.py")

# Sample compounds with known hepatotoxicity
compounds <- data.frame(
  ID = c("ACETA", "TROG", "DICLO", "INDO", "IBUPR", "CAFF", "ASPIR", "ETOH"),
  SMILES = c(
    "CC(=O)NC1=CC=C(C=C1)O",  # Acetaminophen (hepatotoxic)
    "Cc1c(C)c2c(c(C)c1O)CCC(C)(COc1ccc(CC3SC(=O)NC3=O)cc1)O2",  # Troglitazone (hepatotoxic)
    "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",  # Diclofenac (hepatotoxic)
    "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",  # Indomethacin (hepatotoxic)
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen (safe)
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine (safe)
    "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin (safe)
    "CCO"  # Ethanol (safe)
  ),
  Hepatotoxic = c(1, 1, 1, 1, 0, 0, 0, 0)
)

# Calculate descriptors
calculator <- HepatotoxicityDescriptors()
descriptors <- calculator$calculate_all_descriptors(compounds$SMILES)

# Save in GA-SVM format
desc_data <- cbind(xIDx = compounds$ID, descriptors[, -1])  # Remove SMILES column
write.table(desc_data, "datasets/hepatotoxicity_descriptors.desc", sep = " ", row.names = FALSE, quote = FALSE)

# Save toxicity data
tox_data <- data.frame(xIDx = compounds$ID, Hepatotoxic = compounds$Hepatotoxic)
write.table(tox_data, "datasets/hepatotoxicity_toxicity.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Create train/test split
set.seed(42)
train_ids <- sample(compounds$ID, 5)
test_ids <- setdiff(compounds$ID, train_ids)
val_ids <- sample(test_ids, 1)
test_ids <- setdiff(test_ids, val_ids)

split_data <- data.frame(
  xIDx = compounds$ID,
  SET_hepatotoxicity = ifelse(compounds$ID %in% train_ids, "train",
                             ifelse(compounds$ID %in% test_ids, "test", "val"))
)
write.table(split_data, "datasets/hepatotoxicity_split.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Create input script for real data
real_input_script <- gsub("datasets/test_A","datasets/test_B", "datasets/hepatotoxicity_descriptors", input_script)
real_input_script <- gsub("datasets/test_tox_A.txt\",\"datasets/test_tox_B.txt", "datasets/hepatotoxicity_toxicity.txt", real_input_script)
real_input_script <- gsub("datasets/test_A_SET.txt\",\"datasets/test_B_SET.txt", "datasets/hepatotoxicity_split.txt", real_input_script)
real_input_script <- gsub("DESCRIPTOR.FILES=c\\(\"datasets/test_A\",\"datasets/test_B\"\\)", 
                         "DESCRIPTOR.FILES=\"datasets/hepatotoxicity_descriptors\"", real_input_script)
real_input_script <- gsub("DESCRIPTOR.postfixes=c\\(\"xx\",\"yy\",\"zz\"\\)", 
                         "DESCRIPTOR.postfixes=\"desc\"", real_input_script)
real_input_script <- gsub("TOXICITY.FILES=c\\(\"datasets/test_tox_A.txt\",\"datasets/test_tox_B.txt\"\\)", 
                         "TOXICITY.FILES=\"datasets/hepatotoxicity_toxicity.txt\"", real_input_script)
real_input_script <- gsub("TRAIN.TEST.SPLIT.FILES=c\\(\"datasets/test_A_SET.txt\",\"datasets/test_B_SET.txt\"\\)", 
                         "TRAIN.TEST.SPLIT.FILES=\"datasets/hepatotoxicity_split.txt\"", real_input_script)
real_input_script <- gsub("Y_class", "Hepatotoxic", real_input_script)
real_input_script <- gsub("Y_SET_1", "SET_hepatotoxicity", real_input_script)

writeLines(real_input_script, "run_M1_real_data.R")

cat("Created real data input script: run_M1_real_data.R\n")

# Final instructions
cat("\n=======================================================\n")
cat("Setup Complete - Using EXACT GA-SVM Scripts!\n")
cat("=======================================================\n\n")

cat("To run the M1 model training:\n\n")

cat("1. Test with sample data (recommended first):\n")
cat("   Rscript run_M1_exact.R\n\n")

cat("2. Run with real molecular descriptors:\n")
cat("   Rscript run_M1_real_data.R\n\n")

cat("Required R packages:\n")
cat("install.packages(c('genalg', 'e1071', 'pROC', 'ggplot2', 'reshape', 'plyr', 'grid'))\n\n")

cat("Notes:\n")
cat("- These scripts use the EXACT GA-SVM code from Mulliner et al. (2016)\n")
cat("- No modifications made to the original algorithm\n")
cat("- Output will be M1_hepatotoxicity_model.RData\n")
cat("- For production runs, increase pS to 1000 in the input script\n")
cat("- The scripts will create all output files as described in the original paper\n")