##########################################################################
##########################################################################
#
#   This is a R script. Please change the file extension to .r to use the file properly.
#   On notes how to use this file please see the file GA-SVM_v2.5_main_script.txt (or .r) 
#
#   This script was written by Dr. Denis Mulliner and published as supporting material 
#   to a publication titled "Computational Models for Human and Animal Hepatotoxicity with a Global Application Scope"
#   Author: Mulliner, Denis; Schmidt, Friedemann; Stolte, Manuela; Spirkl, Hans-Peter; Czich, Andreas; Amberg, Alexander
#   published in the Journal of Chemical Research in Toxicology in 2016.
#   Please cite what using this script!
#
##########################################################################
##########################################################################


## GA-SVM_2.5_load_data.r
## by Denis Mulliner
##
## a part of the script loading the data

cat(" ****-- LOADING DATA --**** \n\n")
load.folder=""

## LOADNIG TOX DATA
for (tox.file in TOXICITY.FILES)
{
if (file.exists(tox.file)) {
                             cat("Loading toxicity data from:\n  ",tox.file,"\n")
                             toxdata <- read.delim(tox.file,header=TRUE,na.strings="NA")
                             colnames(toxdata)[1] <- "ID"
                    } else { cat("ERROR - Tox data file does not exist:\n",tox.file,"\n");q(save="no")}

tox.col <- grep(ENDPOINT, colnames(toxdata))[1]
if (exists("tox.col")) {if (length(tox.col)==0|is.na(tox.col)) {
                                                                 cat("ERROR - tox data column ",tox.col.search," not found!\n");q(save="no")
                                         } else {cat("Reading tox data from column ",colnames(toxdata)[tox.col],"\n") }
                } else {cat("ERROR - tox data column with ",ENDPOINT," not found!\n");q(save="no")}
if (!exists("all.toxdata")) {all.toxdata <- toxdata[,c(1,tox.col)]
                     } else {all.toxdata <- rbind(all.toxdata,toxdata[,c(1,tox.col)])}
}; toxdata <- all.toxdata; rm(all.toxdata)

## LOADNIG Train/Test/Validation DATA
for (ttv.file in TRAIN.TEST.SPLIT.FILES)
{
if (file.exists(ttv.file)) {
                             cat("Loading train/test/val split information from:\n  ",ttv.file,"\n")
                             ttvdata <- read.delim(ttv.file,header=TRUE,na.strings="NA")
                             colnames(ttvdata)[1] <- "ID"
                    } else { cat("Train/test/validation data file does not exist!!\n");q(save="no")}

ttv.col <- grep(SET.column.header,colnames(ttvdata))
if (exists("ttv.col")) {if (length(ttv.col)==0|is.na(ttv.col))
                                         {#cat("ERROR - train/test/val column with ",PHASE,".",tt," not found!\n",sep="")#;q(save="no")
                                           cat("ERROR - SET.column.header",SET.column.header,"not found!\n",sep="");q(save="no")
                                         } else {cat("Reading SET information from column ",colnames(ttvdata)[ttv.col],"\n") }
                } else {cat("ERROR - SET.column.header",SET.column.header,"not found!\n",sep="")}#;q(save="no")}
if (!exists("all.ttvdata")) {all.ttvdata <- ttvdata[,c(1,ttv.col)]
                     } else {colnames(ttvdata)[ttv.col] <- colnames(all.ttvdata)[2]
                             all.ttvdata <- rbind(all.ttvdata,ttvdata[,c(1,ttv.col)])
                            }
}; ttvdata <- all.ttvdata; rm(all.ttvdata)
if (length(which(ttvdata[,2]==""))>0) ttvdata[which(ttvdata[,2]==""),2] <- NA

## preparing dataset
if (exists("dataset")) {rm(dataset)}
if (exists("all.DATA")) {rm(all.DATA)}
dataset <- ttvdata
rm(ttvdata)
##--print(head(dataset[,c(1,ncol(dataset))]))
## LOADING Descriptor files
#
# the process goes more or less like this:
#
#   1. descriptors A are loaded for data set a
#   2. descriptors A are loaded for data set b
#      .. same for dataset c, d, ...
#   3. loaded descriptors are connected by rbind into a data.frame 1
#   4. descriptors B are loaded for data set a
#   5. descriptors B are loaded for data set b
#   6. loaded descriptors are connected by rbind into a data.frame 2
#   7. dataframe 1 and 2 are joined with merge command
#      .. same for descriptors C, D, E ...
#

#for (DesC in unlist(strsplit(desc,"[.]")))
for (DesC in DESCRIPTOR.postfixes)
{
 for (dataset.prefix in DESCRIPTOR.FILES)
 {
#  if (DesC != "cur")
#  {
      descfile <- paste(load.folder,dataset.prefix,".",DesC,sep="")
      if (file.exists(descfile))
      {
        cat("loading",descfile,"\n")
        DATA <- read.delim(file=descfile,header=TRUE,sep=" ")
        #cat("ncol(DATA) =",ncol(DATA),"\n")
        colnames(DATA)[1] <- "ID"
        if (!exists("all.DATA")) {all.DATA <- DATA
                                  #cat("ncol(all.DATA) =",ncol(all.DATA),"\n")
                          } else {all.DATA <- rbind.fill(all.DATA,DATA)
                                  #cat("ncol(all.DATA) =",ncol(all.DATA),"\n")
                                  }
      } else { cat(paste("ERROR - descriptor file ",descfile," does not exist!\n"))}#;q(save="no")}
  #}
} ## END OF: for (dataset.prefix in DATASETS)
    if (exists("all.DATA")) {
                             dataset <- merge(dataset,all.DATA,by="ID",all=TRUE)
                             rm(all.DATA)
                            }
} # END OF: for (DesC in DESCRIPTOR.postfixes)

dataset <- merge(dataset,toxdata,by="ID",all="TRUE")

#tty <- gsub("RTECS[^.]*[.]","",tt)
#tty <- gsub("HTdl[1-9]00","",tty)
#tty <- gsub("HTl*o*g*MDZ*","",tty)

cat("\n\nwriting all files to",getwd(),"\n")


if (FALSE)
{
# load extra TTV files
for (EX in extraTTV)
{
 if (!is.na(EX)) {
                  cat("\nextra train/test/val split data\nLoading extra TTV data to replace the ttv data of some compounds\n")
                  cat("  loading:",EX,"\n")
                  if (file.exists(EX)) {
                                        ex.ttv <- read.delim(EX,head=TRUE,sep="\t")
                                        ex.ttv[,1] <- as.character(ex.ttv[,1])
                                        cat("head of loaded file:\n")
                                        print(head(ex.ttv))
                                        cat("\nidentifying compounds by ID in column:",colnames(ex.ttv)[1],"\n")
                                        cat("using information from column:",colnames(ex.ttv)[2],"\n")
                                        comps.in <- ex.ttv[which(ex.ttv[,1] %in% dataset[,1]),1]
                                        cat("adding ttv data\n")
                                        for (i in comps.in) {dataset[which(dataset[,1]==i),2] <- ex.ttv[which(ex.ttv[,1]==i),2]}
                                } else {cat(" ERROR - ",EX," not found\nextra TTV info can't be used!\n")}
                 }
}
}

# make class a factor!
if (length(grep("regression",svm.type))==0) dataset[,ncol(dataset)] <- as.factor(dataset[,ncol(dataset)])


cat("\n ****-- PREPROCESSING DATA --**** \n")
######## REMOVE DESCRIPTORS WITH NAs
#-----------print(dataset[which(!is.na(dataset[,2])),c(1,2)])
cat("---\n")
not.NA.in.ttv <- which(!is.na(dataset[,2]))
take.out <- rep(1,ncol(dataset))
if (!LOAD) {
     for (i in 3:(ncol(dataset)-1)) {
         if (NROW(unique(dataset[which(dataset[,2]=="train"),i]))==1) {take.out[i]=0
                                                                       #cat("Descriptor with only one unique value in training set:",colnames(dataset)[i],"\n")
                                                                      }
         NA.in.desc <- not.NA.in.ttv[which(is.na(dataset[not.NA.in.ttv,i]))]
         if (NROW(which(is.na(dataset[which(!is.na(dataset[,2])),i]))) > 0) {take.out[i]=0
                                                                             #cat("Descriptor with NAs:",colnames(dataset)[i],"\n")
                                                                             #cat("Compounds with NAs:",paste(dataset[NA.in.desc,1],collapse=" "),"\n")
                                                                            }
        }
cat("removing",NROW(which(take.out==0)),"descriptors with only one unique value or one or more NA in training, test, or validation set\n")
dataset <- dataset[,which(take.out==1)]
if (ncol(dataset) < 4) {cat("ERROR - all descriptors have been removed! Model building is not possible.\n        please make sure you do not have compounds without any descriptors in your training set!\n");q(save="no")}
}

if (NUMERIC.ENDPOINT) {cat(" The Endpoint to model is numeric and SVM does not tolerate\n")
                       cat(" NAs or missing values. So for all compounds with missing\n")
                       cat(" endpoint values the second column (train/test/val spliting)\n")
                       cat(" is set to NA!\n")
                       dataset[which(is.na(dataset[,ncol(dataset)])),2] <- NA}

if (DEBUG){ # check remaining columns for uniqueness and NAs
             cat("DEBUG: ----------------- summary of train data set:\n")
             cat("calling:   summary(dataset[which(dataset[,2]==\"train\"),])")
             cat(summary(dataset[which(dataset[,2]=="train"),]),"\n")
             cat("---------------------------------------------\n\n")
           }

### check if validation set is used
if (NROW(which(dataset[,2]=="val"))>0) {DoVal <- 1; cat("performing additional validation\n")}




########   ---> load population

cat("\nrunning",iter,"iterations with a population sitze of",pS,"\n")

if (!is.na(N.desc)) {
cat("creating start population of size",pS,"with a maximum of",N.desc,"positives\n")
## make small suggestions
first.guess <- matrix(rep(0,pS*(ncol(dataset)-3)),pS)
#head(first.guess)
for (i in 1:pS) {for (j in 1:N.desc) {first.guess[i,as.integer(runif(1)*(ncol(dataset)-3))] <- 1}}
}

### LOAD population to evaluate or to continue optimization
#LOAD = FALSE # TRUE #FALSE
if (LOAD) {
  for (l in 1:500) {
  if (l == 1) {LOAD.file.this <-  paste(FILENAME,"_runs_",l,"_population.txt",sep="")
               LOAD.file.next <-  paste(FILENAME,"_runs_",l+1,"_population.txt",sep="")
       } else {LOAD.file.this <- LOAD.file.next
               LOAD.file.next <-  paste(FILENAME,"_runs_",l+1,"_population.txt",sep="")
               }
 ### #  LOAD.file.2 <-  paste(func.prefix,"_popSize_",pS,"_iter_",iter,postfix,"_runs_",l+1,"_population.txt",sep="")
  if (file.exists(LOAD.file.this)&!file.exists(LOAD.file.next)) {
    cat(paste("loading data from run",l,"\n"))
    popul <- read.delim(LOAD.file.this,header=TRUE,sep=" ")
    cat(paste("Descriptors not used in the loaded population:",
                paste(colnames(dataset)[which(!colnames(dataset)[c(-1,-2,-ncol(dataset))] %in% colnames(popul))],collapse=", "),"\n"))
    dataset <- dataset[,c(1,2,which(colnames(dataset) %in% colnames(popul)),ncol(dataset))]
    cat("removing unused descriptors\n")
    bm.file <- paste(FILENAME,"_runs_",l,"_mean_best.txt",sep="")
    OUT.bm <- read.delim(file=bm.file,header=TRUE,sep=" ")
    colnames(popul) <- 0:(ncol(popul)-1)
    load.iter <- l
  break
  }
  }
}
if (exists("l")) {if (l==500) {LOAD=FALSE;print(paste("WARNING !!!! ---- nothing to load from file ",FILENAME));print("->starting new!")}}

