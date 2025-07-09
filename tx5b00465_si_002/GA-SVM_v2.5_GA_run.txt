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



## --- do a single run without optimization?
if (!is.na(use.these.descriptors)) # using only these given descriptors
{
 if (NROW(use.these.descriptors)==1)
 {
  if (file.exists(use.these.descriptors))
  {
   cat("loading predefined names of descriptors to use from file\n   the file is expected to contain the name of a descriptor on every line\n")
   use.D <- read.delim(use.these.descriptors,header=FALSE)
   if (exists("use.D")) {cat("... loading names succesful!\n")}
   use.these.descriptors <- as.character(use.D[,1])
  }
 }
 cat("using only these descriptors:\n")
 print(use.these.descriptors)
 if (!all(use.these.descriptors %in% colnames(dataset)))
 {
  cat("ERROR - specified descriptors not among loaded descriptors!\n")
  cat("        missing:\n")
  print(use.these.descriptors[which(!use.these.descriptors %in% colnames(dataset))])
  q(save="no")
 }
 
 Chromosom <- rep(0,ncol(dataset)-3)
 Chromosom[which(colnames(dataset) %in% use.these.descriptors)-2] <- 1
# cat("\nThe Chromosom:\n")
# print(Chromosom)
# cat("\nThe Descriptors:\n")
# print(colnames(dataset)[which(Chromosom==1)+2])
# cat("evaluating performance\n")
 eval.these.descriptors <- evalFunc(Chromosom)
# cat(eval.these.descriptors,"\n")
 GAmodel <- list("type"="binary chromosom","size"=ncol(dataset)-3,"popSize"=1,"iters"=1,"suggestions"=NA,
                 "population"=matrix(Chromosom,nrow=1),
                 "elitism"=NA,
                 "mutationChance"=0,
                 "evaluations"=eval.these.descriptors,
                        "best"=eval.these.descriptors,
                        "mean"=eval.these.descriptors)

 OUT.bm <- data.frame("best"=eval.these.descriptors,"mean"=eval.these.descriptors)
 ii=1 # this is the number of the current itteration!

 plot.GA.output(GAmodel)

 bm.file <- paste(FILENAME,"_runs_",ii,"_mean_best.txt",sep="")
 write.table(OUT.bm,file=bm.file,row.names=FALSE)
} else {
#############################
## ready to otimize

# setting the iterations:
#     iter are the internal iterations. After these a plot is saved
#     Nruns ar the external iterations, i.e. the repetition of iter internal iterations
if (LOAD) {
   ii.min <- load.iter + 1
   ii.max <- load.iter + Nruns
} else {
   ii.min <- 1
   ii.max <- Nruns
}

if (DEBUG)
{
 cat("\nDEBUG:---------------\ndata set overview:\nwith data:\n")
 cat("N(train) = ",NROW(which(dataset[,2]=="train"&!is.na(dataset[,ncol(dataset)]))),"\n")
 cat("N(test)  =",NROW(which(dataset[,2]=="test"&!is.na(dataset[,ncol(dataset)]))),"\n")
 cat("N(val)   =",NROW(which(dataset[,2]=="val"&!is.na(dataset[,ncol(dataset)]))),"\n")
 cat("without data:\n")
 cat("N(train) = ",NROW(which(dataset[,2]=="train"&is.na(dataset[,ncol(dataset)]))),"\n")
 cat("N(test)  =",NROW(which(dataset[,2]=="test"&is.na(dataset[,ncol(dataset)]))),"\n")
 cat("N(val)   =",NROW(which(dataset[,2]=="val"&is.na(dataset[,ncol(dataset)]))),"\n")
}


cat("\n")
cat("\n")
cat("-------------\n")
cat("optimizing...\n")
cat("-------------\n")
cat(paste("Output will be written to:",FILENAME,"...\n"))
chrom.size = ncol(dataset)-3
if (is.na(MUTATION)&!is.na(increaseMUTATION)) {MUTATION=(1/(chrom.size +1))*increaseMUTATION} ## this is the default value of rbga
if (is.na(MUTATION)) {MUTATIONtext=paste(round(1/(chrom.size+1),digits=4))} else {MUTATIONtext=paste(round(MUTATION,digits=4))}

# print out used parameters before going into optimization phase
cat("\nUsed data:\n")
cat("   N(train) = ",NROW(which(dataset[,2]=="train"&!is.na(dataset[,ncol(dataset)]))),"\n")
cat("   N(test)  =",NROW(which(dataset[,2]=="test"&!is.na(dataset[,ncol(dataset)]))),"\n")
cat("   N(val)   =",NROW(which(dataset[,2]=="val"&!is.na(dataset[,ncol(dataset)]))),"\n")
cat("\nUsed parameters for the SVM process:\n")
cat("   kernel = ",kernel.type,"\n   SVM = ",svm.type,"\n")
if (!is.na(NU)) {cat("   nu = ",NU,"\n")} else {cat("   Using default nu value\n")}
if (!is.na(GAMMA)) {cat("   gamma = ",GAMMA,"\n")} else {cat("   Using default gamma value\n")}
if (!is.na(epsilon)) {cat("   epsilon = ",epsilon,"\n")} else {cat("   Using default epsilon value\n")}
cat("   tolerance = ",tolerance,"\n")
cat("\nUsed parameters for the genetic algorithm:\n   popSize = ",pS,"\n   mutation = ",MUTATION,"\n   Number of iterations before generating output = ",iter,"\n")
cat("   Number of descriptors to choose from = ",chrom.size,"\n")
cat("   Number of descriptors chosen for each individuum in first run = ")
if (is.na(N.desc)) {cat("default (approx. 10%)\n")} else {cat(N.desc,"\n")}
if (PROBABILITY) {FIT.txt <- "AUC"} else {FIT.txt <- "ACC"}
if (NUMERIC.ENDPOINT) {FIT.txt <- fitness}
cat("   Using weighting for test set fitness: ",testweight,"\n")
if (length(grep("regression",svm.type))==0) {
                                              if (all(class.weight==1)) {cat("   using same weight (1) for all classes\n")
                                                                 } else {cat("   using the following class weights:\n")
                                                                         cat(round(class.weight,digits=3),"\n")
                                                                        }
                                            }
cat("\n----->>>>> ... running\n")
for (ii in ii.min:ii.max){
cat("  Iterations ",iter*(ii-1)+1," to ",iter*ii,"\n")
if (ii == ii.min) {# TRUE if it is the first run of the GA-optimization
                   if (LOAD) {# TRUE when loading a population
                    if (DEBUG) {print("DEBUG: rbga.bin -- first run loading population")}
                    GAmodel <- rbga.bin(size = chrom.size,
                                        popSize = pS,
                                        iters = iter,
                                        elitism = T,
                                        evalFunc = evalFunc,
                                        mutationChance=MUTATION,
                                        suggestions=matrix(as.matrix(popul[order(popul[,1])[1:(nrow(popul)-1)],-1]),nrow(popul)-1))
                     } else {# TRUE if no population is loaded
                     if (DEBUG) {print("DEBUG: rbga.bin -- first run not loading a population")}
                     if (is.na(N.desc)) {# TRUE if nothign is loaded and a standard run is done
                     if (DEBUG) {print("DEBUG: rbga.bin -- first run: standard")}
                     GAmodel <- rbga.bin(size = chrom.size,
                                         popSize = pS,
                                         iters = iter,
                                         mutationChance=MUTATION,
                                         elitism = T,
                                         evalFunc = evalFunc)
                     } else {# TRUE if N.desc is specified, no population is loaded but a predefined population is used
                     if (DEBUG) {print("DEBUG: rbga.bin -- first run not loading a population with user def. pop (N.desc)")}
                     GAmodel <- rbga.bin(size = chrom.size,
                                         popSize = pS,
                                         iters = iter,
                                         elitism = T,
                                         evalFunc = evalFunc,
                                         mutationChance=MUTATION,
                                         suggestions=first.guess[-1,]) # the [-1,] is necessary because for some strange reason
                                                                       # the suggestions must be smaller, than the population.
                                    }}
 } else {# TRUE if this is not the first run of the GA-optimization
         if (DEBUG) {print(paste("DEBUG: rbga.bin -- run number ",ii))}
         GAmodel <- rbga.bin(size = chrom.size,
                             popSize = pS,
                             iters = iter,
                             mutationChance=MUTATION,
                             elitism = T,
                             evalFunc = evalFunc,
                             suggestions=GAmodel$population[order(GAmodel$evaluations)[1:(pS-1)],])} # use all but the last
                                                                                                     # individual from previous
                                                                                                     # population. This is for
                                                                                                     # technical reasons only.
                                                                                                     # using all somehow didn't work.
                                                                                                     # A long while I used it with
                                                                                                     # taking only the 9 best models,
                                                                                                     # it was good but lacks any
                                                                                                     # reasonable basis!
if (LOAD | ii != ii.min) {
                          OUT.bm <-  rbind(OUT.bm,as.data.frame(cbind("best"=GAmodel$best,"mean"=GAmodel$mean)))
                  } else {
                          OUT.pop <- list("X" = GAmodel$population)
                          OUT.bm <-  as.data.frame(cbind("best"=GAmodel$best,"mean"=GAmodel$mean))
                         }
if (DEBUG) {print("DEBUG: rbga.bin -- done")}
# -- model bulding is over!

## ploting model and optimization status
## this function also saves the population to a text file!
plot.GA.output(GAmodel)

bm.file <- paste(FILENAME,"_runs_",ii,"_mean_best.txt",sep="")
write.table(OUT.bm,file=bm.file,row.names=FALSE)
}
} # END OF ELSE OF: if (!is.na(use.these.descriptors))
