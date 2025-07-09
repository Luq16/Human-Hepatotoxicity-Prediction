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



#
#  defining SVM settings and the fittness function ...
#   Nov 2014:  including (alpha-version) a cross validation (train + test) in stead of test set
#              evaluation during optimization. Thi is an additional option! Not replacing anythig!
#
#
## creating formula for SVM classification!
form <- paste(colnames(dataset)[ncol(dataset)],"~.")
form
# CLASS.WEIGHT
# if class.weight is defined it could be a numeric vector or a string
# stating that proportion according to the actual distribution should
# be taken
# Class weights are only calculated if svm type is not a regression
cat(svm.type,"\n",grep("regression",svm.type),"\n",length(grep("regression",svm.type)),"\n")
if (length(grep("regression",svm.type))==0)
{
weightTAB <- table(dataset[which(dataset[,2]=="train"&!is.na(dataset[,ncol(dataset)])),ncol(dataset)])
cat("distribution of compound classes in training set:\n",weightTAB,"\n")
cat("in percent:\n",round(weightTAB/sum(weightTAB),digits=3)*100,"\n")
if (any(!is.na(class.weight))) {
    if (all(!is.numeric(class.weight))) {
                  if (class.weight=="prop") {   # class.weight = prop  -- use weights proportional to number of compounds in
                                                #                         training set
                                                class.weight <- sum(weightTAB)/weightTAB
                                             } else {cat("ERROR - no such option as class.weight=",class.weight,"\n")
                                                     cat("        using same weight for all classes\n")
                                                     class.weight <- weightTAB/weightTAB
                                                    }
                  } else { # class.weight is obviously numeric! the given values will be used for the classes
                          if (NROW(class.weight)==NROW(weightTAB)) {
                                                                    class.weight <- (weightTAB/weightTAB)*class.weight
                           } else {cat("ERROR - number of given class weights is not equal to numeber of classes!\n")
                                   cat("      number of classes:",NROW(weightTAB),"     and number of specified weights:",NROW(class.weight),"\n")
                                   cat("         using default (all weights = 1) instead!\n")
                                   class.weight <- weightTAB/weightTAB
                           }
                  }
    } else {class.weight <- weightTAB/weightTAB}
}


################################################ ---- evaluation function
if (length(grep("regression",svm.type))==0) {
if (PROBABILITY) {
## Classification usting the AUC as measure of quality/fitness
cat("Using AUC=(area under ROC curve) to evaluate model quality\n")
if (CROSSVAL) {cat("fittnes of individual is AUC + AUC(",Nfold,"-fold cross val)\n",sep="")
       } else {cat("fittnes of individual is AUC(train) + AUC(test)\n",sep="")}
evalFunc <- function(x) {
         if (sum(x) > 0 ) {
         if (CROSSVAL) {TT.comps=which(dataset[,2]=="train"|dataset[,2]=="test")} else {TT.comps=which(dataset[,2]=="train")}
if (!is.na(NU)) {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,gamma=GAMMA,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),             error=function(x) x="ERROR")}
         } else {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,gamma=GAMMA,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),       error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),                   error=function(x) x="ERROR")}}
         if (!is.character(model)) {
         pred.train <- tryCatch(predict(model,dataset[TT.comps,which(x == 1)+2],probability=T), error=function(x) x="ERROR")
         if (is.character(pred.train)) {cat("ERROR in training set prediction during evaluation of an individual with following descriptors:\n")
                                        print(colnames(dataset)[which(x == 1)+2]); return(0)}
         roc.train <- roc(dataset[TT.comps,ncol(dataset)],attributes(pred.train)$probabilities[,"1"])
         AUC.train <- as.numeric(roc.train$auc)
         
         if (CROSSVAL) {
                        pred.test <- do.SVM.crossval(TT.comps,fold=Nfold,selected.descs=which(x == 1)+2)
                        #print(pred.test)
                } else {
                        pred.test <- tryCatch(predict(model,dataset[which(dataset[,2]=="test"),which(x == 1)+2],probability=T), error=function(x) x="ERROR")
                       }
         if (is.character(pred.test)) {cat("ERROR in test set prediction during evaluation of an individual with following descriptors:\n")
                                       print(colnames(dataset)[which(x == 1)+2]); return(0)}
         roc.test <- roc(dataset[which(dataset[,2]=="test"),ncol(dataset)],attributes(pred.test)$probabilities[,"1"])
         AUC.test <- as.numeric(roc.test$auc)
         if (all(is.na(penalty))) {PEN <- 1} else {
                                                   if (sum(x) > penalty[1]) { PEN <- 1 - ((sum(x)-penalty[1])*penalty[2])^2
                                                                              if (PEN < 0) {PEN=0}
                                                                            } else {PEN <- 1}
                                                   if (DEBUG) {cat("DEBUG: using penalty for model evaluation:\n")
                                                                print(paste("        descriptors:",sum(x),"   penalty PEN:",PEN))
                                                               }
                                              }
         } else {AUC.train <- 0;AUC.test <- 0
                 cat("ERROR while constructing a model with the following descriptors:\n")
                 print(colnames(dataset)[which(x == 1)+2])
                 return(0)}
         if (DEBUG) {cat("DEBUG: returning ",-1*(AUC.train*(2-testweight)+AUC.test*testweight)," * ",PEN,"    descriptors:",sum(x),"\n")}
         return(-1*(AUC.train*(2-testweight)+AUC.test*testweight)*PEN)
         } else {return(0)}
} ## END of evalFunc
} else {
### Classification using ACC as measure of quality/fitness
cat("Using accuracy ACC=#(of correctly classified compounds)/#(compounds) to evaluate model quality\n")
evalFunc <- function(x) {
         if (sum(x) > 0 ) {
if (!is.na(NU)) {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,gamma=GAMMA,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),             error=function(x) x="ERROR")}
         } else {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,gamma=GAMMA,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),       error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),                   error=function(x) x="ERROR")}}
## old version without gamma         if (!is.na(NU)) {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU), error=function(x) x="ERROR")
#                          #                  } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type), error=function(x) x="ERROR")}
         if (!is.character(model)) {
         pred.train <- tryCatch(predict(model,dataset[which(dataset[,2]=="train"),which(x == 1)+2]), error=function(x) x="ERROR")
         if (is.character(pred.train)) {cat("ERROR in training set prediction during evaluation of an individual with following descriptors:\n")
                                        print(colnames(dataset)[which(x == 1)+2]); return(0)}
         Ttrain <- table(dataset[which(dataset[,2]=="train"),ncol(dataset)],pred.train)
         Diagonal.sum <- 0; for (d.sum in 1:ncol(Ttrain)) {Diagonal.sum <- Diagonal.sum + Ttrain[d.sum, d.sum]}
         ACC.train <- Diagonal.sum/sum(Ttrain)
         pred.test <- tryCatch(predict(model,dataset[which(dataset[,2]=="test"),which(x == 1)+2]), error=function(x) x="ERROR")
         if (is.character(pred.test)) {cat("ERROR in test set prediction during evaluation of an individual with following descriptors:\n")
                                       print(colnames(dataset)[which(x == 1)+2]); return(0)}
         Ttest <- table(dataset[which(dataset[,2]=="test"),ncol(dataset)],pred.test)
         Diagonal.sum <- 0; for (d.sum in 1:ncol(Ttest)) {Diagonal.sum <- Diagonal.sum + Ttest[d.sum, d.sum]}
         ACC.test <- Diagonal.sum/sum(Ttest)
         if (all(is.na(penalty))) {PEN <- 1} else {
                                                   if (sum(x) > penalty[1]) { PEN <- 1 - ((sum(x)-penalty[1])*penalty[2])^2
                                                                              if (PEN < 0) {PEN=0}
                                                                            } else {PEN <- 1}
                                                   if (DEBUG) {cat("DEBUG: using penalty for model evaluation:\n")
                                                               cat("        descriptors:",sum(x),"   penalty PEN:",PEN,"\n")
                                                               }
                                              }

#         ACC.test <- (Ttest[1,1]+Ttest[2,2])/sum(Ttest)
         } else {ACC.train <- 0;ACC.test <- 0
                 cat("ERROR while constructing a model with the following descriptors:\n")
                 print(colnames(dataset)[which(x == 1)+2])
                 return(0)}
         if (DEBUG) {cat("DEBUG: returning ",-1*(ACC.train*(2-testweight)+ACC.test*testweight)," * ",PEN,"    descriptors:",sum(x),"\n")}
         return(-1*(ACC.train*(2-testweight)+ACC.test*testweight)*PEN)
         } else {return(0)}
} ### END of function
} ### END of if (PROBABILITY)
} ### END of if (length(grep("regression",svm.type))==0)
#################################################################################################################
#################################################################################################################
if (length(grep("regression",svm.type))>0) {
if (fitness=="RSQUARE") {
cat("Using pearson correlation coefficitent R-squared (R^2) to evaluate model quality\n")
cat("\n")
cat("           /       sum[(x-mean(x))*(y-mean(y))]     \\ 2\n ")
cat("   R^2  = | -------------------------------------- |  \n")
cat("           \\ sqrt(sum[(x-mean(x))^2*(y-mean(y))^2]) /  \n")
cat("\n with x the predicted values and y the true values\n")
                       }
if (fitness=="RMSE") {
cat("Using root mean squared error RMSE to evaluate model quality\n")
cat("\n")
cat("                ------------------\n")
cat("               / ( sum[x-y] )^2   \\ \n ")
cat("   RMSE = _ 2/ ----------------   \n")
cat("            \\/         N    \n")
cat("\n with x the predicted values and y the true values and N the number of value pairs\n")
                       }
cat("\n")

evalFunc <- function(x) {
         if (sum(x) > 0 ) {


if (!is.na(NU)) {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"&!is.na(dataset[,ncol(dataset)])),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,gamma=GAMMA,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"&!is.na(dataset[,ncol(dataset)])),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,tolerance=tolerance,epsilon=epsilon),             error=function(x) x="ERROR")}
         } else {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"&!is.na(dataset[,ncol(dataset)])),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,gamma=GAMMA,tolerance=tolerance,epsilon=epsilon),       error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"&!is.na(dataset[,ncol(dataset)])),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,tolerance=tolerance,epsilon=epsilon),                   error=function(x) x="ERROR")}}

         if (!is.character(model)) {

         pred.train <- tryCatch(predict(model,dataset[which(dataset[,2]=="train"),which(x == 1)+2]), error=function(x) x="ERROR")
         if (is.character(pred.train)) {cat("ERROR in training set prediction during evaluation of an individual with following descriptors:\n")
                                        print(colnames(dataset)[which(x == 1)+2]); return(0)}
         RMS.train=sqrt(sum((dataset[which(dataset[,2]=="train"),ncol(dataset)]-pred.train)^2)/NROW(pred.train))
         RMS.train
         pred.diff <- pred.train-mean(pred.train)
         true.diff <- dataset[which(dataset[,2]=="train"),ncol(dataset)]- mean(dataset[which(dataset[,2]=="train"),ncol(dataset)])
         Rsquare.train=(sum(pred.diff*true.diff)/sqrt(sum(pred.diff^2)*sum(true.diff^2)))^2
         Rsquare.train

         pred.test <- tryCatch(predict(model,dataset[which(dataset[,2]=="test"),which(x == 1)+2]), error=function(x) x="ERROR")
         if (is.character(pred.test)) {print("ERROR in test set prediction during evaluation of an individual with following descriptors:")
                                        print(colnames(dataset)[which(x == 1)+2]); return(0)}
         RMS.test=sqrt(sum((dataset[which(dataset[,2]=="test"),ncol(dataset)]-pred.test)^2)/NROW(pred.test))
         RMS.test
         pred.diff <- pred.test-mean(pred.test)
         true.diff <- dataset[which(dataset[,2]=="test"),ncol(dataset)]- mean(dataset[which(dataset[,2]=="test"),ncol(dataset)])
         Rsquare.test=(sum(pred.diff*true.diff)/sqrt(sum(pred.diff^2)*sum(true.diff^2)))^2
         Rsquare.test

         if (all(is.na(penalty))) {PEN <- 1} else {
                                                   if (sum(x) > penalty[1]) { PEN <- 1 - ((sum(x)-penalty[1])*penalty[2])^2
                                                                              if (PEN < 0) {PEN=0}
                                                                            } else {PEN <- 1}
                                                   if (DEBUG) {print("DEBUG: using penalty for model evaluation:")
                                                                print(paste("        descriptors:",sum(x),"   penalty PEN:",PEN))
                                                               }
                                              }
#         ACC.test <- (Ttest[1,1]+Ttest[2,2])/sum(Ttest)
         } else {Rsquare.train <- 0;Rsquare.test <- 0; RMS.train <- NA; RMS.test <- NA;
                 print("ERROR while constructing a model with the following descriptors:")
                 print(colnames(dataset)[which(x == 1)+2])
                 return(0)}

         if (DEBUG) {print(paste("DEBUG: r^2(train) =",Rsquare.train,"; r^2(test) =",Rsquare.test,"; RMS(train) =",RMS.train,"; RMS(test) =",RMS.test))
                      print(paste("DEBUG: returning ",-1*(Rsquare.train*(2-testweight)+Rsquare.test*testweight)," * ",PEN,"    descriptors:",sum(x)))}
         if (fitness=="RMSE") return((RMS.train*(2-testweight)+RMS.test*testweight)*PEN)
         return(-1*(Rsquare.train*(2-testweight)+Rsquare.test*testweight)*PEN)
         } else {; return(0)}
} ### END of function
} ### END of if (length(grep("regression",svm.type))>0)
#PROBABILITY=FALSE
#evalFunc(c(1,1,1,1,1,1))
#evalFuncReg(c(rep(0,100),1,1,1,1,1,1,1,1,1,1))

if (PROBABILITY) {
evalFunc.name <- "AUC(train)+AUC(test)"; func.prefix <- "tpt"; N.chrom.corr <- "unrestricted"
} else {
evalFunc.name <- "ACC(train)+ACC(test)"; func.prefix <- "tpt"; N.chrom.corr <- "unrestricted"
}
if (length(grep("regression",svm.type))>0) {if (fitness=="RMSE") {evalFunc.name <- "RMSE(train)+RMSE(test)"
                                                          } else {evalFunc.name <- "R^2(train)+R^2(test)"}}
if (testweight != 1) {
                      func.prefix <- "tptw"
                      evalFunc.name <- paste("ACC(train)*",round(2-testweight,digits=2),"+ACC(test)*",round(testweight,digits=2),sep="")
                     }
#evalFunc.name <- "ACC(train)+ACC(test)"; func.prefix <- "tptN50"; N.chrom.corr <- "<=50"


