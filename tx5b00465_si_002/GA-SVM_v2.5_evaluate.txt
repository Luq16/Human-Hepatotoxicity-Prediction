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

library(gridExtra)

if (OPT) {l=ii}

cat("\n\n***** --------------     EVALUATION of OPTIMIZED model ------***\n\n")
cat("sorting population according to score\n")
# sorting population
if (exists("GAmodel")) {popul <- data.frame(GAmodel$evaluations,GAmodel$population)}
popul <- popul[order(popul[,1]),]
rownames(popul) <- 1:NROW(popul)
cat("extracting model",MODEL,"\n")
## deleting all descriptors not used in the extracted model
## the extracted model is defined by the variable MODEL, giving the number
## in the sorted/ordered population. Default is 1 (best model)
#selected.descriptors.names <- colnames(dataset)[which(popul[MODEL,]==1)+1]
dataset <- dataset[,c(1,2,which(popul[MODEL,]==1)+1,ncol(dataset))]
selected.descriptors.names <- colnames(dataset)[-c(1,2,ncol(dataset))]
ncol(dataset)
cat("\nuse descriptors :----------------------------\n")
print(selected.descriptors.names)

print("building model")

if (CROSSVAL) {TT.comps=which(dataset[,2]=="train"|dataset[,2]=="test")} else {TT.comps=which(dataset[,2]=="train")}
#model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type), error=function(x) x="ERROR")
#if (is.na(NU)) {model <- svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),-1*c(1,2)], type=svm.type,kernel=kernel.type)
                   #        } else {model <-          svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),-1*c(1,2)], type=svm.type,kernel=kernel.type,nu=NU)}
if (!is.na(NU)) {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,-1*c(1,2)], type=svm.type,kernel=kernel.type,nu=NU,gamma=GAMMA,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,-1*c(1,2)], type=svm.type,kernel=kernel.type,nu=NU,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),             error=function(x) x="ERROR")}
         } else {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,-1*c(1,2)], type=svm.type,kernel=kernel.type,gamma=GAMMA,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),       error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,-1*c(1,2)], type=svm.type,kernel=kernel.type,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),                   error=function(x) x="ERROR")}}



cat("saving model\n")
model.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,".RData",sep="") ## old file name definition: paste(func.prefix,"_popSize_",pS,"_iter_",iter,postfix,"_runs_",l,"_model_",MODEL,".RData",sep="")
model[["DESCRIPTORS"]] <- selected.descriptors.names
model[["DATASET"]] <- dataset
save(model,file=model.file)

cat("caluculation statistics and saving to _stats.txt file\n")
stats.out <- model.stats(model,dataset)
stats.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,"_stats.txt",sep="") ## old file name definition: paste(func.prefix,"_popSize_",pS,"_iter_",iter,postfix,"_runs_",l,"_model_",MODEL,"_stats.txt",sep="")
write.table(stats.out,file=stats.file,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

#---------------------------------------------------------------------------------------------------------
## constructing data.frames with all predictions...
pred <- predict(model,dataset[TT.comps,c(-1,-2,-ncol(dataset))],probability=T)
TRAIN <- data.frame(dataset[TT.comps,c(1:2,ncol(dataset))],"pred"=as.data.frame(pred),"positive.probability"=attributes(pred)$probabilities[,"1"])
if (CROSSVAL) {test.title = paste(Nfold,"-fold cross val.",sep="")
               TEST.comps <- TT.comps
               pred <- do.SVM.crossval(TEST.comps,fold=Nfold,selected.descs=3:(ncol(dataset)-1))
       } else {test.title = "test set"
               TEST.comps <- which(dataset[,2]=="test")
               pred <- predict(model,dataset[TEST.comps,c(-1,-2,-ncol(dataset))],probability=T)
              }
TEST <- data.frame(dataset[TEST.comps,c(1:2,ncol(dataset))],"pred"=as.data.frame(pred),"positive.probability"=attributes(pred)$probabilities[,"1"])
comps.tt.NA <- which(is.na(dataset[,2]))
REST.usable <- as.numeric(rownames(na.omit(dataset[comps.tt.NA,c(1,3:(ncol(dataset)-1))])))

if (NROW(REST.usable) > 0)
{
 pred <- predict(model,dataset[REST.usable,c(-1,-2,-ncol(dataset))],probability=T)
 REST <- data.frame(dataset[REST.usable,c(1:2,ncol(dataset))],"pred"=as.data.frame(pred),"positive.probability"=attributes(pred)$probabilities[,"1"])
# REST
        if (DoVal == 0) {
                         ALL <- data.frame(rbind(TRAIN,TEST,REST))
         } else {        pred <- predict(model,dataset[which(dataset[,2]=="val"),c(-1,-2,-ncol(dataset))],probability=T)
                         VAL <- data.frame(dataset[which(dataset[,2]=="val"),c(1:2,ncol(dataset))],"pred"=as.data.frame(pred),"positive.probability"=attributes(pred)$probabilities[,"1"])
                         ALL <- data.frame(rbind(TRAIN,TEST,VAL,REST))}
} else {
        if (DoVal == 0) {ALL <- data.frame(rbind(TRAIN,TEST))
         } else {        pred <- predict(model,dataset[which(dataset[,2]=="val"),c(-1,-2,-ncol(dataset))],probability=T)
                         VAL <- data.frame(dataset[which(dataset[,2]=="val"),c(1:2,ncol(dataset))],"pred"=as.data.frame(pred),"positive.probability"=attributes(pred)$probabilities[,"1"])
                         ALL <- data.frame(rbind(TRAIN,TEST,VAL))}
}
ALL <- ALL[order(ALL[,1]),]
ALL[,1] <- as.character(ALL[,1])
if (FILL.WITH.DUMMY.COMPOUNDS)
{
TT.DATA[,1] <- as.character(TT.DATA[,1])
FULL <- merge(TT.DATA[,1],ALL,by=1,all=TRUE)
} else {FULL <- ALL}
#head(FULL,n=100)
colnames(FULL)[1] <- "ID"
colnames(FULL)[4] <- paste(ENDPOINT,".SVM.",SET.column.header,".",paste(DESCRIPTOR.postfixes,collapse="."),".",sum(popul[1,-1]),"vars.model.",MODEL,sep="")

head(FULL)



## writing out the full prediction table was moved to after the domain analysis to include it !!
###write.table(FULL,file=paste(func.prefix,"_popSize_",pS,"_iter_",iter,postfix,"_runs_",l,"_model_",MODEL,"_pred.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

cat("...done prediction\n")

###################### -------------------- ROC ROC ROC ROC ROC ROC
###################### -------------------- ROC ROC ROC ROC ROC ROC
###################### -------------------- ROC ROC ROC ROC ROC ROC
if (CROSSVAL) {crossval.text=paste("\n",Nfold,"-fold cross val.",sep="")} else {crossval.text=""}
ROC.plot.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,"_ROC.png",sep="")
print("saving ROC plots to")
print(ROC.plot.file)
png(file=ROC.plot.file,width=2000,height=2000,res=250)
par(mfrow=c(2,2))
plot(0,0,bty="n" ,xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),col="white",xaxt="n",yaxt="n")
text(0.5,0.5,paste("models evaluated:\n",FILENAME,"\n run ",l,
                   "\n\n",ENDPOINT,"\n",paste(DESCRIPTOR.postfixes,collapse=" "),"\nNU = ",NU,"\nGAMMA = ",GAMMA,"\n",ncol(dataset)-3," descriptors\n",NROW(model$index)," support vectors\nclasses:",paste(names(class.weight),collapse="/"),
                   "   weights:",paste(round(class.weight,digits=2),collapse="/"),crossval.text,sep=""))
# TRAIN
roc.train <- plot.roc(TRAIN[,3],TRAIN[,5])
title(main="training set\n")
thr <- ci(roc.train, of="thresholds", thresholds="best")
thr.numeric.train <- round(as.numeric(rownames(thr$sensitivity)),digits=4)
text(0.5,0.2,paste("AUC =",round(roc.train$auc,digits=3),"\nbest threshold =",round(thr.numeric.train,digits=1)))
plot(thr)
# TEST
roc.test <- plot.roc(TEST[,3],TEST[,5])
thr <- ci(roc.test, of="thresholds", thresholds="best")
thr.numeric.test <- round(as.numeric(rownames(thr$sensitivity)),digits=4)
text(0.5,0.2,paste("AUC =",round(roc.test$auc,digits=3),"\nbest threshold =",round(thr.numeric.test,digits=1)))
plot(thr)
title(main=paste(test.title,"\n"))
## VAL
roc.val <- plot.roc(VAL[,3],VAL[,5])
text(0.5,0.2,paste("AUC =",round(roc.val$auc,digits=3)))
title(main="validation set set\n")
thr <- ci(roc.val, of="thresholds", thresholds="best")
thr.numeric.val <- round(as.numeric(rownames(thr$sensitivity)),digits=4)
dev.off()

print("save ROC statistics ROC_stats.txt file")
stats.ROC <- data.frame("ENDP"=ENDPOINT,"TT.split"=SET.column.header,"desc"=paste(DESCRIPTOR.postfixes,collapse=" "),
                        "status"=NA,"popSize"=pS,"iter"=iter,"runs"=l,
                        "nu"=NU,"gamma"=GAMMA,"tolerance"=tolerance,"epsilon"=epsilon,
                        "testweight"=testweight,"mutationChance"=MUTATION,
                        "N.desc"=sum(popul[MODEL,-1]),"N.sv"=NROW(model$index),
                        "AUC.train"=round(roc.train$auc,digits=4),"AUC.test"=round(roc.test$auc,digits=4))
if (DoVal==1) {stats.ROC[,ncol(stats.ROC)+1] <- round(roc.val$auc,digits=4)
               colnames(stats.ROC) <- c(colnames(stats.ROC)[1:(ncol(stats.ROC)-1)],"AUC.val")
              }
stats.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,"_ROC_stats.txt",sep="") ## old file name definition: paste(func.prefix,"_popSize_",pS,"_iter_",iter,postfix,"_runs_",l,"_model_",MODEL,"_stats.txt",sep="")
write.table(stats.ROC,file=stats.file,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)


#-----------------------------------------------------------------------------------------------------------------
## SE/SP vs. thresholds PLOT
if (PRINT.SE.SP.THRESHOLDS)
{
ROC.thr.plot.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,"_ROC_thresholds.png",sep="")
print("saving plots threshold vs SE/SP to")
print(ROC.thr.plot.file)
png(file=ROC.thr.plot.file,width=2000,height=2000,res=250)
par(mfrow=c(2,2))
plot(0,0,bty="n" ,xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),col="white",xaxt="n",yaxt="n")
text(0.5,0.5,paste("models evaluated:\n",FILENAME,"\n run ",l,
                   "\n\n",ENDPOINT,"\n",paste(DESCRIPTOR.postfixes,collapse=" "),"\nNU = ",NU,"\nGAMMA = ",GAMMA,"\n",ncol(dataset)-3," descriptors\n",NROW(model$index)," support vectors\nclasses:",paste(names(class.weight),collapse="/"),
                   "   weights:",paste(round(class.weight,digits=2),collapse="/"),crossval.text,sep=""))
## TRAIN
plot(roc.train$thresholds,roc.train$sensitivities,lty=1,type="l",pch=20,col="red",lwd=1.2,ylab="Sensitivity/Speceficity",xlab="threshold")
points(roc.train$thresholds,roc.train$specificities,lty=1,type="l",pch=20,col="blue",lwd=1.2)
SE.SP.dif.min <- order(abs(roc.train$sensitivities - roc.train$specificities))[1:2]
SE.SP.eq.thr <- mean(roc.train$thresholds[SE.SP.dif.min])
SE.SP.eq.thr.train <- SE.SP.eq.thr
SE.SP.eq     <- mean(roc.train$sensitivities[SE.SP.dif.min])
abline(v=SE.SP.eq.thr,lty=3); points(SE.SP.eq.thr,SE.SP.eq,pch=8)
legend("topleft",legend=c("Sensitivity","Speceficity"),col=c("red","blue"),lty=1,bg="white",cex=0.8)
text(SE.SP.eq.thr,0.1,paste("SE = SP =",round(SE.SP.eq,digits=3),"\nat threshold =",round(SE.SP.eq.thr,digits=3)))
title(main="training set")
## TEST
plot(roc.test$thresholds,roc.test$sensitivities,lty=1,type="l",pch=20,col="red",lwd=1.2,ylab="Sensitivity/Speceficity",xlab="threshold")
points(roc.test$thresholds,roc.test$specificities,lty=1,type="l",pch=20,col="blue",lwd=1.2)
SE.SP.dif.min <- order(abs(roc.test$sensitivities - roc.test$specificities))[1:2]
SE.SP.eq.thr <- mean(roc.test$thresholds[SE.SP.dif.min])
SE.SP.eq.thr.test <- SE.SP.eq.thr
SE.SP.eq     <- mean(roc.test$sensitivities[SE.SP.dif.min])
abline(v=SE.SP.eq.thr,lty=3); points(SE.SP.eq.thr,SE.SP.eq,pch=8)
legend("topleft",legend=c("Sensitivity","Speceficity"),col=c("red","blue"),lty=1,bg="white",cex=0.8)
text(SE.SP.eq.thr,0.1,paste("SE = SP =",round(SE.SP.eq,digits=3),"\nat threshold =",round(SE.SP.eq.thr,digits=3)))
title(main=test.title)
## VAL
plot(roc.val$thresholds,roc.val$sensitivities,lty=1,type="l",pch=20,col="red",lwd=1.2,ylab="Sensitivity/Speceficity",xlab="threshold")
points(roc.val$thresholds,roc.val$specificities,lty=1,type="l",pch=20,col="blue",lwd=1.2)
SE.SP.dif.min <- order(abs(roc.val$sensitivities - roc.val$specificities))[1:2]
SE.SP.eq.thr <- mean(roc.val$thresholds[SE.SP.dif.min])
SE.SP.eq.thr.val <- SE.SP.eq.thr
SE.SP.eq     <- mean(roc.val$sensitivities[SE.SP.dif.min])
abline(v=SE.SP.eq.thr,lty=3); points(SE.SP.eq.thr,SE.SP.eq,pch=8)
legend("topleft",legend=c("Sensitivity","Speceficity"),col=c("red","blue"),lty=1,bg="white",cex=0.8)
text(SE.SP.eq.thr,0.1,paste("SE = SP =",round(SE.SP.eq,digits=3),"\nat threshold =",round(SE.SP.eq.thr,digits=3)))
title(main="validation set set")
dev.off()
} else { #else of:   if (PRINT.SE.SP.THRESHOLDS)
        SE.SP.eq.thr.train=NA
        SE.SP.eq.thr.test=NA
        SE.SP.eq.thr.val=NA
       }
#-----------------------------------------------------------------------------------------------------------------
# DENSITY PLOT
empty <- data.frame()
empty.plot <- ggplot (empty) + geom_blank() + xlim(0, 10) + ylim(0, 10) +
                annotate("text", x = 5, y = 5, label = paste("models evaluated:\n",FILENAME,"\n run ",l,
                   "\n\n",ENDPOINT,"\n",paste(DESCRIPTOR.postfixes,collapse=" "),"\nNU = ",NU,"\nGAMMA = ",GAMMA,"\n",ncol(dataset)-3," descriptors\n",NROW(model$index)," support vectors\nclasses:",paste(names(class.weight),collapse="/"),
                   "   weights:",paste(round(class.weight,digits=2),collapse="/"),crossval.text,sep=""), size=3.5)   +
                theme(panel.background = theme_rect(fill='white', colour='black'))

## TRAIN
dens <- na.omit(data.frame(cbind("class"=TRAIN[,3],"pred"=TRAIN[,5])))
dens[,1] <- as.factor(as.character(dens[,1]))
train.1 <- ggplot(dens,aes(x=pred,fill=class)) + geom_density(alpha=0.3) +  xlab("Poritive Prediciton Probability") +  xlim(min(dens[,2]),max(dens[,2])) +
           ggtitle("training set")
## TEST
dens <- na.omit(data.frame(cbind("class"=TEST[,3],"pred"=TEST[,5])))
dens[,1] <- as.factor(as.character(dens[,1]))
test.1 <- ggplot(dens,aes(x=pred,fill=class)) + geom_density(alpha=0.3) +  xlab("Poritive Prediciton Probability") +  xlim(min(dens[,2]),max(dens[,2])) +
          ggtitle(test.title)
## VAL
dens <- na.omit(data.frame(cbind("class"=VAL[,3],"pred"=VAL[,5])))
dens[,1] <- as.factor(as.character(dens[,1]))
val.1 <- ggplot(dens,aes(x=pred,fill=class)) + geom_density(alpha=0.3) +  xlab("Poritive Prediciton Probability") +  xlim(min(dens[,2]),max(dens[,2])) +
          ggtitle("validation set")
##
dens.plot.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,"_ROC_density.png",sep="")
print("saving density plot to")
print(dens.plot.file)
png(file=dens.plot.file,width=2500,height=2000,res=250)
grid.arrange(arrangeGrob(empty.plot,train.1,test.1,val.1,heights=c(1/2,1/2)),main="denisity plot")
dev.off()

## <------------------------------------------------------------------------------------ work in progress!!!!!
## <------------------------------------------------------------------------------------ work in progress!!!!!
## <------------------------------------------------------------------------------------ work in progress!!!!!
## <------------------------------------------------------------------------------------ work in progress!!!!!
## <------------------------------------------------------------------------------------ work in progress!!!!!
## <------------------------------------------------------------------------------------ work in progress!!!!!
## <------------------------------------------------------------------------------------ work in progress!!!!!
## <------------------------------------------------------------------------------------ work in progress!!!!!
## <------------------------------------------------------------------------------------ work in progress!!!!!
if (DoVal == 1 & PROBABILITY)
{
## plotting densities and enrichment .... !!! this is new since 23. Apr. 2014
dens.enrich.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,"_ROC_dens_enrich.png",sep="")
png(file=dens.enrich.file,height=2000,width=1200,res=250)
dens.prob(ALL[,2],ALL[,3],ALL[,5])
dev.off()
}

################ ---- >> saving model with info on ROC
print("saving model with threshold info")
auc.model <- list("train"=as.numeric(roc.train$auc),"test"=as.numeric(roc.test$auc),"val"=as.numeric(roc.val$auc))
thr.boot.model <- list("train"=thr.numeric.train,"test"=thr.numeric.test,"val"=thr.numeric.val)
thr.SE.SP.eq <- list("train"=SE.SP.eq.thr.train,"test"=SE.SP.eq.thr.test,"val"=SE.SP.eq.thr.val)
model[["ROC"]] <- list("auc"=auc.model,"thr.bootstrap"=thr.boot.model,"thr.SE.SP.eq"=thr.SE.SP.eq)
model[["DESCRIPTORS"]] <- selected.descriptors.names
model.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,".RData",sep="") ## old file name definition: paste(func.prefix,"_popSize_",pS,"_iter_",iter,postfix,"_runs_",l,"_model_",MODEL,".RData",sep="")
save(model,file=model.file)

if (exists("short")) {if (short) {print("only runnig short evaluation!!!");APPLICABILITY=FALSE;PRUNING=FALSE}}



######################################################################################################################
######################################################################################################################
if (APPLICABILITY)
{
print("")
print("            **- assesing applicability domain -**")
print("")
######################################################################################################################
######################################################################################################################
# save the dataset with only the used descriptors so it later can be loaded and a
# applicability domain analysis can be performed!
 print("writing table of compound data and used descriptors for later applicability domain analysis to:")
 AD.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,"_descriptors_for_domain.txt",sep="")
 print(AD.file)
 write.table(dataset,file=AD.file,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

### testing applicability domain
comps.all <- as.numeric(rownames(na.omit(dataset[,-1*c(1,2,ncol(dataset))]))) # comps contains the rownumbers
                                                                              # in the dataframe dataset
                                                                              # of all compounds for which all
                                                                              # used descriptors are defined
comps.train <- which((dataset[,2]=="train") & !is.na(dataset[,ncol(dataset)])) # comps.train contains the rownumbers
                                                                               # in the data frame dataset of all compounds
                                                                               # of the training set used to construct the
                                                                               # model
comps.test <- which((dataset[,2]=="test") & !is.na(dataset[,ncol(dataset)])) # comps.test contains the rownumbers
                                                                             # in the data frame dataset of all compounds
                                                                             # of the test set with endpoint data
comps.non.train <- comps.all[which(!comps.all %in% comps.train)]    # comps.non.train contains the rownumbers
                                                                    # of compounds in comp.all but not in comp.train
if (DoVal == 1) {comps.val <- which((dataset[,2]=="val") & !is.na(dataset[,ncol(dataset)]))}
###head(dataset.norm[comps.non.train,c(1,2,3,4,ncol(dataset.norm))])
###tail(comps.train)
Ntrain <- NROW(comps.train)
Nall <- NROW(comps.all)

#normalize dataset for applicability domain analysis
#normalizing with min and max values from the training set only!!!
#this is necessary to make the normalization consistent!!
dataset.norm <- dataset
norm.min <- apply(dataset.norm[comps.train,-1*c(1,2,ncol(dataset))],2,function(x) min(x,na.rm=TRUE))
norm.max <- apply(dataset.norm[comps.train,-1*c(1,2,ncol(dataset))],2,function(x) max(x,na.rm=TRUE))
print("normalizing descriptors for applicability domain description")
for (i in 3:(ncol(dataset.norm)-1)) {dataset.norm[,i] <- (dataset.norm[,i]-norm.min[i-2])/(norm.max[i-2]-norm.min[i-2])}
# this older version is faster but uses the min and max of the whole data range
# for normalization and not onyl the range from the training set:
#### dataset.norm[,3:(ncol(dataset)-1)] <- apply(dataset[,-1*c(1,2,ncol(dataset))],2,function(x) (x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))
## head(dataset.norm)


DISTANCE <- dist(dataset.norm[comps.all,-1*c(1,2,ncol(dataset.norm))],method="euclidean") # this is the distance matrix
                                                                                # of all compounds for which a distance can be
                                                                                # calculated, i.e. all descriptors are defined
                                                                                # for these compounds

DM.train <- as.matrix(DISTANCE)[1:Ntrain,1:Ntrain]        # DM.train is the part of the distance matrix
                                                          # containing only the distances from compounds inside
                                                          # the trainig set to eachother

DM.rest  <- as.matrix(DISTANCE)[1:Ntrain,(Ntrain+1):Nall] # DM.rest is the part of the distance matrix
                                                          # containing all distances from compounds not in the training set to
                                                          # all compounds in the training set, excluding distances from
                                                          # compounds outside to eachother as well as from compounds in
                                                          # the training to eachother


# calculating the threshold for minimal distance to applicability domain
DM.train.vector <- DM.train[which(DM.train!=0)]
DM.train.mean <- mean(DM.train[which(DM.train!=0)])
D.train.little.mean <- mean(DM.train.vector[which(DM.train.vector < DM.train.mean)])
D.train.little.sd   <-   sd(DM.train.vector[which(DM.train.vector < DM.train.mean)])
APD <- D.train.little.mean + 0.5*D.train.little.sd
APD
print("Applicability domain threshold for minimum distance to a compound in training set:")
print(paste("APD = <d> + Z*s =",round(D.train.little.mean,digits=2)," + 0.5 * ",round(D.train.little.sd,digits=2)," = ",round(APD,digits=2)))
print("with ")
print(" <d> = mean(D_little)  and s = sd(D_little) and Z = 0.5")
print(" and D_little is the set of distances D with D < mean(D), mean(D) is the mean of all distances")
print(" All distances are euclidean distances using normalized descriptors - normalized with min and ")
print(" max values from the training set only, so they still can be <0 or >1. This")
print(" might differ from the originally proposed distances when using this threshold.")
print("See: Tropsha et al. J. Chem. Inf. Model. 46 (2006), pp. 1984-1995")

## check which of the compounds is outside the applicability domain
DD <- data.frame("ID"=dataset[,1],"tt"=dataset[,2],"ENDP"=dataset[,ncol(dataset)],"Dmin"=0,"Dmax"=0,"Dmean"=0)
DD[comps.non.train,"Dmin"]  <- apply(DM.rest,2,function(x)  min(x[x!=0]))
DD[comps.non.train,"Dmax"]  <- apply(DM.rest,2,function(x)  max(x[x!=0]))
DD[comps.non.train,"Dmean"] <- apply(DM.rest,2,function(x) mean(x[x!=0]))
DD[comps.train,"Dmin"]  <- apply(DM.train,2,function(x)  min(x[x!=0]))
DD[comps.train,"Dmax"]  <- apply(DM.train,2,function(x)  max(x[x!=0]))
DD[comps.train,"Dmean"] <- apply(DM.train,2,function(x) mean(x[x!=0]))
DD[,ncol(DD)+1] <- NA
AD.col.name <- paste("applicability.domain.APD.",round(APD,digits=2),sep="")
colnames(DD)[ncol(DD)] <- AD.col.name
DD[which(DD[,"Dmin"] < APD),AD.col.name] <- "in"
DD[which(DD[,"Dmin"] >= APD),AD.col.name] <- "out"


N.train.out <- NROW(which(DD[comps.train,"Dmin"]>APD))
N.test.out  <- NROW(which(DD[comps.test,"Dmin"]>APD))
if (DoVal == 1) {N.val.out  <- NROW(which(DD[comps.val,"Dmin"]>APD))}
if (N.train.out > 0) {
                      print(paste(N.train.out,"compounds in the training set are actually outside it's applicability domain"))
                      print("and have a minimal distance Dmin >= APD")
                      print("these compounds are:")
                      tmp.OUT <- merge(DD[comps.train[which(DD[comps.train,"Dmin"]>APD)],-3],FULL[comps.train[which(DD[comps.train,"Dmin"]>APD)],-2])
                      colnames(tmp.OUT)[ncol(tmp.OUT)] <- "predicted"
                      print(tmp.OUT)
              } else {
                      print("All compounds in the training set have a minimal distance Dmin < APD")
                     }
if (N.test.out > 0) {
                      print(paste(N.test.out,"compounds in the test set are outside the applicability domain"))
                      print("and have a minimal distance Dmin >= APD")
                      print("these compounds are:")
                      tmp.OUT <- merge(DD[comps.test[which(DD[comps.test,"Dmin"]>APD)],-3],FULL[comps.test[which(DD[comps.test,"Dmin"]>APD)],-2])
                      colnames(tmp.OUT)[ncol(tmp.OUT)] <- "predicted"
                      print(tmp.OUT)
              } else {
                      print("All compounds in the test set are in the applicability domain and have a minimal distance Dmin < APD")
                     }
if (DoVal == 1) {
        if (N.val.out > 0) {
                      print(paste(N.val.out,"compounds in the validation set are outside the applicability domain"))
                      print("and have a minimal distance Dmin >= APD")
                      print("these compounds are:")
                      tmp.OUT <- merge(DD[comps.val[which(DD[comps.val,"Dmin"]>APD)],-3],FULL[comps.val[which(DD[comps.val,"Dmin"]>APD)],-2])
                      colnames(tmp.OUT)[ncol(tmp.OUT)] <- "predicted"
                      print(tmp.OUT)
              } else {
                      print("All compounds in the validation set are in the applicability domain and have a minimal distance Dmin < APD")
                     }
                }
# ---> only save applicability domain in/out
if(FALSE) {FULL <- merge(FULL,DD[,c(1,ncol(DD))])}
# ---> save applicability domain in/out and distance minimum mean and maximum
FULL <- merge(FULL,DD[,c("ID","Dmin","Dmax","Dmean",AD.col.name)])
head(FULL)

} ## END of if (APPLICABILITY)

print("writing table with all predicitons and applicability domain analysis to :")
FULL.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,"_pred.txt",sep="") ## old file name definition: paste(func.prefix,"_popSize_",pS,"_iter_",iter,postfix,"_runs_",l,"_model_",MODEL,"_pred.txt",sep="")
print(FULL.file)
write.table(FULL,file=FULL.file,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)


##
##
## how to calculat fast the distance between training set and all test set compounds??
## something like a distance matrix with all train and test data. Train and test are
## ordere and 1:N is train, and (N+1):M test. and in the distance matrix only
## the elements [(N+1):M,1:N] are needed!
## test:
#X <- data.frame("A"=1:10,"B"=rep(1,10))
#DX <- dist(X,method="euclidean")
#DX
#use <- c(3,4,1,2,5,6,7,8,9,10)
#use
#X[use,]
#DX <- dist(X[use,],method="euclidean")
#DXm <- as.matrix(DX)
#DXm[1:2,3:10]

comps.train <- which((dataset[,2]=="train") & !is.na(dataset[,ncol(dataset)])) # comps.train contains the rownumbers
                                                                               # in the data frame dataset of all compounds
                                                                               # of the training set used to construct the
                                                                               # model
if (CROSSVAL) {comps.train <- which((dataset[,2]=="train"|dataset[,2]=="test")& !is.na(dataset[,ncol(dataset)]))}
comps.test <- which((dataset[,2]=="test") & !is.na(dataset[,ncol(dataset)])) # comps.test contains the rownumbers
                                                                             # in the data frame dataset of all compounds
                                                                             # of the test set with endpoint data
if (CROSSVAL) {comps.test <- which((dataset[,2]=="train"|dataset[,2]=="test")& !is.na(dataset[,ncol(dataset)]))}
if (DoVal) {
comps.val <- which((dataset[,2]=="val") & !is.na(dataset[,ncol(dataset)])) # comps.val contains the rownumbers
                                                                             # in the data frame dataset of all compounds
                                                                             # of the validation set with endpoint data
           } else {comps.val=NULL}

if (CROSSVAL) {TEST.name <- " CV "} else {TEST.name <- "test"}
      cat("... evaluation of first model done!\n")
      cat("#########################################\n")
      cat("#                                       #\n")
      cat("#         the Model                     #\n")
      cat("#                                       #\n")
      cat("#########################################\n")
cat("#   Settings:","\n")
cat("#        nu        =",NU,"\n")
cat("#        gamma     =",GAMMA,"\n")
cat("#        tolerance =",tolerance,"\n")
cat("#        epsilon   =",epsilon,"\n")
cat("#   Dataset:","\n")
cat("#        N(train)  =",NROW(comps.train),"\n")
cat("#          N(pos)  =",NROW(which(dataset[comps.train,ncol(dataset)]==1)),"\n")
cat("#          N(neg)  =",NROW(which(dataset[comps.train,ncol(dataset)]==0)),"\n")
if (CROSSVAL) {cat("#      ",Nfold,"-fold cross validation\n",sep="")}
cat("#        N(",TEST.name,")   = ",NROW(TEST.comps),"\n",sep="")
cat("#        N(val)    =",NROW(comps.val),"\n")
cat("#   Model:","\n")
cat("#      General characteristics:","\n")
cat("#        N(descriptors)     =",sum(popul[MODEL,-1]),"\n")
cat("#        N(support vectors) =",NROW(model$index),"\n")
if (PROBABILITY) {
cat(      "#      ROC analysis:\n")
cat("#        AUC(train)         =",round(roc.train$auc,digits=3),"\n")
cat("#        AUC(",TEST.name,")          = ",round(roc.test$auc,digits=3),"\n",sep="")
if (DoVal==1) {
cat("#        AUC(val)           =",round(roc.val$auc,digits=3),"\n")
}}
if (PROBABILITY) {
cat(      "#      Bootstrap thresholds from ROC curves:\n")
cat("#        threshold(train)   =",thr.numeric.train,"\n")
cat("#        threshold(",TEST.name,")   = ",thr.numeric.test,"\n")
if (DoVal==1) {
cat("#        threshold(val)     =",thr.numeric.val,"\n")
}}
cat("\n")
      cat("#########################################\n")






if (PRUNING)
{
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
print("")
print("            **- Starting model pruning -**")
print("")
#############################################################################################################
#############################################################################################################
## FUNCTION FOR EVALUATION - FUNCTION FOR EVALUATION - FUNCTION FOR EVALUATION - FUNCTION FOR EVALUATION
## FUNCTION FOR EVALUATION - FUNCTION FOR EVALUATION - FUNCTION FOR EVALUATION - FUNCTION FOR EVALUATION
############################################### ---- evaluation function
if (PROBABILITY) {
print("Using AUC=(area under ROC curve) to evaluate model quality")
evalFunc.2 <- function(x) {
         if (sum(x) > 0 ) {
if (!is.na(NU)) {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,gamma=GAMMA,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),             error=function(x) x="ERROR")}
         } else {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,gamma=GAMMA,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),       error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),                   error=function(x) x="ERROR")}}
         if (!is.character(model)) {
         pred.train <- tryCatch(predict(model,dataset[which(dataset[,2]=="train"),which(x == 1)+2],probability=T), error=function(x) x="ERROR")
         if (is.character(pred.train)) {print("ERROR in training set prediction during evaluation of an individual with following descriptors:")
                                        print(colnames(dataset)[which(x == 1)+2]); return(0)}
         roc.train <- roc(dataset[which(dataset[,2]=="train"),ncol(dataset)],attributes(pred.train)$probabilities[,"1"])
         AUC.train <- as.numeric(roc.train$auc)
         pred.test <- tryCatch(predict(model,dataset[which(dataset[,2]=="test"),which(x == 1)+2],probability=T), error=function(x) x="ERROR")
         if (is.character(pred.test)) {print("ERROR in test set prediction during evaluation of an individual with following descriptors:")
                                        print(colnames(dataset)[which(x == 1)+2]); return(0)}
         roc.test <- roc(dataset[which(dataset[,2]=="test"),ncol(dataset)],attributes(pred.test)$probabilities[,"1"])
         AUC.test <- as.numeric(roc.test$auc)
         if (DoVal == 1) {
                          pred.val <- tryCatch(predict(model,dataset[which(dataset[,2]=="val"),which(x == 1)+2],probability=T), error=function(x) x="ERROR")
                          if (is.character(pred.val)) {print("ERROR in test set prediction during evaluation of an individual with following descriptors:")
                                                       print(colnames(dataset)[which(x == 1)+2]); return(0)}
                          roc.val <- roc(dataset[which(dataset[,2]=="val"),ncol(dataset)],attributes(pred.val)$probabilities[,"1"])
                          AUC.val <- as.numeric(roc.val$auc)
                         }
         } else {AUC.train <- 0;AUC.test <- 0
                 print("ERROR while constructing a model with the following descriptors:")
                 print(colnames(dataset)[which(x == 1)+2])
                 return(0)}
         if (DoVal == 1) {
                          return(c(-1*(AUC.train*(2-testweight)+AUC.test*testweight),AUC.train,AUC.test,AUC.val))
                   } else {
                          return(c(-1*(AUC.train*(2-testweight)+AUC.test*testweight),AUC.train,AUC.test))
                          }
         } else {return(0)}
}
} else {
print("Using accuracy ACC=#(of correctly classified compounds)/#(compounds) to evaluate model quality")
evalFunc.2 <- function(x) {
         if (sum(x) > 0 ) {
if (!is.na(NU)) {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,gamma=GAMMA,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")}
         } else {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,gamma=GAMMA,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")}}
## old version without gamma         if (!is.na(NU)) {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU), error=function(x) x="ERROR")
#                          #                  } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type), error=function(x) x="ERROR")}
         if (!is.character(model)) {
         pred.train <- tryCatch(predict(model,dataset[which(dataset[,2]=="train"),which(x == 1)+2]), error=function(x) x="ERROR")
         if (is.character(pred.train)) {print("ERROR in training set prediction during evaluation of an individual with following descriptors:")
                                        print(colnames(dataset)[which(x == 1)+2]); return(0)}
         Ttrain <- table(dataset[which(dataset[,2]=="train"),ncol(dataset)],pred.train)
         Diagonal.sum <- 0; for (d.sum in 1:ncol(Ttrain)) {Diagonal.sum <- Diagonal.sum + Ttrain[d.sum, d.sum]}
         ACC.train <- Diagonal.sum/sum(Ttrain)
         pred.test <- tryCatch(predict(model,dataset[which(dataset[,2]=="test"),which(x == 1)+2]), error=function(x) x="ERROR")
         if (is.character(pred.test)) {print("ERROR in test set prediction during evaluation of an individual with following descriptors:")
                                        print(colnames(dataset)[which(x == 1)+2]); return(0)}
         Ttest <- table(dataset[which(dataset[,2]=="test"),ncol(dataset)],pred.test)
         Diagonal.sum <- 0; for (d.sum in 1:ncol(Ttest)) {Diagonal.sum <- Diagonal.sum + Ttest[d.sum, d.sum]}
         ACC.test <- Diagonal.sum/sum(Ttest)
         ACC.test <- (Ttest[1,1]+Ttest[2,2])/sum(Ttest)
         if (DoVal == 1)
         {
           Tval <- table(dataset[which(dataset[,2]=="val"),ncol(dataset)],predict(model,dataset[which(dataset[,2]=="val"),which(x == 1)+2]))
           Diagonal.sum <- 0; for (d.sum in 1:ncol(Tval)) {Diagonal.sum <- Diagonal.sum + Tval[d.sum, d.sum]}
           ACC.val <- (Diagonal.sum/sum(Tval))
         }
         } else {ACC.train <- 0;ACC.test <- 0
                 print("ERROR while constructing a model with the following descriptors:")
                 print(colnames(dataset)[which(x == 1)+2])
                 return(0)}
         if (DoVal == 1) {
                          return(c(-1*(ACC.train*(2-testweight)+ACC.test*testweight),ACC.train,ACC.test,ACC.val))
                  } else {
                          return(c(-1*(ACC.train*(2-testweight)+ACC.test*testweight),ACC.train,ACC.test))
                         }
         } else {return(0)}
}
} ### END of if (PROBABILITY)
#############################################################################################################
evalFunc.pruned <- function(x) {
         if (sum(x) > 0 ) {
if (!is.na(NU)) {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,gamma=GAMMA,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),             error=function(x) x="ERROR")}
         } else {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,gamma=GAMMA,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),       error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),                   error=function(x) x="ERROR")}}
         if (!is.character(model)) {
         Ttrain <- table(dataset[which(dataset[,2]=="train"),ncol(dataset)],predict(model,dataset[which(dataset[,2]=="train"),which(x == 1)+2]))
         Ttest <- table(dataset[which(dataset[,2]=="test"),ncol(dataset)],predict(model,dataset[which(dataset[,2]=="test"),which(x == 1)+2]))
         if (DoVal == 1)
         {
           Tval <- table(dataset[which(dataset[,2]=="val"),ncol(dataset)],predict(model,dataset[which(dataset[,2]=="val"),which(x == 1)+2]))
         }
         } else {if (DEBUG) {print("DEBUG: error while building the model with the gen:");print(x)} }
         if (DoVal == 1) {
                          return(list("Ttrain"=Ttrain,"Ttest"=Ttest,"Tval"=Tval))
                  } else {
                          return(list("Ttrain"=Ttrain,"Ttest"=Ttest))
                         }
         } else {return(0)}
}
################################################ ---- evaluation function
if (FALSE) {
if (PROBABILITY) {
print("Using AUC=(area under ROC curve) to evaluate model quality")
evalFunc.pruned <- function(x) {
         if (sum(x) > 0 ) {
if (!is.na(NU)) {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,gamma=GAMMA,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),             error=function(x) x="ERROR")}
         } else {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,gamma=GAMMA,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),       error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),                   error=function(x) x="ERROR")}}
         if (!is.character(model)) {
         pred.train <- tryCatch(predict(model,dataset[which(dataset[,2]=="train"),which(x == 1)+2],probability=T), error=function(x) x="ERROR")
         if (is.character(pred.train)) {print("ERROR in training set prediction during evaluation of an individual with following descriptors:")
                                        print(colnames(dataset)[which(x == 1)+2]); return(0)}
         roc.train <- roc(dataset[which(dataset[,2]=="train"),ncol(dataset)],attributes(pred.train)$probabilities[,"1"])
         AUC.train <- as.numeric(roc.train$auc)
         pred.test <- tryCatch(predict(model,dataset[which(dataset[,2]=="test"),which(x == 1)+2],probability=T), error=function(x) x="ERROR")
         if (is.character(pred.test)) {print("ERROR in test set prediction during evaluation of an individual with following descriptors:")
                                        print(colnames(dataset)[which(x == 1)+2]); return(0)}
         roc.test <- roc(dataset[which(dataset[,2]=="test"),ncol(dataset)],attributes(pred.test)$probabilities[,"1"])
         AUC.test <- as.numeric(roc.test$auc)
         } else {AUC.train <- 0;AUC.test <- 0
                 print("ERROR while constructing a model with the following descriptors:")
                 print(colnames(dataset)[which(x == 1)+2])
                 return(0)}
         return(-1*(AUC.train*(2-testweight)+AUC.test*testweight))
         } else {return(0)}
}
} else {
print("Using accuracy ACC=#(of correctly classified compounds)/#(compounds) to evaluate model quality")
evalFunc.pruned <- function(x) {
         if (sum(x) > 0 ) {
if (!is.na(NU)) {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,gamma=GAMMA,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")}
         } else {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,gamma=GAMMA,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")}}
## old version without gamma         if (!is.na(NU)) {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU), error=function(x) x="ERROR")
#                          #                  } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(x == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type), error=function(x) x="ERROR")}
         if (!is.character(model)) {
         pred.train <- tryCatch(predict(model,dataset[which(dataset[,2]=="train"),which(x == 1)+2]), error=function(x) x="ERROR")
         if (is.character(pred.train)) {print("ERROR in training set prediction during evaluation of an individual with following descriptors:")
                                        print(colnames(dataset)[which(x == 1)+2]); return(0)}
         Ttrain <- table(dataset[which(dataset[,2]=="train"),ncol(dataset)],pred.train)
         Diagonal.sum <- 0; for (d.sum in 1:ncol(Ttrain)) {Diagonal.sum <- Diagonal.sum + Ttrain[d.sum, d.sum]}
         ACC.train <- Diagonal.sum/sum(Ttrain)
         pred.test <- tryCatch(predict(model,dataset[which(dataset[,2]=="test"),which(x == 1)+2]), error=function(x) x="ERROR")
         if (is.character(pred.test)) {print("ERROR in test set prediction during evaluation of an individual with following descriptors:")
                                        print(colnames(dataset)[which(x == 1)+2]); return(0)}
         Ttest <- table(dataset[which(dataset[,2]=="test"),ncol(dataset)],pred.test)
         Diagonal.sum <- 0; for (d.sum in 1:ncol(Ttest)) {Diagonal.sum <- Diagonal.sum + Ttest[d.sum, d.sum]}
         ACC.test <- Diagonal.sum/sum(Ttest)
#         ACC.test <- (Ttest[1,1]+Ttest[2,2])/sum(Ttest)
         } else {ACC.train <- 0;ACC.test <- 0
                 print("ERROR while constructing a model with the following descriptors:")
                 print(colnames(dataset)[which(x == 1)+2])
                 return(0)}
         return(-1*(ACC.train*(2-testweight)+ACC.test*testweight))
         } else {return(0)}
}
} ### END of if (PROBABILITY)
} ## END if (FALSE)
###################################################################################################################
## FUNCTION FOR PLOTING - FUNCTION FOR PLOTING - FUNCTION FOR PLOTING - FUNCTION FOR PLOTING - FUNCTION FOR PLOTING
## FUNCTION FOR PLOTING - FUNCTION FOR PLOTING - FUNCTION FOR PLOTING - FUNCTION FOR PLOTING - FUNCTION FOR PLOTING
###################################################################################################################
plot.eval <- function(Tab,plot.title="") {
  NR <- nrow(Tab)
  NC <- ncol(Tab)
  # adjust size of Text to not overlap so much
  Text.size <- 1.6
  if (NC > 3 | NR > 3) {Text.size <- 1.1}
  if (NC > 5 | NR > 5) {Text.size <- 0.8}
#
  plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,0),ylim=c(0,0.5),main=plot.title,cex=Text.size);
  x.pos <- 0.38
  y.pos <- 0.15
  x.size <- 0.20
  y.size <- 0.20
  if (NC > 5 | NR > 5) {x.pos <- 0.4; y.pos <- 0.09
                        x.size <- 0.3; y.size <- 0.3}
  ## SE is a vector of length NR <- nrow(Tab)
  ## for a 2x2 table SE[2] == SP
  SE <- rep(0,NC)
  for (i in 1:NC) {if (i <= NR) {SE[i] <- Tab[i,i]/sum(Tab[i,])}}
  #
  PPV <- rep(0,NR)
  for (i in 1:NR) {if (i <= NC) {PPV[i] <- Tab[i,i]/sum(Tab[,i])}}
  #
  Diagonal.sum <- 0; for (d.sum.r in 1:NR) {for (d.sum.c in 1:NC) {if (d.sum.r == d.sum.c) {Diagonal.sum <- Diagonal.sum + Tab[d.sum.r, d.sum.c]}}}
  ACC <- Diagonal.sum/sum(Tab)
#  text(x.pos-0.5*x.size,y.pos+y.size+0.15,paste("thr. = NA"))
  rect(x.pos,y.pos,x.pos-x.size,y.pos+y.size) # inner box
  rect(x.pos,y.pos+y.size,x.pos-x.size,y.pos+y.size+0.1) # top box  -> "predicte 0 1"
  rect(x.pos,y.pos-0.1,x.pos-x.size,y.pos) # bottom box  -> "NPV/PPV"
  rect(x.pos+0.1,y.pos,x.pos,y.pos+y.size) # left box -> "actual 0 1"
  rect(x.pos-x.size,y.pos,x.pos-x.size-0.1,y.pos+y.size) # right box -> "SP / SE"
  rect(x.pos-x.size,y.pos-0.1,x.pos-x.size-0.1,y.pos) # bottom-right box -> "ACC"
#
  text(x.pos-0.5*x.size,y.pos+y.size+0.08,"predicted",cex=Text.size)
  for (j in 1:NC) {text(x.pos-(1+2*(j-1))*x.size/(2*NC),y.pos+y.size+0.03,colnames(Tab)[j],cex=Text.size)}
#
  text(x.pos+0.075,y.pos+0.5*y.size,"actual",srt=90,cex=Text.size)
  for (i in 1:NR) {text(x.pos+0.025,y.pos+y.size-(1+2*(i-1))*y.size/(2*NR),rownames(Tab)[i],cex=Text.size)}
#
  for (i in 1:NC) {for (j in 1:NR) {text(x.pos-(1+2*(i-1))*x.size/(2*NC),y.pos+y.size-(1+2*(j-1))*y.size/(2*NR),Tab[j,i],cex=Text.size)}}
  for (j in 1:NR) {text(x.pos-x.size-0.05,y.pos+y.size-(1+2*(j-1))*y.size/(2*NR),paste(round(SE[j]*100,digits=0),"%",sep=""),cex=Text.size)}
  for (i in 1:NC) {text(x.pos-(1+2*(i-1))*x.size/(2*NC),y.pos-0.05,paste(round(PPV[i]*100,digits=0),"%",sep=""),cex=Text.size)}
  text(x.pos-x.size-0.05,y.pos-0.05,paste("ACC\n",round(ACC*100,digits=0),"%",sep=""),cex=Text.size)
}
#############################################################################################################
#############################################################################################################
#############################################################################################################


## creating the prune gen
desc.prune <- rep(1,ncol(dataset)-3)
desc.prune

## take out every descriptor once
for (j in 0:NROW(desc.prune))
{
 desc.prune <- rep(1,ncol(dataset)-3)
 if (j > 0) {desc.prune[j] <- 0}
 error <- evalFunc.2(desc.prune)
 if (j==0)
 {
  pruning <- data.frame("total.error"=error[1],"ACC.train"=error[2],"ACC.test"=error[3],"ACC.val"=0,"N"=0,"gens"=paste(desc.prune,collapse=","),"left.out.descriptor"="NONE")
  if (DoVal == 1) {pruning[1,4] <- error[4]}
  pruning[,6] <- as.character(pruning[,6])
  pruning[,7] <- as.character(pruning[,7])
 } else {pruning[j+1,1] <- error[1]
         pruning[j+1,2] <- error[2]
         pruning[j+1,3] <- error[3]
         if (DoVal == 1) {pruning[j+1,4] <- error[4]} else {pruning[j+1,4] <- 0}
         pruning[j+1,5] <- 1
         pruning[j+1,6] <- paste(desc.prune,collapse=",")
         pruning[j+1,7] <- colnames(dataset)[which(desc.prune==0)+2]
 }
}
pruning
print("ranking of selected descriptors")
print("from least to most relevant")
if (PROBABILITY)
       {print(paste("total.error = -1*(",2-testweight,"*AUC(train)+",testweight,"*AUC(test))",sep=""))
} else {print(paste("total.error = -1*(",2-testweight,"*ACC(train)+",testweight,"*ACC(test))",sep=""))}
if (DoVal == 1) {
                 print(pruning[order(pruning[,1]),c(1,2,3,4,7)])
         } else {
                 print(pruning[order(pruning[,1]),c(1,2,3,7)])
                }
take.out.order <- order(pruning[-1,1]) # take.out.order contains
                                       # the index numbers of the descriptors (!! not the columnes of the descriptors!!)
                                       # in the order of their (first order) relevance.
                                       # first order relevance is the influence of leaving
                                       # out just this one descriptor.
                                       # --> take.out.order is in the same system as the gens or
                                       #     then variable desc.prune
if (DEBUG) {print("DEBUG: take.out.order");print(take.out.order)}
if (DEBUG) {print("...check point 1")}

#pruning
#i=27
#colnames(dataset)[which(as.integer(unlist(strsplit(pXX[order(-pXX[,1])[i],5],",")))==0)+2]
#colnames(dataset)
#colnames(dataset)
#colnames(dataset)[order(pruning[-1,1])+1]
#
#pXX <- pruning[which(pruning[,4]==1),]
#colnames(dataset)[which(as.integer(unlist(strsplit(pXX[order(-pXX[,1])[1],5],",")))==0)+2]
#colnames(dataset)[which(as.integer(unlist(strsplit(pXX[order(-pXX[,1])[2],5],",")))==0)+2]
#colnames(dataset)[which(as.integer(unlist(strsplit(pXX[order(-pXX[,1])[3],5],",")))==0)+2]

if (DEBUG) {print(paste("ncol(dataset)",ncol(dataset)))
             print(paste("NROW(desc.prune)",NROW(desc.prune)))
             print(paste("NROW(take.out.order)",NROW(take.out.order)))
             }

dj=nrow(pruning)
for (j in 1:(NROW(desc.prune)-2)) # using -2 so the last model will have two descriptors.
                                  # I'm not sure why, but at least in one cases using
                                  # only one descriptor produced some error. Maybe
                                  # it was the type of descriptor... (vsplus: VSDD7)
                                  # but I'm not sure at all!
{
 jj=j+dj
 if (j == 1) {desc.prune <- rep(1,ncol(dataset)-3)}
 desc.prune[take.out.order[j]] <- 0     # setting the gen desc.prune to zero for the descriptor to
                                        # be taken out at cycle j.
###
 if (DEBUG) {print("...loop for (j in 1:(NROW(desc.prune)-1)) .... ");print(paste("j:",j,"    jj:",jj));print(colnames(dataset)[take.out.order[j]+2])}
 error <- evalFunc.2(desc.prune)
 if (NROW(error)==1 & error == 0) { # model could not be built for some reason
                                   if (DEBUG) {print("DEBUG: the model could not be built for some reason or another ...??!!?!?")}
                                   break
                                  }
 pruning[jj,1] <- error[1]
 pruning[jj,2] <- error[2]
 pruning[jj,3] <- error[3]
 if (DoVal == 1) {pruning[jj,4] <- error[4]} else {pruning[jj,4] <- 0}
 pruning[jj,5] <- NROW(desc.prune)-sum(desc.prune)
 pruning[jj,6] <- paste(desc.prune,collapse=",")
 pruning[jj,7] <- paste(colnames(dataset)[which(desc.prune==0)+2],collapse=",")
}
if (DEBUG) {print("...check point 2")}
#pruning
#pruning[take.out.order[1],]
#pruning[c(take.out.order[1],which(pruning[,5]>1)),]
print("writing table with pruning data")
prune.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,"_pruning.txt",sep="") ## old file name definition: paste(func.prefix,"_popSize_",pS,"_iter_",iter,postfix,"_runs_",l,"_model_",MODEL,"_pruning.txt",sep="")
write.table(pruning[c(take.out.order[1],which(pruning[,5]>1)),],file=prune.file,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

if (DoVal == 1) {ACC.val.all <- pruning[1,4]}
ACC.test.all <- pruning[1,3]
ACC.train.all <- pruning[1,2]
total.error.all <- pruning[1,1]

if (NROW(which(pruning[,1]<total.error.all)) > 0) {
   print("hmm.... ")
  } else {print("all pruned models have a larger absolut total error than the unpruned model")}
if (NROW(which(pruning[,2]>ACC.train.all)) > 0) {
   print("pruned models with a larger ACC(train) than the unpruned model:")
   paste("pruned: ", pruning[which(pruning[,2]>ACC.train.all),7])
  } else {print("all pruned models have a smaller ACC(train) than the unpruned model")}
if (NROW(which(pruning[,3]>ACC.test.all)) > 0) {
   print("hmm.... ")
  } else {print("all pruned models have a smaller ACC(test) than the unpruned model")}


## how many descriptors can be pruned so only 5% of ACC(test) or ACC(train) are lost???
pruning$pp.test <- 0
pruning$pp.train <- 0
models.to.eval <- c(take.out.order[1]+1,which(pruning[,5]>1))
#pruning[models.to.eval,c(1,2,3,4,5,7)]

ACC.test.95p.models <- models.to.eval[which(pruning[models.to.eval,3] >=ACC.test.all*0.95)]
pruning[ACC.test.95p.models,"pp.test"] <- 5
ACC.test.95p.max <- ACC.test.95p.models[which(pruning[ACC.test.95p.models,5]==max(pruning[ACC.test.95p.models,5]))]
ACC.test.95p.max.Npruned <- max(pruning[ACC.test.95p.models,5])
ACC.test.95p.max # this contains the number of the model in the pruning data.frame
                 # specifying the model with the maximum number of pruned descriptors
                 # and at the same time a maximum reduction of ACC.test by 5%
ACC.test.95p.max.gen <- as.integer(unlist(strsplit(pruning[ACC.test.95p.max,6],",")))
ACC.test.95p.max.gen

ACC.train.95p.models <- models.to.eval[which(pruning[models.to.eval,2] >=ACC.train.all*0.95)]
pruning[ACC.train.95p.models,"pp.train"] <- 5
ACC.train.95p.max <- ACC.train.95p.models[which(pruning[ACC.train.95p.models,5]==max(pruning[ACC.train.95p.models,5]))]
ACC.train.95p.max.Npruned <- max(pruning[ACC.train.95p.models,5])
ACC.train.95p.max # this contains the number of the model in the pruning data.frame
                 # specifying the model with the maximum number of pruned descriptors
                 # and at the same time a maximum reduction of ACC.train by 5%
ACC.train.95p.max.gen <- as.integer(unlist(strsplit(pruning[ACC.train.95p.max,6],",")))
ACC.train.95p.max.gen
if (DoVal == 1) {
ACC.val.95p.models <- models.to.eval[which(pruning[models.to.eval,4] >=ACC.val.all*0.95)]
pruning[ACC.val.95p.models,"pp.val"] <- 5
ACC.val.95p.max <- ACC.val.95p.models[which(pruning[ACC.val.95p.models,5]==max(pruning[ACC.val.95p.models,5]))]
ACC.val.95p.max.Npruned <- max(pruning[ACC.val.95p.models,5])
ACC.val.95p.max # this contains the number of the model in the pruning data.frame
                 # specifying the model with the maximum number of pruned descriptors
                 # and at the same time a maximum reduction of ACC.val by 5%
ACC.val.95p.max.gen <- as.integer(unlist(strsplit(pruning[ACC.val.95p.max,6],",")))
ACC.val.95p.max.gen
}

##### 1% reduction in ACC
ACC.test.99p.models <- models.to.eval[which(pruning[models.to.eval,3] >=ACC.test.all*0.99)]
pruning[ACC.test.99p.models,"pp.test"] <- 1
ACC.test.99p.max <- ACC.test.99p.models[which(pruning[ACC.test.99p.models,5]==max(pruning[ACC.test.99p.models,5]))]
ACC.test.99p.max.Npruned <- max(pruning[ACC.test.99p.models,5])
ACC.test.99p.max # this contains the number of the model in the pruning data.frame
                 # specifying the model with the maximum number of pruned descriptors
                 # and at the same time a maximum reduction of ACC.test by 1%
ACC.test.99p.max.gen <- as.integer(unlist(strsplit(pruning[ACC.test.99p.max,6],",")))
ACC.test.99p.max.gen

ACC.train.99p.models <- models.to.eval[which(pruning[models.to.eval,2] >=ACC.train.all*0.99)]
pruning[ACC.train.99p.models,"pp.train"] <- 1
ACC.train.99p.max <- ACC.train.99p.models[which(pruning[ACC.train.99p.models,5]==max(pruning[ACC.train.99p.models,5]))]
ACC.train.99p.max.Npruned <- max(pruning[ACC.train.99p.models,5])
ACC.train.99p.max # this contains the number of the model in the pruning data.frame
                 # specifying the model with the maximum number of pruned descriptors
                 # and at the same time a maximum reduction of ACC.train by 1%
ACC.train.99p.max.gen <- as.integer(unlist(strsplit(pruning[ACC.train.99p.max,6],",")))
ACC.train.99p.max.gen
if (DoVal == 1) {
ACC.val.99p.models <- models.to.eval[which(pruning[models.to.eval,4] >=ACC.val.all*0.99)]
pruning[ACC.val.99p.models,"pp.val"] <- 1
ACC.val.99p.max <- ACC.val.99p.models[which(pruning[ACC.val.99p.models,5]==max(pruning[ACC.val.99p.models,5]))]
ACC.val.99p.max.Npruned <- max(pruning[ACC.val.99p.models,5])
ACC.val.99p.max # this contains the number of the model in the pruning data.frame
                 # specifying the model with the maximum number of pruned descriptors
                 # and at the same time a maximum reduction of ACC.val by 1%
ACC.val.99p.max.gen <- as.integer(unlist(strsplit(pruning[ACC.val.99p.max,6],",")))
ACC.val.99p.max.gen
}

# creat two models
# 1. slightly pruned that ACC(test) >= 99% ACC(test)
# 2. heavily pruned that ACC(test) >= 95% ACC(test)
#----------------------------
N.pruned.models <- 0
if (NROW(ACC.test.99p.models) > 1)
   {
    N.pruned.models <- N.pruned.models +1
    print(pruning[ACC.test.99p.max,])
   } else {print("no slightly pruned model available")}
#----------------------------
# 2.
if (NROW(ACC.test.95p.models) > 1)
   {
    N.pruned.models <- N.pruned.models +1
   } else {print("no heavily pruned model available")}


# there are several possible plots here ...
# the plot matrix should be (x,y) with y=2 or y=3 if val is done
#
# --------------------------------    -------------------------------    -------------------------------
# -                              -    -                             -    -                             -
# -    plot error/ACC vs.        -    -      empty or               -    -                             -
# -    number of pruned          -    -      info on pruned         -    -                             -
# -    descriptors               -    -      descriptors            -    -                             -
# -                              -    -                             -    -                             -
# --------------------------------    -------------------------------    -------------------------------
#  original model !!!!
# --------------------------------    -------------------------------    -------------------------------
# -                              -    -                             -    -                             -
# -       training               -    -           test              -    -         val                 -
# -                              -    -                             -    -                             -
# -                              -    -                             -    -                             -
# -                              -    -                             -    -                             -
# --------------------------------    -------------------------------    -------------------------------
#  slighetly pruned model (99p)
# --------------------------------    -------------------------------    -------------------------------
# -                              -    -                             -    -                             -
# -     training                 -    -          test               -    -         val                 -
# -                              -    -                             -    -                             -
# -                              -    -                             -    -                             -
# -                              -    -                             -    -                             -
# --------------------------------    -------------------------------    -------------------------------
#  heavily pruned model (95p)
# --------------------------------    -------------------------------    -------------------------------
# -                              -    -                             -    -                             -
# -     training                 -    -          test               -    -         val                 -
# -                              -    -                             -    -                             -
# -                              -    -                             -    -                             -
# -                              -    -                             -    -                             -
# --------------------------------    -------------------------------    -------------------------------

ACC.train.min <- min(pruning[c(take.out.order[1],which(pruning[,5]>1)),2])
ACC.train.max <- max(pruning[c(take.out.order[1],which(pruning[,5]>1)),2])
ACC.test.min <- min(pruning[c(take.out.order[1],which(pruning[,5]>1)),3])
ACC.test.max <- max(pruning[c(take.out.order[1],which(pruning[,5]>1)),3])
ylims=c(min(ACC.test.min,ACC.train.min),max(ACC.test.max,ACC.train.max))
legend.pos.y=(ylims[2]-ylims[1])/2 + ylims[1]

line.w=1.5
png.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,"_pruning.png",sep="")## old file name definition: paste(func.prefix,"_popSize_",pS,"_iter_",iter,postfix,"_runs_",l,"_model_",MODEL,"_pruning.png",sep="")
if (DoVal == 1)  {## plot will be 4X3
                  png(file=png.file,height=2600,width=2000,res=180)
                  par(mfrow=c(4,3))
          } else {png(file=png.file,height=2500,width=1600,res=180)
                  par(mfrow=c(4,2))}
par(mar=c(5,5,3,5))
plot(pruning[c(1,models.to.eval),5],pruning[c(1,models.to.eval),1],type="b",lwd=line.w,xlab="Number of pruned descriptors",main="model pruning",ylab="total error",yaxt="n")
axis(2,c(-2.0,-1.5,-1.0,-0.5,0.0),as.character(c(-2.0,-1.5,-1.0,-0.5,"0.0")),las=2)
par(new=TRUE,ylbias=-0.5)
plot(pruning[c(1,models.to.eval),5],pruning[c(1,models.to.eval),2],col="green",lwd=line.w,type="b",xlab="",xaxt="n",yaxt="n",ylab="",ylim=ylims)
axis(4,c(0.0,0.2,0.4,0.6,0.8,1.0),as.character(c("0.0",0.2,0.4,0.6,0.8,1.0)),las=2)
if (DoVal == 1) {mtext("ACC(train) / ACC(test) / ACC(val)",side=4,line=2)
                 points(pruning[c(1,models.to.eval),5],pruning[c(1,models.to.eval),4],col="orange",lwd=line.w,type="b")
         } else {mtext("ACC(train) / ACC(test)",side=4,line=2)}
points(pruning[c(1,models.to.eval),5],pruning[c(1,models.to.eval),3],col="red",lwd=line.w,type="b")
points(pruning[which(pruning[,"pp.test"]==1),5],pruning[which(pruning[,"pp.test"]==1),3],pch=5,lwd=line.w)
points(pruning[which(pruning[,"pp.train"]==1),5],pruning[which(pruning[,"pp.train"]==1),2],pch=5,lwd=line.w)
points(pruning[which(pruning[,"pp.test"]==5),5],pruning[which(pruning[,"pp.test"]==5),3],pch=3,lwd=line.w)
points(pruning[which(pruning[,"pp.train"]==5),5],pruning[which(pruning[,"pp.train"]==5),2],pch=3,lwd=line.w)
#if (DoVal == 0) {legend(1,legend.pos.y,legend=c("total error","ACC(train)","ACC(test)","1% or less ACC loss","5% or less ACC loss"),pch=c(1,1,1,5,3),col=c("black","green","red","black","black"),lwd=c(2,2,2,-1,-1))
#         } else {legend(1,legend.pos.y,legend=c("total error","ACC(train)","ACC(test)","ACC(val)","1% or less ACC loss","5% or less ACC loss"),pch=c(1,1,1,1,5,3),col=c("black","green","red","orange","black","black"),lwd=c(2,2,2,2,-1,-1))}
if (PROBABILITY) {
if (DoVal == 0) {legend("bottomright",legend=c("total error","AUC(train)","AUC(test)","1% or less AUC loss","5% or less AUC loss"),pch=c(1,1,1,5,3),col=c("black","green","red","black","black"),lwd=c(2,2,2,-1,-1))
         } else {legend("bottomright",legend=c("total error","AUC(train)","AUC(test)","AUC(val)","1% or less AUC loss","5% or less AUC loss"),pch=c(1,1,1,1,5,3),col=c("black","green","red","orange","black","black"),lwd=c(2,2,2,2,-1,-1))}
} else {
if (DoVal == 0) {legend("bottomright",legend=c("total error","ACC(train)","ACC(test)","1% or less ACC loss","5% or less ACC loss"),pch=c(1,1,1,5,3),col=c("black","green","red","black","black"),lwd=c(2,2,2,-1,-1))
         } else {legend("bottomright",legend=c("total error","ACC(train)","ACC(test)","ACC(val)","1% or less ACC loss","5% or less ACC loss"),pch=c(1,1,1,1,5,3),col=c("black","green","red","orange","black","black"),lwd=c(2,2,2,2,-1,-1))}
}
#
# evaluate pruned models
plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
text(0.5,0.5,paste(ENDPOINT,"\n",SET.column.header,"\n",paste(DESCRIPTOR.postfixes,collapse=" "),"\ntestweight = ",testweight,"\nMUTATION = ",MUTATION,"\nNU = ",NU,"\nGAMMA = ",GAMMA,sep=""),cex=1.2)
if (DoVal == 1) {plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))}

T.list <- evalFunc.pruned(rep(1,ncol(dataset)-3))
plot.eval(T.list[[1]],paste("training set\noriginal model with",ncol(dataset)-3,"descriptors"))
plot.eval(T.list[[2]],paste("test set\noriginal model with",ncol(dataset)-3,"descriptors"))
if (DoVal == 1) {plot.eval(T.list[[3]],paste("validation set\noriginal model with",ncol(dataset)-3,"descriptors"))}

if (NROW(ACC.test.99p.max.gen) > 1)
   {
    T.list <- evalFunc.pruned(ACC.test.99p.max.gen)
    plot.eval(T.list[[1]],paste("training set\nusing",ncol(dataset)-3-ACC.test.99p.max.Npruned,"from",ncol(dataset)-3,"descriptors"))
    plot.eval(T.list[[2]],paste("test set\nusing",ncol(dataset)-3-ACC.test.99p.max.Npruned,"from",ncol(dataset)-3,"descriptors"))
    if (DoVal == 1) {plot.eval(T.list[[3]],paste("validation set\nusing",ncol(dataset)-3-ACC.test.99p.max.Npruned,"from",ncol(dataset)-3,"descriptors"))}
   } else {plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
           text(0.5,0.5,"- model not available -",cex=1.2)
           plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
           text(0.5,0.5,"- model not available -",cex=1.2)
           if (DoVal == 1) {plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
                            text(0.5,0.5,"- model not available -",cex=1.2)}
   }
if (NROW(ACC.test.95p.max.gen) > 1)
   {
    T.list <- evalFunc.pruned(ACC.test.95p.max.gen)
    plot.eval(T.list[[1]],paste("training set\nusing",ncol(dataset)-3-ACC.test.95p.max.Npruned,"from",ncol(dataset)-3,"descriptors"))
    plot.eval(T.list[[2]],paste("test set\nusing",ncol(dataset)-3-ACC.test.95p.max.Npruned,"from",ncol(dataset)-3,"descriptors"))
    if (DoVal == 1) {plot.eval(T.list[[3]],paste("validation set\nusing",ncol(dataset)-3-ACC.test.95p.max.Npruned,"from",ncol(dataset)-3,"descriptors"))}
   } else {plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
           text(0.5,0.5,"- model not available -",cex=1.2)
           plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
           text(0.5,0.5,"- model not available -",cex=1.2)
           if (DoVal == 1) {plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
                            text(0.5,0.5,"- model not available -",cex=1.2)}
   }
dev.off()



print("... done pruning")

################################################################################################################################
################################################################################################################################
if (ADDITIONAL.MODELS) {
print("---- building additional models and saving them as .RData file")
N.desc.used <- sum(popul[1,-1])

comps.ttv <- which(!is.na(dataset[,ncol(dataset)]) & (dataset[,2]=="train" | dataset[,2]=="test" | dataset[,2]=="val"))

################################################################################################################
################################################################################################################
##### Additional FUNCTION FOR PLOTING - Additional FUNCTION FOR PLOTING - Additional FUNCTION FOR PLOTING --####
##### Additional FUNCTION FOR PLOTING - Additional FUNCTION FOR PLOTING - Additional FUNCTION FOR PLOTING --####
################################################################################################################
### Function --- saving an arbitrary model
save.model <- function(model.GEN,model.title,model.comps)
{
print(paste("now saving model:",model.title,"---------*"))
## --- building the model
 if (PROBABILITY)
{
if (!is.na(NU)) {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[model.comps,c(which(model.GEN == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,gamma=GAMMA,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[model.comps,c(which(model.GEN == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),             error=function(x) x="ERROR")}
         } else {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[model.comps,c(which(model.GEN == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,gamma=GAMMA,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),       error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[model.comps,c(which(model.GEN == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),                   error=function(x) x="ERROR")}}
} else {
if (!is.na(NU)) {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[model.comps,c(which(model.GEN == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,gamma=GAMMA,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[model.comps,c(which(model.GEN == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")}
         } else {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[model.comps,c(which(model.GEN == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,gamma=GAMMA,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[model.comps,c(which(model.GEN == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")}}
}
selected.descriptors.names <- colnames(dataset)[which(model.GEN == 1)+2]
#------------ evaluating model
if (!is.character(model)) {
if (PROBABILITY) {
                   #--- predicting train, test and validation compounds
                   pred.train <- tryCatch(predict(model,dataset[which(dataset[,2]=="train"),which(model.GEN == 1)+2],probability=T), error=function(x) x="ERROR")
                   if (is.character(pred.train)) {print("ERROR in training set prediction during evaluation of an individual with following descriptors:")
                                                  print(colnames(dataset)[which(x == 1)+2]); return(0)}
                   TRAIN <-  data.frame(dataset[which(dataset[,2]=="train"),c(1,2,ncol(dataset))],as.factor(as.vector(pred.train)),attributes(pred.train)$probabilities[,"1"])
                   colnames(TRAIN) <- c("ID","tt.split",colnames(dataset)[ncol(dataset)],"pred","probability")
                   pred.test <- tryCatch(predict(model,dataset[which(dataset[,2]=="test"),which(model.GEN == 1)+2],probability=T), error=function(x) x="ERROR")
                   if (is.character(pred.test)) {print("ERROR in test set prediction during evaluation of an individual with following descriptors:")
                                                 print(colnames(dataset)[which(model.GEN == 1)+2]); return(0)}
                   TEST <- data.frame(dataset[which(dataset[,2]=="test"),c(1,2,ncol(dataset))],as.factor(as.vector(pred.test)),attributes(pred.test)$probabilities[,"1"])
                   colnames(TEST) <- c("ID","tt.split",colnames(dataset)[ncol(dataset)],"pred","probability")
                   if (DoVal == 1) {
                                    pred.val <- tryCatch(predict(model,dataset[which(dataset[,2]=="val"),which(model.GEN == 1)+2],probability=T), error=function(x) x="ERROR")
                                    if (is.character(pred.val)) {print("ERROR in test set prediction during evaluation of an individual with following descriptors:")
                                                                 print(colnames(dataset)[which(model.GEN == 1)+2]); return(0)}
                                    VAL <- data.frame(dataset[which(dataset[,2]=="val"),c(1,2,ncol(dataset))],as.factor(as.vector(pred.val)),attributes(pred.val)$probabilities[,"1"])
                                    colnames(VAL) <- c("ID","tt.split",colnames(dataset)[ncol(dataset)],"pred","probability")
                                    ALL <- data.frame(rbind(TRAIN,TEST,VAL))
                            } else {ALL <- data.frame(rbind(TRAIN,TEST))}
} else {
                           TRAIN <- data.frame(dataset[which(dataset[,2]=="train"),c(1:2,ncol(dataset))],"pred"=predict(model,dataset[which(dataset[,2]=="train"),c(-1,-2,-ncol(dataset))]))
                           TEST  <- data.frame(dataset[which(dataset[,2]=="test"),c(1:2,ncol(dataset))],"pred"=predict(model,dataset[which(dataset[,2]=="test"),c(-1,-2,-ncol(dataset))]))
                           if (DoVal == 0) {ALL <- data.frame(rbind(TRAIN,TEST))
                                    } else {VAL  <- data.frame(dataset[which(dataset[,2]=="val"),c(1:2,ncol(dataset))],"pred"=predict(model,dataset[which(dataset[,2]=="val"),c(-1,-2,-ncol(dataset))]))
                                            ALL <- data.frame(rbind(TRAIN,TEST,VAL))}
}  # END of if (PROBABILITY)
 #----- this is the same for PRIOBABILITY=TRUE/FALSE
                           ## constructing table with prediciton values and classes and
                           ALL <- ALL[order(ALL[,1]),]
                           ALL[,1] <- as.character(ALL[,1])
                           TT.DATA[,1] <- as.character(TT.DATA[,1])
                           FULL <- merge(TT.DATA[,1],ALL,by=1,all=TRUE)
                           #head(FULL,n=100)
                           colnames(FULL)[1] <- "ID"
                           colnames(FULL)[4] <- paste(ENDPOINT,".SVM.",SET.column.header,".",paste(DESCRIPTOR.postfixes,collapse="."),".",sum(popul[1,-1]),"vars.model_",model.title,".",MODEL,sep="")
                           if (PROBABILITY) {
                                             colnames(FULL)[5] <- paste(ENDPOINT,".SVM.",SET.column.header,".",paste(DESCRIPTOR.postfixes,collapse="."),".",sum(popul[1,-1]),"vars.model_",model.title,".",MODEL,".Positive.Prediciton.Probability",sep="")
                                            }
                           FULL[model.comps,2] <- "train"
 #-------------------------------------------------------------------------------------------------------------------------------
                           FULL.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,model.title,"_pred.txt",sep="")
                           print(FULL.file)
                           write.table(FULL,file=FULL.file,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
                           # save the domain data
                           # ... this is not necessary because its the same data that was saved for the original model!
 #-------------------------------------------------------------------------------------------------------------------------------
### ROC PLOT
ROC.plot.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,model.title,"_ROC.png",sep="")
print("saving additional ROC plots to")
print(ROC.plot.file)
png(file=ROC.plot.file,width=2000,height=2000,res=250)
par(mfrow=c(2,2))
plot(0,0,bty="n" ,xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),col="white",xaxt="n",yaxt="n")
text(0.5,0.5,paste("models evaluated:\n",model.title,"\n",FILENAME,"\n run ",l,
                   "\n\n",ENDPOINT,"\n",paste(DESCRIPTOR.postfixes,collapse=" "),"\nNU = ",NU,"\nGAMMA = ",GAMMA,"\n",ncol(dataset)-3," descriptors","\nclasses:",paste(names(class.weight),collapse="/"),
                   "   weights:",paste(round(class.weight,digits=2),collapse="/"),sep=""))
# TRAIN
roc.train <- plot.roc(TRAIN[,3],TRAIN[,5])
title(main="training set\n")
thr <- ci(roc.train, of="thresholds", thresholds="best")
thr.numeric.train <- round(as.numeric(rownames(thr$sensitivity)),digits=4)
text(0.5,0.2,paste("AUC =",round(roc.train$auc,digits=3),"\nbest threshold =",round(thr.numeric.train,digits=1)))
plot(thr)
# TEST
roc.test <- plot.roc(TEST[,3],TEST[,5])
thr <- ci(roc.test, of="thresholds", thresholds="best")
thr.numeric.test <- round(as.numeric(rownames(thr$sensitivity)),digits=4)
text(0.5,0.2,paste("AUC =",round(roc.test$auc,digits=3),"\nbest threshold =",round(thr.numeric.test,digits=1)))
plot(thr)
title(main="test set\n")
if (DoVal == 1)
{
## VAL
roc.val <- plot.roc(VAL[,3],VAL[,5])
text(0.5,0.2,paste("AUC =",round(roc.val$auc,digits=3)))
title(main="validation set set\n")
thr <- ci(roc.val, of="thresholds", thresholds="best")
thr.numeric.val <- round(as.numeric(rownames(thr$sensitivity)),digits=4)

}
dev.off()

#-----------------------------------------------------------------------------------------------------------------
if (FALSE) {
## SE/SP vs. thresholds PLOT
ROC.thr.plot.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,model.title,"_ROC_thresholds.png",sep="")
print("saving plots threshold vs SE/SP to")
print(ROC.thr.plot.file)
png(file=ROC.thr.plot.file,width=2000,height=2000,res=250)
par(mfrow=c(2,2))
plot(0,0,bty="n" ,xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),col="white",xaxt="n",yaxt="n")
text(0.5,0.5,paste("models evaluated:\n",FILENAME,"\n run ",l,
                   "\n\n",ENDPOINT,"\n",paste(DESCRIPTOR.postfixes,collapse=" "),"\nNU = ",NU,"\nGAMMA = ",GAMMA,"\n",ncol(dataset)-3," descriptors","\nclasses:",paste(names(class.weight),collapse="/"),
                   "   weights:",paste(round(class.weight,digits=2),collapse="/"),sep=""))
## TRAIN
plot(roc.train$thresholds,roc.train$sensitivities,lty=1,type="l",pch=20,col="red",lwd=1.2,ylab="Sensitivity/Speceficity",xlab="threshold")
points(roc.train$thresholds,roc.train$specificities,lty=1,type="l",pch=20,col="blue",lwd=1.2)
SE.SP.dif.min <- order(abs(roc.train$sensitivities - roc.train$specificities))[1:2]
SE.SP.eq.thr <- mean(roc.train$thresholds[SE.SP.dif.min])
SE.SP.eq.thr.train <- SE.SP.eq.thr
SE.SP.eq     <- mean(roc.train$sensitivities[SE.SP.dif.min])
abline(v=SE.SP.eq.thr,lty=3); points(SE.SP.eq.thr,SE.SP.eq,pch=8)
legend("topleft",legend=c("Sensitivity","Speceficity"),col=c("red","blue"),lty=1,bg="white",cex=0.8)
text(SE.SP.eq.thr,0.1,paste("SE = SP =",round(SE.SP.eq,digits=3),"\nat threshold =",round(SE.SP.eq.thr,digits=3)))
title(main="training set")
## TEST
plot(roc.test$thresholds,roc.test$sensitivities,lty=1,type="l",pch=20,col="red",lwd=1.2,ylab="Sensitivity/Speceficity",xlab="threshold")
points(roc.test$thresholds,roc.test$specificities,lty=1,type="l",pch=20,col="blue",lwd=1.2)
SE.SP.dif.min <- order(abs(roc.test$sensitivities - roc.test$specificities))[1:2]
SE.SP.eq.thr <- mean(roc.test$thresholds[SE.SP.dif.min])
SE.SP.eq.thr.test <- SE.SP.eq.thr
SE.SP.eq     <- mean(roc.test$sensitivities[SE.SP.dif.min])
abline(v=SE.SP.eq.thr,lty=3); points(SE.SP.eq.thr,SE.SP.eq,pch=8)
legend("topleft",legend=c("Sensitivity","Speceficity"),col=c("red","blue"),lty=1,bg="white",cex=0.8)
text(SE.SP.eq.thr,0.1,paste("SE = SP =",round(SE.SP.eq,digits=3),"\nat threshold =",round(SE.SP.eq.thr,digits=3)))
title(main="test set")
## VAL
plot(roc.val$thresholds,roc.val$sensitivities,lty=1,type="l",pch=20,col="red",lwd=1.2,ylab="Sensitivity/Speceficity",xlab="threshold")
points(roc.val$thresholds,roc.val$specificities,lty=1,type="l",pch=20,col="blue",lwd=1.2)
SE.SP.dif.min <- order(abs(roc.val$sensitivities - roc.val$specificities))[1:2]
SE.SP.eq.thr <- mean(roc.val$thresholds[SE.SP.dif.min])
SE.SP.eq.thr.val <- SE.SP.eq.thr
SE.SP.eq     <- mean(roc.val$sensitivities[SE.SP.dif.min])
abline(v=SE.SP.eq.thr,lty=3); points(SE.SP.eq.thr,SE.SP.eq,pch=8)
legend("topleft",legend=c("Sensitivity","Speceficity"),col=c("red","blue"),lty=1,bg="white",cex=0.8)
text(SE.SP.eq.thr,0.1,paste("SE = SP =",round(SE.SP.eq,digits=3),"\nat threshold =",round(SE.SP.eq.thr,digits=3)))
title(main="validation set set")
dev.off()
}
#-----------------------------------------------------------------------------------------------------------------
# DENSITY PLOT
empty <- data.frame()
empty.plot <- ggplot (empty) + geom_blank() + xlim(0, 10) + ylim(0, 10) +
                annotate("text", x = 5, y = 5, label = paste("models evaluated:\n",FILENAME,"\n run ",l,
                   "\n\n",ENDPOINT,"\n",paste(DESCRIPTOR.postfixes,collapse=" "),"\nNU = ",NU,"\nGAMMA = ",GAMMA,"\n",ncol(dataset)-3," descriptors\nclasses:",paste(names(class.weight),collapse="/"),
                   "   weights:",paste(round(class.weight,digits=2),collapse="/"),sep=""), size=3.5)   +
                theme(panel.background = theme_rect(fill='white', colour='black'))

## TRAIN
dens <- na.omit(data.frame(cbind("class"=TRAIN[,3],"pred"=TRAIN[,5])))
dens[,1] <- as.factor(as.character(dens[,1]))
train.1 <- ggplot(dens,aes(x=pred,fill=class)) + geom_density(alpha=0.3) +  xlab("Poritive Prediciton Probability") +  xlim(min(dens[,2]),max(dens[,2])) +
           ggtitle("training set")
## TEST
dens <- na.omit(data.frame(cbind("class"=TEST[,3],"pred"=TEST[,5])))
dens[,1] <- as.factor(as.character(dens[,1]))
test.1 <- ggplot(dens,aes(x=pred,fill=class)) + geom_density(alpha=0.3) +  xlab("Poritive Prediciton Probability") +  xlim(min(dens[,2]),max(dens[,2])) +
          ggtitle("test set")
## VAL
dens <- na.omit(data.frame(cbind("class"=VAL[,3],"pred"=VAL[,5])))
dens[,1] <- as.factor(as.character(dens[,1]))
val.1 <- ggplot(dens,aes(x=pred,fill=class)) + geom_density(alpha=0.3) +  xlab("Poritive Prediciton Probability") +  xlim(min(dens[,2]),max(dens[,2])) +
          ggtitle("validation set")
##
dens.plot.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,model.title,"_ROC_density.png",sep="")
print("saving density plot to")
print(dens.plot.file)
png(file=dens.plot.file,width=2500,height=2000,res=250)
grid.arrange(arrangeGrob(empty.plot,train.1,test.1,val.1,heights=c(1/2,1/2)),main="denisity plot")
dev.off()

##------------- now adding the ROC data to the model object
print("Adding ROC data to the model")
auc.model <- list("train"=as.numeric(roc.train$auc),"test"=as.numeric(roc.test$auc),"val"=as.numeric(roc.val$auc))
thr.boot.model <- list("train"=thr.numeric.train,"test"=thr.numeric.test,"val"=thr.numeric.val)
thr.SE.SP.eq <- list("train"=SE.SP.eq.thr.train,"test"=SE.SP.eq.thr.test,"val"=SE.SP.eq.thr.val)
model[["ROC"]] <- list("auc"=auc.model,"thr.bootstrap"=thr.boot.model,"thr.SE.SP.eq"=thr.SE.SP.eq)
model[["DESCRIPTORS"]] <- selected.descriptors.names

##------------- all plotting done...
                           model.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,model.title,".RData",sep="") ## old file name definition: paste(func.prefix,"_popSize_",pS,"_iter_",iter,postfix,"_runs_",l,"_model_",MODEL,model.file.name.postfix,".RData",sep="")
                           save(model,file=model.file)
                           stats.out <-  model.stats(model,dataset)
                           stats.out$N.desc <- sum(GEN)
                           stats.file <- paste(FILENAME,"_runs_",l,"_model_",MODEL,model.title,"_stats.txt",sep="") ## old file name definition: paste(func.prefix,"_popSize_",pS,"_iter_",iter,postfix,"_runs_",l,"_model_",MODEL,model.file.name.postfix,"_stats.txt",sep="")
                           write.table(stats.out,file=stats.file,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
                           # save the prediction
                   } else {
                           print("building model failed!")
                          }
print("Done saving model -------------------------------**")
}
################################################################################################################
################################################################################################################



#-------------------------------------------
## full model with all available compounds
## _all_ttv
print(      "1. unpruned models:")
print(paste(" 1.1. a model built with all available compounds (training + test (+ validation) set) and all optimized",ncol(dataset)-3,"descriptors"))
print(paste("      --> all_ttv_",ncol(dataset)-3,sep=""))
GEN <- rep(1,ncol(dataset)-3)
save.model(GEN,"_all_ttv",comps.ttv)
## _tt
print(paste(" 1.2. a model built with all atrain and test compounds (validation set remains for validation) and all optimized",ncol(dataset)-3,"descriptors"))
print(paste("      --> _tt_",ncol(dataset)-3,sep=""))
GEN <- rep(1,ncol(dataset)-3)
save.model(GEN,"_tt",c(comps.train,comps.test))


###### --- save pruned models with ACC(test) maximum 1% less than in the full model
if (NROW(ACC.test.99p.max.gen) > 1)
{
print(      "2. slightly pruned models:")
print(paste(" 2.1. a model built with all training set compounds and ",ncol(dataset)-3-ACC.test.99p.max.Npruned,"descriptors"))
print(paste("      --> pruned_",ncol(dataset)-3-ACC.test.99p.max.Npruned,sep=""))
### slightly pruned
### _pruned_XX
GEN <- ACC.test.99p.max.gen
model.file.name.postfix <- paste("_pruned_",sum(GEN),sep="")
save.model(GEN,model.file.name.postfix,comps.train)

# _tt_pruned_XX
print(paste(" 2.2. a model built with all train and test compounds (validation remains for validation!) and ",ncol(dataset)-3-ACC.test.99p.max.Npruned,"descriptors"))
print(paste("      --> _tt_pruned_",ncol(dataset)-3-ACC.test.99p.max.Npruned,sep=""))
model.file.name.postfix <- paste("_tt_pruned_",sum(GEN),sep="")
model.file.name.postfix
save.model(GEN,model.file.name.postfix,c(comps.train,comps.test))

# _all_ttv_pruned_XX
print(paste(" 2.3. a model built with all available compounds (training + test (+ validation) set) and ",ncol(dataset)-3-ACC.test.99p.max.Npruned,"descriptors"))
print(paste("      --> all_ttv_pruned_",ncol(dataset)-3-ACC.test.99p.max.Npruned,sep=""))
model.file.name.postfix <- paste("_all_ttv_pruned_",sum(GEN),sep="")
model.file.name.postfix
save.model(GEN,model.file.name.postfix,comps.ttv)
} else {print("no slightly pruned model to save ... sorry!")}

if (NROW(ACC.test.95p.max.gen) > 1)
{
## _pruned_YY
print(      "3. heavily pruned models:")
print(paste(" 3.1. a model built with all training set compounds and ",ncol(dataset)-3-ACC.test.95p.max.Npruned,"descriptors"))
print(paste("      --> pruned_",ncol(dataset)-3-ACC.test.95p.max.Npruned,sep=""))
GEN <- ACC.test.95p.max.gen
model.file.name.postfix <- paste("_pruned_",sum(GEN),sep="")
save.model(GEN,model.file.name.postfix,comps.train)

## _tt_pruned_YY
print(paste(" 3.2. a model built with all train and test compounds (validation set remains for validation) and ",ncol(dataset)-3-ACC.test.95p.max.Npruned,"descriptors"))
print(paste("      --> _tt_pruned_",ncol(dataset)-3-ACC.test.95p.max.Npruned,sep=""))
model.file.name.postfix <- paste("_tt_pruned_",sum(GEN),sep="")
save.model(GEN,model.file.name.postfix,c(comps.train,comps.test))

# _all_ttv_pruned_YY
print(paste(" 3.3. a model built with all available compounds (training + test (+ validation) set) and ",ncol(dataset)-3-ACC.test.95p.max.Npruned,"descriptors"))
print(paste("      --> all_ttv_pruned_",ncol(dataset)-3-ACC.test.95p.max.Npruned,sep=""))
model.file.name.postfix <- paste("_all_ttv_pruned_",sum(GEN),sep="")
save.model(GEN,model.file.name.postfix,comps.ttv)

} else {print("no heavily pruned model to save ... sorry!")}


if (all(!is.na(output.PRUNED)))
                            {
                            print(paste("You requested the output of models",paste(output.PRUNED,collapse=", "),"from the pruned model list"))
                            for (i in 1:NROW(output.PRUNED))
                            {
                            if (!output.PRUNED[i]==0&!output.PRUNED[i]>max(pruning[,5]))
                            {
                            print(paste("Preparing output for pruned model",output.PRUNED[i],", i.e. the model where",output.PRUNED[i],"descriptors were ignored"))
                            if (output.PRUNED[i] == 1) {pp <- take.out.order[1]
                                                } else {pp <- which(pruning[,5]==output.PRUNED[i])}
                            optGEN <- as.integer(unlist(strsplit(pruning[pp,6],",")))
                            print(paste("with",sum(optGEN),"descriptors and the corresponding gen:"))
                            print(optGEN)
                            model.file.name.postfix <- paste("_pruned_",sum(optGEN),sep="")
                            save.model(optGEN,model.file.name.postfix,comps.train)
                            model.file.name.postfix <- paste("_tt_pruned_",sum(optGEN),sep="")
                            save.model(optGEN,model.file.name.postfix,c(comps.train,comps.test))
                            model.file.name.postfix <- paste("_all_ttv_pruned_",sum(optGEN),sep="")
                            save.model(optGEN,model.file.name.postfix,comps.ttv)
                            } else {print(paste("The model",output.PRUNED[i]," you requested does not exist!"))}
                            }
                            }

} else { if (!is.na(output.PRUNED)) {cat("-----WARNING---------------------------------------------------\n")
                                     cat("-----------  You specified output.PRUNED to output additional \n")
                                     cat("-----------  pruned models but ADDITIONAL.MODELS is FALSE. No pruned\n-----------  outout is produced\n")
                                     }
         cat(" to save pruned models and related domain data redo the evaluation with ADDITIONAL.MODELS = TRUE ")}

} #END of if (PRUNING)
