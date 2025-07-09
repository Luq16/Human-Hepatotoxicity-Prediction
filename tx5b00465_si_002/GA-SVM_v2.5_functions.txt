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




#############################################################################################################
# plotting a confusion matrix of any size mxn (tested up to 6x6)
plot.eval <- function(Tab,plot.title="") {
  NR <- nrow(Tab)
  NC <- ncol(Tab)
  # adjust size of Text to not overlap so much
  Text.size <- 1.0
  if (NC > 3 | NR > 3) {Text.size <- 0.8}
  if (NC > 5 | NR > 5) {Text.size <- 0.7}
#
  plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,0),ylim=c(0,0.5),main=plot.title,cex=Text.size,bty="n");
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
###################################################################################################################
plot.GA.output <- function(GAmodel)
{
## write population to file!!!
## this is task one and
#  save(GAmodel,file="GAmodel_for_test.RData")
#  cat("This is plot.GA.output\n")
  popul <- data.frame(cbind("eval"=GAmodel$evaluations,as.data.frame(GAmodel$population)))
  colnames(popul) <- c("eval",colnames(dataset)[3:(ncol(dataset)-1)])
  popul.file <- paste(FILENAME,"_runs_",ii,"_population.txt",sep="")
  write.table(popul,file=popul.file,row.names=FALSE)
  
## evaluate population, count similar individuals etc.
if (POP.EVAL) {
cat("\n\n## evaluation of the current population ##\n")
unique.pop <- ddply(popul,colnames(popul)[-1],nrow)
cat(" ",nrow(unique.pop),"unique chromosoms are in the population of",nrow(popul),"chromosoms\n")

n.chosen.descs <- apply(popul[,-1],1, function(x) sum(x))
desc.min=min(n.chosen.descs)
desc.max=max(n.chosen.descs)
cat("  Between",desc.min,"and",desc.max,"descriptors have bee selected\n")
cat("  Fitness ranges from ",min(popul[,1]),"to",max(popul[,1]),"with a mean of",mean(popul[,1]),"and a standard deviation of",sd(popul[,1]),"\n")

POP.matrix <- as.matrix(popul)
cat("top chromosoms:\n")
similar.descs <- data.frame("n.identical.descs"=0,"number"=0,"n.identical.chroms"=0)


unique.pop <- as.matrix(unique.pop)
for (i in 1:nrow(unique.pop))
{
 XX <- which(apply(POP.matrix[,-1],1,function(x) identical(x,unique.pop[i,-ncol(unique.pop)])))[1]
 if (!exists("unique.pop.fit")) {unique.pop.fit <- POP.matrix[XX,1]
                         } else {unique.pop.fit <- c(unique.pop.fit,POP.matrix[XX,1])}
}

for (i in 1:nrow(unique.pop))
{
 best <- order(unique.pop.fit)[i]
 identical.rows <- apply(POP.matrix[,-1],1,function(x) identical(x,unique.pop[best,-ncol(unique.pop)]))
 if (i==1) {best.of.the.best = best}
 # now match the selected descriptors
 match.with.best <- which(which(POP.matrix[best,-1]==1) %in% which(POP.matrix[best.of.the.best,-1]==1))
 if (length(match.with.best)>0) {n.match.with.best <- NROW(match.with.best  )}
 if (i < 20) {
     cat("  ",i,") fitness=",round(unique.pop.fit[best],digits=3)," chromosom makes up ",round((NROW(which(identical.rows))/nrow(POP.matrix))*100,digits=1),"% of population",sep="")
     cat(" and ",n.match.with.best," descriptors (",round((n.match.with.best/NROW(which(POP.matrix[best,-1]==1)))*100,digits=0),"%) are also selected in best chromosom\n",sep="")
 }
 similar.descs[i,"n.identical.descs"] <- n.match.with.best
 similar.descs[i,"number"] <- best
 similar.descs[i,"n.identical.chroms"] <- NROW(which(identical.rows))

# if (i>50) {break}
}

B <- sum(POP.matrix[best.of.the.best,-1])
cat(round((NROW(which(similar.descs[,1]/B>0.9))/nrow(POP.matrix))*100,digits=1),"% of chromosoms share > 90% of descriptors with the best chromosom\n",sep="")
cat(round((NROW(which(similar.descs[,1]/B>0.8))/nrow(POP.matrix))*100,digits=1),"% of chromosoms share > 80% of descriptors with the best chromosom\n",sep="")
cat(round((NROW(which(similar.descs[,1]/B>0.7))/nrow(POP.matrix))*100,digits=1),"% of chromosoms share > 70% of descriptors with the best chromosom\n",sep="")
cat("## --- end of evaluation of population ##\n")
} else {cat("If you would like an evaluation of the population at this point rerun with option POP.EVAL=T\n")}


original.par <- par()
if (DEBUG) {print("DEBUG: function plot.GA.output -- starting plot")}
png.file <- paste(FILENAME,"_runs_",ii,".png",sep="")
if (DoVal == 0) {
                 if (PROBABILITY) {
                                   png(file=png.file,height=3000,width=2000,res=300)
                                   par(mfrow=c(3,2))
                           } else {
                                   png(file=png.file,height=2000,width=2000,res=300)
                                   par(mfrow=c(2,2))
                                  }
         } else {
                 if (PROBABILITY) {
                                   png(file=png.file,height=2800,width=2200,res=300)
                                   par(mfrow=c(4,3))
                           } else {
                                   png(file=png.file,height=1600,width=2200,res=300)
                                   par(mfrow=c(2,3))
                                  }
                }
model.name <- paste(ENDPOINT,"_",SET.column.header,sep="")
plot(1:nrow(OUT.bm),OUT.bm[,2],type="l",xlab="iteration step",ylab="fitness",ylim=c(min(OUT.bm),max(OUT.bm)),main=model.name);
points(1:nrow(OUT.bm),OUT.bm[,1],type="l",col="blue");
legend("topright",col=c("black","blue"),legend=c("mean","best"),lty=1)
par(mar=c(0,0,0,0))
plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),bty="n");
if (is.na(MUTATION)) {MUTATIONtext="default"} else {MUTATIONtext=paste(round(MUTATION,digits=4))}
if (is.na(NU)) {NU.text="was not defined"} else {NU.text=paste("=",NU)}
if (is.na(GAMMA)) {GAMMA.text="was not defined"} else {GAMMA.text=paste("=",GAMMA)}
if (CROSSVAL) {crossval.text=paste("\n",Nfold,"-fold cross val.",sep="")} else {crossval.text=""}
text(0.5,0.5,paste(model.name,"\n\nSVM with ",kernel.type," kernel\n",svm.type,"\n nu",NU.text,"\n gamma",GAMMA.text,"\niterations =", iter,"\npopSize = ",pS,
                   "\n variables =",sum(GAmodel$population[order(GAmodel$evaluations)[1],]),"/",ncol(GAmodel$population),
                   "\neval. property:\n",evalFunc.name,"\n with N(desc)",N.chrom.corr,"\n mutation chance",MUTATIONtext,"\nclasses:",paste(names(class.weight),collapse="/"),
                   "   weights:",paste(round(class.weight,digits=2),collapse="/"),crossval.text))

## <- select best model
best.chrom <- GAmodel$population[order(GAmodel$evaluations)[1],]

#cat("\nuse descriptors :----------------------------\n")
#print(colnames(dataset)[which(best.chrom == 1)+2])

# --------------- now prediciton and ploting
if (length(grep("regression",svm.type))==0) {
if (!PROBABILITY) {
if (CROSSVAL) {TT.comps <- which(dataset[,2]=="train"|dataset[,2]=="test")} else {TT.comps <- which(dataset[,2]=="train")}
if (!is.na(NU)) {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,c(which(best.chrom == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,gamma=GAMMA,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,c(which(best.chrom == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),             error=function(x) x="ERROR")}
         } else {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,c(which(best.chrom == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,gamma=GAMMA,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),       error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,c(which(best.chrom == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),                   error=function(x) x="ERROR")}}
if (!is.character(model)) { #do this if is.character(model)==FALSE and model was built succesfully
## Train
  if (DEBUG) {print("DEBUG: function plot.GA.output -- model was succesfully built")}
  Ttrain <- table(dataset[TT.comps,ncol(dataset)],predict(model,dataset[TT.comps,which(best.chrom == 1)+2]))
  Diagonal.sum <- 0; for (d.sum in 1:ncol(Ttrain)) {Diagonal.sum <- Diagonal.sum + Ttrain[d.sum, d.sum]}
  ACC.train <- Diagonal.sum/sum(Ttrain)
## Test
  Ttest <- table(dataset[which(dataset[,2]=="test"),ncol(dataset)],predict(model,dataset[which(dataset[,2]=="test"),which(best.chrom == 1)+2]))
  Diagonal.sum <- 0; for (d.sum in 1:ncol(Ttest)) {Diagonal.sum <- Diagonal.sum + Ttest[d.sum, d.sum]}
  ACC.test <- Diagonal.sum/sum(Ttest)
## evaluating the validation set
  if (DoVal == 1) {
     Tval <- table(dataset[which(dataset[,2]=="val"),ncol(dataset)],predict(model,dataset[which(dataset[,2]=="val"),which(best.chrom == 1)+2]))
     Diagonal.sum <- 0; for (d.sum in 1:ncol(Tval)) {Diagonal.sum <- Diagonal.sum + Tval[d.sum, d.sum]}
     ACC.val <- Diagonal.sum/sum(Tval)
    }

  if (DoVal == 1) {plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),bty="n")
                   text(0.5,0.5,paste("tolerance =",tolerance,"\n epsilon (insensitivity-loss) =",epsilon,"\nN(support vectors) =",NROW(model$index)))
                  }
# par(original.par)
## plotting the confusion matrices
  plot.eval(Ttrain,"training set")
  plot.eval(Ttest,"test set")
  if (DoVal == 1) {plot.eval(Tval,"validation set")}
## end of png plot!
} else { # this is done when is.character(model)==TRUE
  if (DEBUG) {print("DEBUG: function plot.GA.output -- model could not be built")}
  #
  plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),main="training set")
  text(0.5,0.5,"ERROR")
  plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),main="test set")
  text(0.5,0.5,"ERROR")
  if (DoVal == 1) {plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),main="validation set")
                   text(0.5,0.5,"ERROR")}
}# END of if (!is.character(model)) {
}else {# ELSE of if (!PROBABILITY)     <----------------  PROBABILITY USING AUC!!!
#-----------------------------------------------------------------------------------------------

if (CROSSVAL) {TT.comps <- which(dataset[,2]=="train"|dataset[,2]=="test")} else {TT.comps <- which(dataset[,2]=="train")}
if (!is.na(NU)) {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,c(which(best.chrom == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,gamma=GAMMA, probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,c(which(best.chrom == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU, probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),             error=function(x) x="ERROR")}
         } else {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,c(which(best.chrom == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,gamma=GAMMA, probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),       error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[TT.comps,c(which(best.chrom == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type, probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),                   error=function(x) x="ERROR")}}
if (DoVal == 1){plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),bty="n")
       text(0.5,0.5,paste("tolerance =",tolerance,"\n epsilon (insensitivity-loss) =",epsilon,"\nN(support vectors) =",NROW(model$index)))}
if (!is.character(model)) { #do this if is.character(model)==FALSE and model was built succesfully
### ROC curves!!!
    pred.train <- tryCatch(predict(model,dataset[TT.comps,which(best.chrom == 1)+2],probability=T), error=function(x) x="ERROR")
    roc.train <- plot.roc(dataset[TT.comps,ncol(dataset)],attributes(pred.train)$probabilities[,"1"])
    thr <- ci(roc.train, of="thresholds", thresholds="best")
    thr.numeric.train <- round(as.numeric(rownames(thr$sensitivity)),digits=4)
    text(0.4,0.2,paste("AUC =",round(roc.train$auc,digits=3),"\n best thr. =",round(thr.numeric.train,digits=1)))
    plot(thr)
    title(main="trainings set\n")
    if (!CROSSVAL) {
                    Test.comps=which(dataset[,2]=="test")
                    pred.test <- tryCatch(predict(model,dataset[Test.comps,which(best.chrom == 1)+2],probability=T), error=function(x) x="ERROR")
            } else {
                    Test.comps=which(dataset[,2]=="train"|dataset[,2]=="test")
                    pred.test <- do.SVM.crossval(Test.comps,fold=Nfold,selected.descs=which(best.chrom == 1)+2)
                   }
    roc.test <- plot.roc(dataset[Test.comps,ncol(dataset)],attributes(pred.test)$probabilities[,"1"])
    if (CROSSVAL) {title(main=paste(Nfold,"-fold cross validation\n",sep=""))} else {title(main="test set\n")}
    thr <- ci(roc.test, of="thresholds", thresholds="best")
    thr.numeric.test <- round(as.numeric(rownames(thr$sensitivity)),digits=4)
    text(0.4,0.2,paste("AUC =",round(roc.test$auc,digits=3),"\n best thr. =",round(thr.numeric.test,digits=1)))
    plot(thr)
if (DoVal == 1)
 {
    pred.val <- tryCatch(predict(model,dataset[which(dataset[,2]=="val"),which(best.chrom == 1)+2],probability=T), error=function(x) x="ERROR")
    roc.val <- plot.roc(dataset[which(dataset[,2]=="val"),ncol(dataset)],attributes(pred.val)$probabilities[,"1"])
    text(0.4,0.2,paste("AUC =",round(roc.val$auc,digits=3)))
    title(main="validation set\n")
  }
#### now plot confusion matrix using threshold mean(train;test)
thr.use <- mean(thr.numeric.train,thr.numeric.test)
plot.eval(table(dataset[TT.comps,ncol(dataset)],as.factor(sapply(attributes(pred.train)$probabilities[,"1"],function(x) if(x>=thr.use){1}else{0}))),
          paste("\n\ntraining set\nthr. = mean(train,test) =",round(thr.use,digits=3)))
if (CROSSVAL) {test.plot.title="\n\ncross validation"} else {test.plot.title=paste("\n\ntest set\nthr. = mean(train,test) =",round(thr.use,digits=3))}
plot.eval(table(dataset[Test.comps,ncol(dataset)],as.factor(sapply(attributes(pred.test)$probabilities[,"1"],function(x) if(x>=thr.use){1}else{0}))),
          test.plot.title)
if (DoVal == 1) {plot.eval(table(dataset[which(dataset[,2]=="val"),ncol(dataset)],as.factor(sapply(attributes(pred.val)$probabilities[,"1"],function(x) if(x>=thr.use){1}else{0}))),
                           paste("\n\nvalidation set\nthr. = mean(train,test) =",round(thr.use,digits=3)))}
### and density plots....
dens <- na.omit(data.frame(cbind("class"=dataset[TT.comps,ncol(dataset)],"pred"=attributes(pred.train)$probabilities[,"1"])))
#print(colnames(dens))
dens[,1] <- as.factor(dens[,1])
vp.train <- viewport(x=0.2,y=0.13,width=0.39,height=0.27,clip="on")
gg.train <- ggplot(dens,aes(x=pred,fill=class)) + geom_density(alpha=0.3) +  xlab("Pos. Pred. Probability") +
            xlim(min(dens[,2]),max(dens[,2])) +  theme(axis.title.x=element_text(size=10),axis.title.y=element_text(size=10))+
            ggtitle("training set")

plot(gg.train, vp=vp.train)

dens <- na.omit(data.frame(cbind("class"=dataset[Test.comps,ncol(dataset)],"pred"=attributes(pred.test)$probabilities[,"1"])))
dens[,1] <- as.factor(dens[,1])
if (CROSSVAL) {test.plot.title=paste(Nfold,"-fold cross val.",sep="")} else {test.plot.title="test set"}
vp.test <- viewport(x=0.52,y=0.13,width=0.39,height=.27  ,clip="on")
gg.test <- ggplot(dens,aes(x=pred,fill=class)) + geom_density(alpha=0.3) +  xlab("Pos. Pred. Probability") +
            xlim(min(dens[,2]),max(dens[,2])) + theme(axis.title.x=element_text(size=10),axis.title.y=element_text(size=10))+
            ggtitle(test.plot.title)
plot(gg.test, vp=vp.test)
plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),bty="n")

if (DoVal == 1)
{
dens <- na.omit(data.frame(cbind("class"=dataset[which(dataset[,2]=="val"),ncol(dataset)],"pred"=attributes(pred.val)$probabilities[,"1"])))
dens[,1] <- as.factor(dens[,1])
vp.val <- viewport(x=0.87,y=0.13,width=.39,height=.27,clip="on")
gg.val <- ggplot(dens,aes(x=pred,fill=class)) + geom_density(alpha=0.3) +  xlab("Pos. Pred. Probability") +
            xlim(min(dens[,2]),max(dens[,2])) + theme(axis.title.x=element_text(size=10),axis.title.y=element_text(size=10))+
            ggtitle("validation set")
plot(gg.val, vp=vp.val)
plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),bty="n")
}

} # END of if (!is.character(model)) {
} # END of if (!PROBABILITY)
} else {# ELSE of if (length(grep("regression",svm.type))==0)
if (!is.na(NU)) {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(best.chrom == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,gamma=GAMMA,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(best.chrom == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,tolerance=tolerance,epsilon=epsilon),             error=function(x) x="ERROR")}
         } else {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(best.chrom == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,gamma=GAMMA,tolerance=tolerance,epsilon=epsilon),       error=function(x) x="ERROR")
                            } else {model <- tryCatch(svm(as.formula(form),data=dataset[which(dataset[,2]=="train"),c(which(best.chrom == 1)+2,ncol(dataset))], type=svm.type,kernel=kernel.type,tolerance=tolerance,epsilon=epsilon),                   error=function(x) x="ERROR")}}

if (DoVal == 1){plot(0,0,col=0,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),bty="n")
       text(0.5,0.5,paste("tolerance =",tolerance,"\n epsilon (insensitivity-loss) =",epsilon,"\nN(support vectors) =",NROW(model$index)))}

if (!is.character(model)) { #do this if is.character(model)==FALSE and model was built succesfully

    par(mar=c(5.1,4.1,4.1,2.1))
    pred.train <- tryCatch(predict(model,dataset[which(dataset[,2]=="train"),which(best.chrom == 1)+2]), error=function(x) x="ERROR")
# statistics .... man, that's the real stuff!!
         RMS.train=sqrt(sum((dataset[which(dataset[,2]=="train"),ncol(dataset)]-pred.train)^2)/NROW(pred.train))
         RMS.train
         pred.diff <- pred.train-mean(pred.train)
         true.diff <- dataset[which(dataset[,2]=="train"),ncol(dataset)]- mean(dataset[which(dataset[,2]=="train"),ncol(dataset)])
         Rsquare.train=(sum(pred.diff*true.diff)/sqrt(sum(pred.diff^2)*sum(true.diff^2)))^2
         Rsquare.train
    plot(pred.train,dataset[which(dataset[,2]=="train"),ncol(dataset)],ylab=ENDPOINT)
    title(main=paste("training set\nR^2 =",round(Rsquare.train,digits=2),"; RMS =",round(RMS.train,digits=2)))
    abline(a=0,b=1,lty=2)
    pred.test <- tryCatch(predict(model,dataset[which(dataset[,2]=="test"),which(best.chrom == 1)+2],probability=T), error=function(x) x="ERROR")
# statistics .... oh, so beautiful!!
         RMS.test=sqrt(sum((dataset[which(dataset[,2]=="test"),ncol(dataset)]-pred.test)^2)/NROW(pred.test))
         RMS.test
         pred.diff <- pred.test-mean(pred.test)
         true.diff <- dataset[which(dataset[,2]=="test"),ncol(dataset)]- mean(dataset[which(dataset[,2]=="test"),ncol(dataset)])
         Rsquare.test=(sum(pred.diff*true.diff)/sqrt(sum(pred.diff^2)*sum(true.diff^2)))^2
         Rsquare.test
    plot(pred.test,dataset[which(dataset[,2]=="test"),ncol(dataset)],ylab=ENDPOINT)
    title(main=paste("test set\nR^2 =",round(Rsquare.test,digits=2),"; RMS =",round(RMS.test,digits=2)))
    abline(a=0,b=1,lty=2)

###     for parameter space scan ourput the model statistics in a form identical to ROC_stats
stats.Rsquare <- data.frame("ENDP"=ENDPOINT,"TT.split"=SET.column.header,"desc"=paste(DESCRIPTOR.postfixes,collapse=" "),
                        "status"=NA,"popSize"=pS,"iter"=iter,"runs"=ii,
                        "nu"=NU,"gamma"=GAMMA,"tolerance"=tolerance,"epsilon"=epsilon,
                        "testweight"=testweight,"mutationChance"=MUTATION,
                        "N.desc"=sum(popul[1,-1]),"N.sv"=NROW(model$index),
                        "Rsquare.train"=round(Rsquare.train,digits=4),"Rsquare.test"=round(Rsquare.test,digits=4),
                        "Rsquare.val"=NA)
if (DoVal == 1)
 {
    pred.val <- tryCatch(predict(model,dataset[which(dataset[,2]=="val"),which(best.chrom == 1)+2],probability=T), error=function(x) x="ERROR")
# statistics .... and again, just incredible!!
         RMS.val=sqrt(sum((dataset[which(dataset[,2]=="val"),ncol(dataset)]-pred.val)^2)/NROW(pred.val))
         RMS.val
         pred.diff <- pred.val-mean(pred.val)
         true.diff <- dataset[which(dataset[,2]=="val"),ncol(dataset)]- mean(dataset[which(dataset[,2]=="val"),ncol(dataset)])
         Rsquare.val=(sum(pred.diff*true.diff)/sqrt(sum(pred.diff^2)*sum(true.diff^2)))^2
         Rsquare.val
    plot(pred.val,dataset[which(dataset[,2]=="val"),ncol(dataset)],ylab=ENDPOINT)
    abline(a=0,b=1,lty=2)
    title(main=paste("validation set\nR^2 =",round(Rsquare.val,digits=2),"; RMS =",round(RMS.val,digits=2)))
    stats.Rsquare$Rsquare.val = Rsquare.val
  }

  stats.file <- paste(FILENAME,"_runs_",ii,"_model_1_Rsquare_stats.txt",sep="")
  write.table(stats.Rsquare,file=stats.file,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

  }
} # END of if (length(grep("regression",svm.type))==0)
#
graphics.off()
} ## end of function!!!!!!
#############################################################################################################
#############################################################################################################


my.density <- function(x, bw = 'nrd0', at) {
    x<-na.omit(x)

    #####
    #Borrowed from density.default for compatibility

    if (is.character(bw)) {
        if (length(x) < 2)
            stop("need at least 2 points to select a bandwidth automatically")
        bw <- switch(tolower(bw), nrd0 = bw.nrd0(x), nrd = bw.nrd(x),
            ucv = bw.ucv(x), bcv = bw.bcv(x), sj = , `sj-ste` = bw.SJ(x,
                method = "ste"), `sj-dpi` = bw.SJ(x, method = "dpi"),
            stop("unknown bandwidth rule"))
    }
    ######
    at <- matrix(at, ncol=1)
    y <- apply(at, 1, FUN=function(a, x, bw) sum(dnorm(a, x, bw)/length(x)), x=x, bw=bw )
    return(list(x=at, y=y))

}


### Function for density plots and enrichment in validation set
dens.prob <- function(TT,ENDP,PROB,my.grid=seq(0,1,0.005),THR=-1)
{
 print("---------------------------------------")
 print("- VALIDATION SET --- density analysis -")
 print("---------------------------------------")

 my.lwd=1.5

DS <- data.frame(TT,ENDP,PROB)
# my.grid=seq(0,1,0.005)
# THR=-1
# DS <- data.frame(ALL[,2],ALL[,3],ALL[,5])
##


 par(mfrow=c(3,1),oma=c(6,5,4,2.5),mar=c(0.1,0.1,0.1,0.1),las=1)
 abline.at  <- THR
 abline.col <- rep("grey",NROW(abline.at))

 val.0 <- which(DS[,1]=="val"&DS[,2]==0)
 val.1 <- which(DS[,1]=="val"&DS[,2]==1)

# adjust grid!
 print(paste("Probability values reach from",min(DS[c(val.0,val.1),3]),"to",max(DS[c(val.0,val.1),3])))
 span <- max(DS[c(val.0,val.1),3]) - min(DS[c(val.0,val.1),3])
 if (span < 0.5) {
                  if (DEBUG) {print("adjusting grid to probability values!")}
                  grid.start <- min(DS[c(val.0,val.1),3]) - span*0.1
                  if (grid.start < 0) {grid.start = 0
                               } else { grid.start = grid.start - (grid.start %% 0.005)}
                  if (DEBUG) {print(paste("new grid starts at", grid.start))}
                  grid.end <- max(DS[c(val.0,val.1),3]) + span*0.1
                  if (grid.end > 1) {grid.end = 1
                             } else {grid.end = grid.end - (grid.end %% 0.005) + 0.005}
                  if (DEBUG) {print(paste("new grid ends at", grid.end))}
                  my.grid <- seq(grid.start,grid.end,(grid.end-grid.start)/200)
                  if (DEBUG) {print(paste("new grid step size ",(grid.end-grid.start)/200))}
                  }



 N.1 <- NROW(val.1)
 N.0 <- NROW(val.0)
 R.1   <- N.1/(N.1+N.0)
 R.0   <- 1 - R.1
 # relative density
 bwd=0.025
 if (!is.na(bwd)) {
                    D.0 <- my.density(DS[val.0 ,3],at=my.grid,bw=bwd)
                    D.1 <- my.density(DS[val.1 ,3],at=my.grid,bw=bwd)
           } else {
                    D.0 <- my.density(DS[val.0 ,3],at=my.grid)
                    D.1 <- my.density(DS[val.1 ,3],at=my.grid)
                  }
 max.D.0 <- D.0$x[which(D.0$y==max(D.0$y))]
 max.D.1 <- D.1$x[which(D.1$y==max(D.1$y))]
 print(paste("maximum density of negatives at: ",max.D.0))
 print(paste("maximum density of positives at: ",max.D.1))

 D.tot <- D.1$y + D.0$y           # total density:  D.tot = D.1 + D.0
 ## cut low end of total density
 D.tot[which(D.tot<0.001)] <- 10^10
 T.0 <- D.0$y/D.tot
 T.1 <- D.1$y/D.tot

 ## cut low end of densities
 D.0$y[which(D.0$y<0.001)] <- 0
 D.1$y[which(D.1$y<0.001)] <- 0
 ## normalized enrichment
 E.0 <- (D.0$y/D.1$y)/(N.0/N.1)
 E.1 <- (D.1$y/D.0$y)/(N.1/N.0)

 ### the data.frame EE.0 contains NROW(E.0)-1 rows
 ### in row 1 we find E.0[1] in column 1 and E.0[2] in column 2
 ### in row 2 we find E.0[2] in column 1 and E.0[3] in column 2
 ### ...
 ### this can be used to find the points where the plot is crossing
 ### the 2 (3) lines, e.i. value were there are double the amount
 ### of positives than there would be in a random distribution (enrichment=2)
 EE.0 <- data.frame(E.0[1:(NROW(E.0)-1)],E.0[2:NROW(E.0)])
 ### same for E.1
 EE.1 <- data.frame(E.1[1:(NROW(E.1)-1)],E.1[2:NROW(E.1)])

##  find the middle
 j=1
 EE.0.up <- which(EE.0[,1]<=j&EE.0[,2]>=j)
 EE.0.down <- which(EE.0[,1]>=j&EE.0[,2]<=j)
 middle <- which(D.0$y==max(D.0$y[c(EE.0.up,EE.0.down)]))

 abline.at <- D.0$x[middle]
 abline.col <- "cyan"

 ## E crossing j
 j=2
 EE.0.up   <- which(EE.0[,1]<=j&EE.0[,2]>=j)
 EE.0.down <- which(EE.0[,1]>=j&EE.0[,2]<=j)
 # remove crossings if density D.1 is very small at that point!
 # very smal is here < 0.1! this might have to be adjusted
 useful.dens=0.5
 if (any(D.0$y[EE.0.up]<useful.dens)) {EE.0.up <- EE.0.up[which(D.0$y[EE.0.up]>useful.dens)]}
 if (any(D.0$y[EE.0.down]<useful.dens)) {EE.0.down <- EE.0.down[which(D.0$y[EE.0.down]>useful.dens)]}
 abline.at <- c(abline.at,D.0$x[c(EE.0.up,EE.0.down)])
 abline.col = c(abline.col,rep("grey",NROW(c(EE.0.up,EE.0.down))))
 ##
 EE.1.up <- which(EE.1[,1]<=j&EE.1[,2]>=j)
 EE.1.down <- which(EE.1[,1]>=j&EE.1[,2]<=j)
 # remove crossings if density D.1 is very small at that point!
 # very smal is here < 0.1! this might have to be adjusted
 if (any(D.1$y[EE.1.up]<useful.dens)) {EE.1.up <- EE.1.up[which(D.1$y[EE.1.up]>useful.dens)]}
 if (any(D.1$y[EE.1.down]<useful.dens)) {EE.1.down <- EE.1.down[which(D.1$y[EE.1.down]>useful.dens)]}
 abline.at <- c(abline.at,D.0$x[c(EE.1.up,EE.1.down)])
 abline.col = c(abline.col,rep("orange",NROW(c(EE.1.up,EE.1.down))))

 ##
 print("checking enrichment E.0(x) = (D.0(x)/D.1(x)) / (N.0(x)/N.1(x))")
 print("checking enrichment E.1(x) = (D.1(x)/D.0(x)) / (N.1(x)/N.0(x))")
 print("with D.1 the probability density of positive compounds and N.1 the total number of positive compounds")
 print("with D.0 the probability density of negative compounds and N.0 the total number of negative compounds")
 print("  RESULTS:")
 print(paste("central point with no enrichment and maximum density: E.0(c) = E.1(c) = 1 at c =",D.0$x[middle]))

 if (length(c(EE.1.up,EE.1.down))>0) {
      if (any(c(EE.1.up,EE.1.down) -middle < 0)) {
                                                  print("WARNING: for positive compounds there is an area with enrichment > 2 left of the central point")
                                                 }
## calculate the point where E.1 goes up over 2 on the right side of middle and closest to middle
      if (NROW(EE.1.up)> 1) {
                             EE.1.x2 <- EE.1.up[which(EE.1.up > middle)[1]]
                     } else {if (NROW(EE.1.up)==1) {EE.1.x2 <- EE.1.up
                                                   } else {EE.1.x2 <- middle}}
                                     } else {EE.1.x2 <- middle}
# EE.1.up - middle
 if (length(c(EE.0.up,EE.0.down))>0) {
     if (any(middle - c(EE.0.up,EE.0.down) < 0)) {
                                                  print("WARNING: for negative compounds there is an area with enrichment > 2 right of the central point")
                                                 }

     ## calculate the point where E.0 goes below 2 on the left side of middle and closest to middle
      if (NROW(EE.0.down)> 1) {
                               EE.0.x2 <- rev(EE.0.down[which(EE.0.down < middle)])[1]
                       } else {if (NROW(EE.0.down)==1) {EE.0.x2 <- EE.0.down
                                                       } else {EE.0.x2 <- middle}}
                              } else {EE.0.x2 <- middle}


 if (length(middle)==0)
 {
 print("WARNING: no point with E.0 = E.1 = 1 found!")
 abline.at=-1
 abline.col="black"
 }

 if (DEBUG) {print(paste("middle    = ",D.0$x[middle]))
              print(paste("EE.0.up   = ",paste(D.0$x[EE.0.up],collapse=", ")))
              print(paste("EE.0.down = ",paste(D.0$x[EE.0.down],collapse=", ")))
              print(paste("EE.0.x2   = ",D.0$x[EE.0.x2]))
              print(paste("EE.1.up   = ",paste(D.0$x[EE.1.up],collapse=", ")))
              print(paste("EE.1.down = ",paste(D.0$x[EE.1.down],collapse=", ")))
              print(paste("EE.1.x2   = ",D.0$x[EE.1.x2]))
             }


 plot(D.0,col=1,type="l",ylim=c(0,max(c(D.1$y,D.0$y))),xlab="",ylab="",xaxt="n",lwd=my.lwd)
 if (length(EE.0.x2)>0&length(EE.1.x2)>0&all(!is.na(EE.0.x2:EE.1.x2)))
       {ploygon.y <- apply(data.frame(D.0$y[EE.0.x2:EE.1.x2],D.1$y[EE.0.x2:EE.1.x2]),1, function(x) max(x))
        polygon(c(D.0$x[EE.0.x2],D.0$x[EE.0.x2:EE.1.x2],D.0$x[EE.1.x2]),c(0,ploygon.y,0),col="grey")
       }
 points(D.1,col=2,type="l",lwd=my.lwd)
 points(D.0,col=1,type="l",lwd=my.lwd)
 points(DS[val.0,3],rep(0,NROW(val.0)),pch="|")
 points(DS[val.1,3],rep(0,NROW(val.1)),pch="|",col="red")

 leg.txt <- c("0","1","class size","  0 / 1",paste(N.0,"/",N.1),paste(round(R.0*100,digits=0),"% / ",round(R.1*100,digits=0),"%",sep=""))
 if (R.1 < R.0) {legend("topright",legend=leg.txt,col=c(1,2,0,0,0,0),lty=1)
         } else {legend("topleft",legend=leg.txt,col=c(1,2,0,0,0,0),lty=1)}
 for ( j in 1:NROW(abline.at)) {abline(v=abline.at[j],col=abline.col[j])}


 plot(D.0$x,T.0,type="l",ylim=c(0,1),xlab="",ylab="",xaxt="n",lwd=my.lwd)
 if (length(EE.0.x2)>0&length(EE.1.x2)>0)
          {ploygon.y <- apply(data.frame(T.0[EE.0.x2:EE.1.x2],T.1[EE.0.x2:EE.1.x2]),1, function(x) max(x))
           polygon(c(D.0$x[EE.0.x2],D.0$x[EE.0.x2:EE.1.x2],D.0$x[EE.1.x2]),c(0,ploygon.y,0),col="grey")
          }
 points(D.0$x,T.1,type="l",col="red",lwd=my.lwd)
 points(D.0$x,T.0,type="l",col="black",lwd=my.lwd)
 points(DS[val.0,3],rep(0,NROW(val.0)),pch="|")
 points(DS[val.1,3],rep(0,NROW(val.1)),pch="|",col="red")
 abline(h=0.5,col="grey",lty=3)
 abline(h=0.66,col="violet",lty=3)
 abline(h=0.8,col="purple",lty=5)
 for ( j in 1:NROW(abline.at)) {abline(v=abline.at[j],col=abline.col[j])}

 ylim.max <- 8
 plot(D.0$x,E.0,type="l",ylim=c(0,ylim.max),xlab="",ylab="",lwd=my.lwd)
 if (length(EE.0.x2)>0&length(EE.1.x2)>0)
           {ploygon.y <- apply(data.frame(E.0[EE.0.x2:EE.1.x2],E.1[EE.0.x2:EE.1.x2]),1, function(x) max(x))
            polygon(c(D.0$x[EE.0.x2],D.0$x[EE.0.x2:EE.1.x2],D.0$x[EE.1.x2]),c(0,ploygon.y,0),col="grey")
           }
 points(D.0$x,E.1,type="l",col="red",lwd=my.lwd)
 points(D.0$x,E.0,type="l",col="black",lwd=my.lwd)
 points(DS[val.0,3],rep(0,NROW(val.0)),pch="|")
 points(DS[val.1,3],rep(0,NROW(val.1)),pch="|",col="red")
 abline(h=1,col="grey",lty=3)
 abline(h=2,col="lightgreen",lty=3)
 abline(h=3,col="green",lty=5)
 for ( j in 1:NROW(abline.at)) {abline(v=abline.at[j],col=abline.col[j])}

 mtext("Prediction Probability",1,3, outer=TRUE)
 mtext("density               ",2,3,adj=1,outer=TRUE, las=0)
 mtext("relative density",2,3, outer=TRUE,las=0)
 mtext("   norm. relative enrichment",2,3,adj=0, outer=TRUE,las=0)
 mtext("VALIDATION SET",3,2, outer=TRUE,las=0,cex=1.2)

 print("")
 print("")
}



######################################################################################################################
######################################################################################################################
model.stats <- function(Model,Dataset) {
# --- calculate a stats data.frame
PRED.train <- predict(Model,Dataset[which(Dataset[,2]=="train"),c(-1,-2,-ncol(Dataset))])
Ttrain <- table(Dataset[which(Dataset[,2]=="train"),ncol(Dataset)],PRED.train)
ACC.train <- (Ttrain[1,1]+Ttrain[2,2])/sum(Ttrain)
SP.train <- Ttrain[1,1]/(Ttrain[1,1]+Ttrain[1,2])
SE.train <- Ttrain[2,2]/(Ttrain[2,1]+Ttrain[2,2])
PPV.train <- Ttrain[2,2]/(Ttrain[1,2]+Ttrain[2,2])
NPV.train <- Ttrain[1,1]/(Ttrain[1,1]+Ttrain[2,1])
PRED.test <- predict(Model,Dataset[which(Dataset[,2]=="test"),c(-1,-2,-ncol(Dataset))])
Ttest <- table(Dataset[which(Dataset[,2]=="test"),ncol(Dataset)],PRED.test)
ACC.test <- (Ttest[1,1]+Ttest[2,2])/sum(Ttest)
SP.test <- Ttest[1,1]/(Ttest[1,1]+Ttest[1,2])
SE.test <- Ttest[2,2]/(Ttest[2,1]+Ttest[2,2])
PPV.test <- Ttest[2,2]/(Ttest[1,2]+Ttest[2,2])
NPV.test <- Ttest[1,1]/(Ttest[1,1]+Ttest[2,1])
ACC.train
ACC.test
if (DoVal == 1) {
PRED.val <- predict(Model,Dataset[which(Dataset[,2]=="val"),c(-1,-2,-ncol(Dataset))])
Tval <- table(Dataset[which(Dataset[,2]=="val"),ncol(Dataset)],PRED.val)
ACC.val <- (Tval[1,1]+Tval[2,2])/sum(Tval)
SP.val <- Tval[1,1]/(Tval[1,1]+Tval[1,2])
SE.val <- Tval[2,2]/(Tval[2,1]+Tval[2,2])
PPV.val <- Tval[2,2]/(Tval[1,2]+Tval[2,2])
NPV.val <- Tval[1,1]/(Tval[1,1]+Tval[2,1])
ACC.val
}
if (DoVal == 0) {
Stats.Out <- data.frame("ENDP"=ENDPOINT,"TT.split"=SET.column.header,"desc"=paste(DESCRIPTOR.postfixes,collapse=" "),
                        "status"=NA,"popSize"=pS,"iter"=iter,"runs"=l,
                        "testweight"=testweight,"mutationChance"=MUTATION,
                        "N.desc"=ncol(Dataset)-3,
                        "ACC.train"=round(ACC.train,digits=4),"ACC.test"=round(ACC.test,digits=4),
                        "SE.train"=round(SE.train,digits=4),"SP.train"=round(SP.train,digits=4),
                        "PPV.train"=round(PPV.train,digits=4),"NPV.train"=round(NPV.train,digits=4),
                        "SE.test"=round(SE.test,digits=4),"SP.test"=round(SP.test,digits=4),
                        "PPV.test"=round(PPV.test,digits=4),"NPV.test"=round(NPV.test,digits=4),
                        "Descriptors"=paste(colnames(Dataset)[-1*c(1,2,ncol(Dataset))],collapse=",")
                        )
} else {
Stats.Out <- data.frame("ENDP"=ENDPOINT,"TT.split"=SET.column.header,"desc"=paste(DESCRIPTOR.postfixes,collapse=" "),
                        "status"=NA,"popSize"=pS,"iter"=iter,"runs"=l,
                        "testweight"=testweight,"mutationChance"=MUTATION,
                        "N.desc"=ncol(Dataset)-3,
                        "ACC.train"=round(ACC.train,digits=4),"ACC.test"=round(ACC.test,digits=4),
                        "SE.train"=round(SE.train,digits=4),"SP.train"=round(SP.train,digits=4),
                        "PPV.train"=round(PPV.train,digits=4),"NPV.train"=round(NPV.train,digits=4),
                        "SE.test"=round(SE.test,digits=4),"SP.test"=round(SP.test,digits=4),
                        "PPV.test"=round(PPV.test,digits=4),"NPV.test"=round(NPV.test,digits=4),
                        "SE.val"=round(SE.val,digits=4),"SP.val"=round(SP.val,digits=4),
                        "PPV.val"=round(PPV.val,digits=4),"NPV.val"=round(NPV.val,digits=4),
                        "Descriptors"=paste(colnames(Dataset)[-1*c(1,2,ncol(Dataset))],collapse=",")
                        )
}
return(Stats.Out)
}


## cross validation
do.SVM.crossval <- function(comps,          # row numbers of compounds to use
                            fold=10,        # fold of cross validation
                            selected.descs  # column numbers to take as descriptors in model building
                            )
{

# small example how to combine prediction objects to a single prediciton object!
# this is just a test
if (FALSE)
{
library(e1071)
A <- data.frame("A"=c(rep(1,15),rep(0,15))+runif(30),"B"=runif(30),"C"=runif(30),"D"=c(rep(1,15),rep(0,15)))
A[,"D"] <- as.factor(A[,"D"])
A
B  <- data.frame("A"=c(1,0)+runif(2),"B"=runif(2),"C"=runif(2),"D"=c(0,1))
B
C  <- data.frame("A"=c(0,1)+runif(2),"B"=runif(2),"C"=runif(2),"D"=c(0,1))
C
model <- svm(y=A[,4],x=A[,-4], type="nu-classification",kernel="radial", probability=T,nu=0.1,gamma=0.01)
#
pred.A <- predict(model,A[,-4],probability=T)
pred.A
#
pred.B <- predict(model,B[,-4],probability=T)
class.B <- as.numeric(as.character(pred.B))
prob.B <- attributes(pred.B)$probabilities
names(attributes(pred.B))
pred.B
#
pred.C <- predict(model,C[,-4],probability=T)
class.C <- as.numeric(as.character(pred.C))
prob.C <- attributes(pred.C)$probabilities
pred.C
#
pred.D <- factor(c(class.B,class.C),level=c("0","1"))
attributes(pred.D)$names <- as.character(1:NROW(pred.D))
attributes(pred.D)$probabilities <- rbind(prob.B,prob.C)
colnames(attributes(pred.D)$probabilities)
pred.D
#
pred.D[1:2] <- pred.C
pred.D
attributes(pred.D)$probabilities[1:2,1] <- attributes(pred.C)$probabilities[,1]
attributes(pred.D)$probabilities[1:2,2] <- attributes(pred.C)$probabilities[,2]
pred.D
}
#### --- start of the main cross validation script

 #cat("doing ",fold,"-fold cross validation with ",NROW(comps)," compounds: ", round(NROW(comps)/fold,digits=1)," compounds per fold!\n",sep="")

 # build a random return prediction if error occurs
 fake.return <- factor(c(1,rep(0,NROW(comps)-1)),levels=c("0","1"))
 attributes(fake.return)$names <- as.character(1:NROW(comps))
 X <- runif(NROW(comps))
 attributes(fake.return)$probabilities <- data.frame("A"=X,"B"=1-X)
 colnames(attributes(fake.return)$probabilities) <- c("1","0")

 # comps.CV is a scrambled order of comps
 # if comps contains 100 rownumbers comps.CV contains the numbers 1 to 100 in scrambled order
 comps.CV <- order(runif(NROW(comps)))


 # class.pred.CV will contain the classes predicted by the svm model
 class.pred.CV <- rep(2,NROW(comps))
 # prob.pred.CV will contain the probabilities predicted by the svm model
 prob.pred.CV  <- data.frame("X"=rep(0,NROW(comps)),"Y"=rep(0,NROW(comps)))
# test.comps <- data.frame("comps"=comps,"match.exclude"=rep(F,NROW(comps)),"match.predict"=rep(F,NROW(comps)))


fold.size <- as.integer(NROW(comps)/fold)
if (LOOCV) {fold.size=1}
for (i in 1:fold)
{
  if (i < fold) {#cat("comps ",(i-1)*fold.size +1," to ",i*fold.size," -- ",i*fold.size - ((i-1)*fold.size),"\n")
                 comps.CV.exclude <- (((i-1)*fold.size +1):(i*fold.size))
         } else {#cat("comps ",(i-1)*fold.size +1," to ",NROW(comps)," -- ",NROW(comps) - ((i-1)*fold.size),"\n")
                 comps.CV.exclude <-(((i-1)*fold.size +1):(NROW(comps)))
                 }
  # compounds used for this training run
  CV.train.comps <- comps[comps.CV[-comps.CV.exclude]]
  
  if (!is.na(NU)) {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[CV.train.comps,c(selected.descs,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,gamma=GAMMA,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon), error=function(x) x="ERROR")
                              } else {model <- tryCatch(svm(as.formula(form),data=dataset[CV.train.comps,c(selected.descs,ncol(dataset))], type=svm.type,kernel=kernel.type,nu=NU,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),             error=function(x) x="ERROR")}
           } else {if (!is.na(GAMMA)){model <- tryCatch(svm(as.formula(form),data=dataset[CV.train.comps,c(selected.descs,ncol(dataset))], type=svm.type,kernel=kernel.type,gamma=GAMMA,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),       error=function(x) x="ERROR")
                              } else {model <- tryCatch(svm(as.formula(form),data=dataset[CV.train.comps,c(selected.descs,ncol(dataset))], type=svm.type,kernel=kernel.type,probability=T,class.weights=class.weight,tolerance=tolerance,epsilon=epsilon),                   error=function(x) x="ERROR")}}
 # check if model was built
 if (is.character(model)) {cat("ERROR in cross validation prediction during evaluation of an individual with following descriptors:\n")
                           print(colnames(dataset)[selected.descs]); return(fake.return)}
 # predict excluded compounds
 pred.CV.tmp <- tryCatch(predict(model,dataset[comps[comps.CV[comps.CV.exclude]],selected.descs],probability=T), error=function(x) x="ERROR")
         if (is.character(pred.CV.tmp)) {cat("ERROR in cross validation prediction during evaluation of an individual with following descriptors:\n")
                                         cat(colnames(dataset)[selected.descs]); return(fake.return)}

  class.pred.CV[comps.CV[comps.CV.exclude]] <- as.numeric(as.character(pred.CV.tmp))
  prob.pred.CV[comps.CV[comps.CV.exclude],1] <- attributes(pred.CV.tmp)$probabilities[,1]
  prob.pred.CV[comps.CV[comps.CV.exclude],2] <- attributes(pred.CV.tmp)$probabilities[,2]
}

# construct the output to be an object similar to a SVM prediciton object
pred.CV <- factor(class.pred.CV,levels=c("0","1"))
attributes(pred.CV)$names <- as.character(1:NROW(comps))
attributes(pred.CV)$probabilities <- prob.pred.CV
colnames(attributes(pred.CV)$probabilities) <- c("1","0")
return(pred.CV)
}

