# Install the ROCR package
install.packages('ROCR')
# loading ROCR library
library(ROCR)
# load data
data(ROCR.hiv)
attach(ROCR.hiv)

# ROCR.hiv$hiv.svm contains svm classification data
pred.svm <- prediction(hiv.svm$predictions, hiv.svm$labels)

# calculating the values for ROC curve,auc,accuracy,specificity,sensitivity
pred <- prediction(target_pred, target_class)
perf.1 <- performance(pred.svm,measure="tpr",x.measure="fpr")
perf.2 <- performance(pred.svm,measure="acc",x.measure="cutoff")
perf.3 <- performance(pred.svm,measure="auc")
perf.4 <- performance(pred.svm,measure="spec")
perf.5 <- performance(pred.svm,measure="sens")

# changing params for the ROC plot - width, etc
par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=1.3,cex.lab=1.4)

# plotting the ROC curve
plot(perf.1,col="red",lty=3, lwd=3)

## Accuracy plot vs cutoff plot
# plot(perf.2,lty=3,lwd=0.5,col="lightpink")
# abline(v=0.2213)
# plot(perf.2,avg="vertical",col="lightcoral",spread.estimate="boxplot",box.lty=7, box.lwd=5, box.col="slategray4",lwd=3,add=T)

# converting S4 class to vector
perf.3  <- unlist(slot(perf.3 , "y.values"))
# adding min and max ROC AUC to the center of the plot
minauc<-min(round(auc, digits = 2))
maxauc<-max(round(auc, digits = 2))
minauct <- paste(c("min(AUC)  = "),minauc,sep="")
maxauct <- paste(c("max(AUC) = "),maxauc,sep="")
legend(0.3,0.6,c(minauct,maxauct,"\n"),border="white",cex=1.7,box.col = "white")
