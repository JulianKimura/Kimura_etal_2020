#import table with mean tpm values of each stage
mean_TPM<-read.table("/Users/JulianKimura/Documents/Lab/Broseq/Broseq_2019_analysis/Relative_Expression_TPM/mean_TPM_stage.txt",header = TRUE, stringsAsFactors=FALSE, row.names = 1)

df_mean_TPM<- as.data.frame(mean_TPM)

opar=par(ps=18)
label = c('EP','M', 'D', '56-75', '76-94', 'Pill', 'PH', 'PPH')
a <- c("44.69349736",	"48.60122075",	"28.04272022",	"31.28284198",	"24.83665524",	"24.38990051",	"25.85694314",	"37.53332395")
plot(a,axes=F,xlab="Stages",ylab="MeanTPM",type="b",col="violet")
axis(2)
axis(1,at=1:length(label),labels=label)
title(main = "", xlab="Stages", ylab = "MeanTPM")
legend(3,45,c("gata6-3"),col=c("violet"),lty=c(1,1,1,1,1,1,1),lwd=c(1,1,1,1,1,1,1))