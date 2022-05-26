# ptetPCP_Properties.R

setwd("/ptetPCP/")
svmProps = read.csv("./ptetPCP_csv/ptetPCP_Properties.csv", header= T)
svmProps$svm.pred = gsub(pattern = " ", replacement = "\n", svmProps$svm.pred)
svmProps$svm.pred = gsub(pattern = "-", replacement = "\n", svmProps$svm.pred)

# Enrichment of genomic features
ggplot(svmProps) + geom_count(aes(x=svm.pred, y=Chromosome))
sort(table(svmProps$Chromosome), decreasing = T)

ggplot(svmProps) + geom_violin(aes(x=svm.pred, y=mRNA_VEG)) +   scale_y_continuous(trans = "log10") + 
  xlab("Predicted Compartment") + ylab("Log mRNA Expression (VEG Growth)")
summary(aov(mRNA_VEG ~ svm.pred, data = svmProps))

ggplot(svmProps) + geom_boxplot(aes(x=svm.pred, y=nExons)) + ylim(c(0,10)) + 
  xlab("Predicted Compartment") + ylab("Exons/gene")

ggplot(svmProps) + geom_boxplot(aes(x=svm.pred, y=nIntrons)) + ylim(c(0,10)) + 
  xlab("Predicted Compartment") + ylab("Introns/gene")

ggplot(svmProps) + geom_count(aes(x=svm.pred, y=DE_Reciliation))
ggplot(svmProps, aes(x=svm.pred, y=DE_Reciliation)) + geom_bar(aes(y = (..count..)/sum(..count..)))
table(svmProps$DE_Reciliation)

ggplot(svmProps) + geom_count(aes(x=svm.pred, y=DE_TrichocystDischarge))
ggplot(svmProps, aes(x=svm.pred, y=DE_TrichocystDischarge)) + geom_bar(aes(y = (..count..)/sum(..count..)))
table(svmProps$DE_TrichocystDischarge)

ggplot(svmProps) + geom_count(aes(x=svm.pred, y=DE_Autogamy))
ggplot(svmProps, aes(x=svm.pred, y=DE_Autogamy)) + geom_bar(aes(y = (..count..)/sum(..count..)))
table(svmProps$DE_Autogamy)

# mRNA vs PSMs
svmProps$mRNA_VEG_log = log(svmProps$mRNA_VEG)
svmProps$nPSM_log = log(svmProps$nPSM)
ggplot(svmProps, aes(x=mRNA_VEG, y=nPSM)) + geom_point() + geom_smooth(method = "lm")
ggplot(svmProps, aes(x=mRNA_VEG_log, nPSM_log)) + geom_point() + geom_smooth(method = "lm") + 
  xlab("Log mRNA Expression (VEG Growth") + ylab("Log PSMs/protein")



# Trichocysts need more exploration
tmps = svmProps[grep(pattern = "*TMP*", x = svmProps$Synonyms),]
table(tmps$svm)
table(tmps$svm.pred)
table(tmps$DE_TrichocystDischarge)

# Word Clouds for Descriptions
#
# Word Clouds
library("tm")
library('wordcloud')
library(tidyverse)
source("./ptetPCP_scripts/makeWordCloud.R")

commonWords = c("protein", "domain", "family", "function", "containing", "profile", "factor", " protein", "protein ")
# makeWordCloud() make word cloud given dataframe and vector of common words
for(orgs in unique(svmProps$svm.pred)){
  png(file=paste("./ptetPCP_plots/", orgs, ".png", sep=""))
  makeWordCloud(svmProps[which(svmProps$svm.pred == orgs),], commonWords)
  graphics.off()
}


# Enrichment of protein features
ggplot(svmProps) + geom_boxplot(aes(x=svm.pred, y=nPSM))  + ylim(0,750) +
  xlab("Predicted Compartment")  + ylab("Protein Length (N amino acids)")

ggplot(svmProps) + geom_boxplot(aes(x=svm.pred, y=pI))  + ylim(0,12) +
  xlab("Predicted Compartment")  + ylab("Isoelectric Point")

ggplot(svmProps) + geom_violin(aes(x=svm.pred, y=nPSM)) + scale_y_continuous(trans = "log") + 
  xlab("Predicted Compartment")  + ylab("Log PSMs/protein")

ggplot(svmProps) + geom_count(aes(x=svm.pred, y=TransmembraneHelix, color = ..n..)) + scale_size(range=c(0,20)) + 
  xlab("Predicted Compartment") + ylab("Predicted Transmembrane Domain")
table(svmProps$TransmembraneHelix)[1] / (table(svmProps$TransmembraneHelix)[1] + table(svmProps$TransmembraneHelix)[2])

ggplot(svmProps) + geom_count(aes(x=svm.pred, y=SignalPeptide, color = ..n..)) + scale_size(range=c(0,10)) + 
  xlab("Predicted Compartment") + ylab("Predicted Signal Peptide")
table(svmProps$TransmembraneHelix)[1] / (table(svmProps$TransmembraneHelix)[1] + table(svmProps$TransmembraneHelix)[2])

ggplot(svmProps[grep("PTET*", svmProps$Accession),]) + geom_count(aes(x=svm.pred, y=TargetP, color = ..n..)) + 
  scale_size(range=c(0,20)) + 
  xlab("Predicted Compartment") + ylab("Signal Peptide, Mitochondrial Targetting Sequence, Other")



