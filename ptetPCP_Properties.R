# ptetPCP_Properties.R

library(tidyverse)
setwd("/ptetPCP/")
svmProps = read.csv("./ptetPCP_csv/ptetPCP_Properties.csv", header= T)
svmProps$svm.pred = gsub(pattern = " ", replacement = "\n", svmProps$svm.pred)
svmProps$svm.pred = gsub(pattern = "-", replacement = "\n", svmProps$svm.pred)

# Quick summary of all  features
allSummaries = by(svmProps[grep(pattern = "PTET*", x = svmProps$Accession),], factor(svmProps[grep(pattern = "PTET*", x = svmProps$Accession),'svm.pred']), summary)
allSummariesDF = do.call(rbind.data.frame, allSummaries)
write.csv(allSummariesDF[grep("Median*", allSummariesDF$Freq),], 
          file = "/ptetPCP/ptetPCP_csv/ptetPCP_Properties-Summary-Median.csv", quote = F)
write.csv(merge(merge(data.frame(svmProps[grep(pattern = "PTET*", x = svmProps$Accession),] %>% group_by(svm.pred) %>% summarise(TMdomain = sum(TransmembraneHelix==T)/(sum(TransmembraneHelix==T) + sum(TransmembraneHelix==F)))),
                      data.frame(svmProps[grep(pattern = "PTET*", x = svmProps$Accession),] %>% group_by(svm.pred) %>% summarise(SignalP = sum(TargetP=="SP")/(sum(TargetP=="SP") + sum(TargetP=="OTHER")))), by="svm.pred"), 
                data.frame(svmProps[grep(pattern = "PTET*", x = svmProps$Accession),] %>% group_by(svm.pred) %>% summarise(mTP = sum(TargetP=="mTP")/(sum(TargetP=="mTP") + sum(TargetP=="OTHER")))), by="svm.pred"), 
          file = "/ptetPCP/ptetPCP_csv/ptetPCP_Properties-Summary-TrueFalse.csv", quote = F, row.names = F)


##### Enrichment of Features
#  Genomic features
ggplot(svmProps) + geom_count(aes(x=svm.pred, y=Chromosome))
sort(table(svmProps$Chromosome), decreasing = T)
  # Uninteresting

ggplot(svmProps) + geom_violin(aes(x=svm.pred, y=mRNA_VEG)) +   scale_y_continuous(trans = "log10") + 
  xlab("Predicted Compartment") + ylab("Log mRNA Expression (VEG Growth)")
summary(aov(mRNA_VEG ~ svm.pred, data = svmProps))
TukeyHSD(aov(mRNA_VEG ~ svm.pred, data = svmProps))
write.csv(data.frame(TukeyHSD(aov(mRNA_VEG ~ svm.pred, svmProps))$'svm.pred'), 
          file = "/ptetPCP/ptetPCP_csv/ptetPCP_mRNA_Tukey.csv", quote = F, row.names = T)

ggplot(svmProps) + geom_boxplot(aes(x=svm.pred, y=nExons)) + ylim(c(0,10)) + 
  xlab("Predicted Compartment") + ylab("Exons/gene")
    # meh

ggplot(svmProps) + geom_boxplot(aes(x=svm.pred, y=nIntrons)) + ylim(c(0,10)) + 
  xlab("Predicted Compartment") + ylab("Introns/gene")
    # meh

ggplot(svmProps) + geom_count(aes(x=svm.pred, y=DE_Reciliation))
write.csv(data.frame(svmProps[which(complete.cases(svmProps[,'DE_Reciliation'])),] %>% group_by(svm.pred) %>% summarise(DErecilliation = sum(DE_Reciliation==TRUE)/(sum(DE_Reciliation==TRUE) + sum(DE_Reciliation==FALSE)))),
          file="/ptetPCP/ptetPCP_csv/ptetPCP_Properties_DE-recilia.csv", quote = F, row.names = F)
  # ~5% of unknown: Basal Body 2 is ~16%

ggplot(svmProps) + geom_count(aes(x=svm.pred, y=DE_TrichocystDischarge))
write.csv(data.frame(svmProps[which(complete.cases(svmProps[,'DE_TrichocystDischarge'])),] %>% group_by(svm.pred) %>% summarise(DEtrich = sum(DE_TrichocystDischarge==TRUE)/(sum(DE_TrichocystDischarge==TRUE) + sum(DE_TrichocystDischarge==FALSE)))),
        file="/ptetPCP/ptetPCP_csv/ptetPCP_Properties_DE-trich.csv", quote = F, row.names = F)
  # ~2.4% of unknown: 6% of Trichocyst Matrix

ggplot(svmProps) + geom_count(aes(x=svm.pred, y=DE_Autogamy))
data.frame(svmProps[which(complete.cases(svmProps[,'DE_Autogamy'])),] %>% group_by(svm.pred) %>% summarise(DE_Auto = sum(DE_Autogamy==TRUE)/(sum(DE_Autogamy==TRUE) + sum(DE_Autogamy==FALSE))))
    # meh

ggplot(svmProps) + geom_count(aes(x=svm.pred, y=ExpressionGroup))
data.frame(svmProps[which(complete.cases(svmProps[,'ExpressionGroup'])),] %>% group_by(svm.pred) %>% summarise(DE_Auto = sum(DE_Autogamy==TRUE)/(sum(DE_Autogamy==TRUE) + sum(DE_Autogamy==FALSE))))
  # meh


#  Proteomic features
ggplot(svmProps) + geom_boxplot(aes(x=svm.pred, y=nAA))  + ylim(0,3500) +
  xlab("Predicted Compartment")  + ylab("Protein Length (N amino acids)")
summary(aov(nAA ~ svm.pred, data = svmProps))
TukeyHSD(aov(nAA ~ svm.pred, data = svmProps))
write.csv(data.frame(TukeyHSD(aov(nAA ~ svm.pred, data = svmProps))$'svm.pred'), 
          file = "/ptetPCP/ptetPCP_csv/ptetPCP_nAA.csv", quote = F, row.names = T)

ggplot(svmProps) + geom_boxplot(aes(x=svm.pred, y=pI))  + ylim(0,12) +
  xlab("Predicted Compartment")  + ylab("Isoelectric Point")
summary(aov(pI ~ svm.pred, svmProps))
TukeyHSD(aov(pI ~ svm.pred, svmProps))
write.csv(data.frame(TukeyHSD(aov(pI ~ svm.pred, svmProps))$'svm.pred'), 
          file = "/ptetPCP/ptetPCP_csv/ptetPCP_pI_Tukey.csv", quote = F, row.names = T)

ggplot(svmProps) + geom_violin(aes(x=svm.pred, y=nPSM)) + scale_y_continuous(trans = "log") + 
  xlab("Predicted Compartment")  + ylab("Log PSMs/protein")
summary(aov(nPSM ~ svm.pred, svmProps))
TukeyHSD(aov(nPSM ~ svm.pred, svmProps))
write.csv(data.frame(TukeyHSD(aov(nPSM ~ svm.pred, svmProps))$'svm.pred'), 
          file = "/ptetPCP/ptetPCP_csv/ptetPCP_nPSM_Tukey.csv", quote = F, row.names = T)

ggplot(svmProps) + geom_violin(aes(x=svm.pred, y=mRNA_VEG)) + scale_y_continuous(trans = "log") + 
  xlab("Predicted Compartment")  + ylab("Log mRNA Expression (VEG)")
summary(aov(mRNA_VEG ~ svm.pred, svmProps))
TukeyHSD(aov(mRNA_VEG ~ svm.pred, svmProps))
write.csv(data.frame(TukeyHSD(aov(mRNA_VEG ~ svm.pred, svmProps))$'svm.pred'), 
          file = "/ptetPCP/ptetPCP_csv/ptetPCP_mRNAveg_Tukey.csv", quote = F, row.names = T)

ggplot(svmProps) + geom_count(aes(x=svm.pred, y=TransmembraneHelix, color = ..n..)) + scale_size(range=c(0,20)) + 
  xlab("Predicted Compartment") + ylab("Predicted Transmembrane Domain")
table(svmProps$TransmembraneHelix)

ggplot(svmProps) + geom_count(aes(x=svm.pred, y=SignalPeptide, color = ..n..)) + scale_size(range=c(0,10)) + 
  xlab("Predicted Compartment") + ylab("Predicted Signal Peptide")
table(svmProps$TransmembraneHelix)[1] / (table(svmProps$TransmembraneHelix)[1] + table(svmProps$TransmembraneHelix)[2])

ggplot(svmProps[grep("PTET*", svmProps$Accession),]) + geom_count(aes(x=svm.pred, y=TargetP, color = ..n..)) + 
  scale_size(range=c(0,20)) + 
  xlab("Predicted Compartment") + ylab("Signal Peptide, Mitochondrial Targetting Sequence, Other")


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
length(which(tmps$DE_TrichocystDischarge == T)) / (length(which(tmps$DE_TrichocystDischarge == T)) + length(which(tmps$DE_TrichocystDischarge == F)))
  # 4.6%

#####
##### Word Clouds
library("tm")
library('wordcloud')
library(tidyverse)
source("./ptetPCP_scripts/makeWordCloud.R")

commonWords = c("coil", " coil", "coil ", 
                "coiled",  " coiled", "coiled ", "  coiled", "coiled  ", " coiled ", 
                "  coiled  ", "  coiled  ", "Coiled", "Coiled ", " Coiled",
                "fold", "subunit", "protein", "domain", "family", "function", "containing", "profile", "factor", 
                " protein", "protein ", "protein  ", "  protein", "Protein")
# makeWordCloud() make word cloud given dataframe and vector of common words
for(orgs in unique(svmProps$svm.pred)){
  png(file=paste("/ptetPCP/ptetPCP_plots/wordcloud/", orgs, ".png", sep=""))
  makeWordCloud(svmProps[which(svmProps$svm.pred == orgs),], commonWords)
  graphics.off()
}


