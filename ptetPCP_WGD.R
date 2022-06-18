  # ptetPCP_WGD.R

setwd("/ptetPCP/")
load("./ptetPCP_rdata/ptetPCP_PSS.RData")


ptetWGD = read.table("/Paramecium_Paraorthologs/ptetraurelia_mac_51_annotation_v2.0.WGD.tree", header=T) 
ptetWGD$NB = NULL

for(q in 1:ncol(ptetWGD)){  # make dots into NA
  ptetWGD[,q] = sub("^[.]$", NA, x = ptetWGD[,q])
}
head(ptetWGD)

# source all pairwise combinations
source("./ptetPCP_scripts/WGDpairs.R")
source("./ptetPCP_scripts/WGDtestMakeBind")

wgdDF = data.frame(matrix(nrow=1, ncol=4))  # intiate final DF
wgdDF[1,] = c("GeneA", "GeneB", "Type", "PSS")

# create final DF
for(i in 1:nrow(ptetWGD)){
  ptetWGD_i = ptetWGD[i,]
  nProts = length(which(is.na(ptetWGD_i) == F))
  
  if(nProts > 1){
    wgdDF = WGDtestMakeBind(wgd1_1, ptetWGD_i, lDist_NA, wgd = "WGD1")
    wgdDF = WGDtestMakeBind(wgd1_2, ptetWGD_i, lDist_NA, wgd = "WGD1")
    wgdDF = WGDtestMakeBind(wgd1_3, ptetWGD_i, lDist_NA, wgd = "WGD1")
    wgdDF = WGDtestMakeBind(wgd1_4, ptetWGD_i, lDist_NA, wgd = "WGD1")
    
    wgdDF = WGDtestMakeBind(wgd2_1, ptetWGD_i, lDist_NA, wgd = "WGD2")
    wgdDF = WGDtestMakeBind(wgd2_2, ptetWGD_i, lDist_NA, wgd = "WGD2")
    wgdDF = WGDtestMakeBind(wgd2_3, ptetWGD_i, lDist_NA, wgd = "WGD2")
    wgdDF = WGDtestMakeBind(wgd2_4, ptetWGD_i, lDist_NA, wgd = "WGD2")
    wgdDF = WGDtestMakeBind(wgd2_5, ptetWGD_i, lDist_NA, wgd = "WGD2")
    wgdDF = WGDtestMakeBind(wgd2_6, ptetWGD_i, lDist_NA, wgd = "WGD2")
    wgdDF = WGDtestMakeBind(wgd2_7, ptetWGD_i, lDist_NA, wgd = "WGD2")
    wgdDF = WGDtestMakeBind(wgd2_8, ptetWGD_i, lDist_NA, wgd = "WGD2")
    
    wgdDF = WGDtestMakeBind(wgd3_1, ptetWGD_i, lDist_NA, wgd = "WGD3")
    wgdDF = WGDtestMakeBind(wgd3_2, ptetWGD_i, lDist_NA, wgd = "WGD3")
    wgdDF = WGDtestMakeBind(wgd3_3, ptetWGD_i, lDist_NA, wgd = "WGD3")
    wgdDF = WGDtestMakeBind(wgd3_4, ptetWGD_i, lDist_NA, wgd = "WGD3")
    wgdDF = WGDtestMakeBind(wgd3_5, ptetWGD_i, lDist_NA, wgd = "WGD3")
    wgdDF = WGDtestMakeBind(wgd3_6, ptetWGD_i, lDist_NA, wgd = "WGD3")
    wgdDF = WGDtestMakeBind(wgd3_7, ptetWGD_i, lDist_NA, wgd = "WGD3")
    wgdDF = WGDtestMakeBind(wgd3_8, ptetWGD_i, lDist_NA, wgd = "WGD3")
    wgdDF = WGDtestMakeBind(wgd3_9, ptetWGD_i, lDist_NA, wgd = "WGD3")
    wgdDF = WGDtestMakeBind(wgd3_10, ptetWGD_i, lDist_NA, wgd = "WGD3")
    wgdDF = WGDtestMakeBind(wgd3_11, ptetWGD_i, lDist_NA, wgd = "WGD3")
    wgdDF = WGDtestMakeBind(wgd3_12, ptetWGD_i, lDist_NA, wgd = "WGD3")
    wgdDF = WGDtestMakeBind(wgd3_13, ptetWGD_i, lDist_NA, wgd = "WGD3")
    wgdDF = WGDtestMakeBind(wgd3_14, ptetWGD_i, lDist_NA, wgd = "WGD3")
    wgdDF = WGDtestMakeBind(wgd3_15, ptetWGD_i, lDist_NA, wgd = "WGD3")
    wgdDF = WGDtestMakeBind(wgd3_16, ptetWGD_i, lDist_NA, wgd = "WGD3")
  }
}
# Clean it up
wgdDF = wgdDF[complete.cases(wgdDF),]
colnames(wgdDF) = wgdDF[1,]
wgdDF = wgdDF[2:nrow(wgdDF),]

#save.image("./ptetPCP_rdata/ptetPCP_WGD-ProteinSimilarity.RData")
load("./ptetPCP_rdata/ptetPCP_WGD-ProteinSimilarity.RData")

# Make a random DF
table(wgdDF$Type)
    # WGD1 WGD2 WGD3 
    # 195  187   58 

ranDF = data.frame(matrix(nrow=1, ncol=4))
colnames(ranDF) = c("GeneA", "GeneB", "Type", "PSS")
set.seed(42)
ranRows1 = round(runif(58, min = 1, length(lDist_NA)))  # smallest group is 58 (WGD3)
set.seed(666)
ranRows2 = round(runif(58, min = 1, length(lDist_NA)))

for(f in 1:58){
  tmpRanDF = data.frame(matrix(nrow=1, ncol=4))
  colnames(tmpRanDF) = c("GeneA", "GeneB", "Type", "PSS")
  tmpRanDF$GeneA = names(lDist_NA)[f]
  tmpRanDF$GeneB = names(lDist_NA[[f]])[f]
  tmpRanDF$Type = "Random"
  tmpRanDF$PSS = lDist_NA[[ranRows1[f]]][ranRows2[f]]
  ranDF = rbind(ranDF, tmpRanDF)
}
ranDF = ranDF[2:nrow(ranDF),]
pssWGDran = rbind(wgdDF, ranDF)

# Write
write.csv(pssWGDran, "./ptetPCP_csv/ptetPCP_WGD.R", quote = F, row.names = F)   
#save.image("./ptetPCP_rdata/ptetPCP_WGD-ProteinSimilarity.RData")

# Plots
library(ggplot2)
pssWGDran$PSS = as.numeric(pssWGDran$PSS)
ggplot(pssWGDran) + geom_violin(aes(x=Type, y=PSS)) + xlab("Ohnolog Type (or Random Pair)") + ylab("Protein Similarity Score")
ggplot(pssWGDran) + geom_boxplot(aes(x=Type, y=PSS))

load("./ptetPCP_rdata/ptetPCP_SVM.Rdata")
setwd("./ptetPCP_eucDist/NA-norm/WGD/")
for(q in 1:nrow(pssWGDran)){                               # all pairwise distance plots
  if(length(which(row.names(fData(preds2)) %in% as.character(wgdDF[q,1:2]))) == 2){
    png(paste(as.character(pssWGDran[q,1]), "_", as.character(pssWGDran[q,2]), "_hist.png", sep=""))
    plotDist(preds2[as.character(pssWGDran[q,1:2])], pcol = c("black", "red"), ylim = c(0,1))
    dev.off()  # shut off plot saving
  }
}

# Stats
wilcox.test(x=pssWGDran[which(pssWGDran$Type == "WGD1"),"PSS"],
            y=pssWGDran[which(pssWGDran$Type == "Random"),"PSS"])  # p-value = 2.2e-16

wilcox.test(x=pssWGDran[which(pssWGDran$Type == "WGD1"),"PSS"],
            y=pssWGDran[which(pssWGDran$Type == "WGD2"),"PSS"])  # p-value = .6
wilcox.test(x=pssWGDran[which(pssWGDran$Type == "WGD1"),"PSS"],
            y=pssWGDran[which(pssWGDran$Type == "WGD3"),"PSS"])  # p-value = .94
wilcox.test(x=pssWGDran[which(pssWGDran$Type == "WGD2"),"PSS"],
            y=pssWGDran[which(pssWGDran$Type == "WGD3"),"PSS"])  # p-value = .75

# MOST DIVERGENT PAIRS
#WGD1
range(round(pssWGDran[which(pssWGDran$Type == "WGD1"),'PSS'], 5))
  # 0.4613 1.0000
pssWGDran[which(pssWGDran$Type == "WGD1"),][which(round(pssWGDran[which(pssWGDran$Type == "WGD1"),'PSS'], 5) == round(0.4613, 5)),]
pssWGDran[which(pssWGDran$Type == "WGD1"),][which(round(pssWGDran[which(pssWGDran$Type == "WGD1"),'PSS'], 5) == round(1.000, 5)),]
#WGD2
range(round(pssWGDran[which(pssWGDran$Type == "WGD2"),'PSS'], 5))
# 0.39169 0.99188
pssWGDran[which(pssWGDran$Type == "WGD2"),][which(round(pssWGDran[which(pssWGDran$Type == "WGD2"),'PSS'], 5) == round(0.39169, 5)),]
pssWGDran[which(pssWGDran$Type == "WGD2"),][which(round(pssWGDran[which(pssWGDran$Type == "WGD2"),'PSS'], 5) == round(0.99188, 5)),]
#WGD3
range(round(pssWGDran[which(pssWGDran$Type == "WGD3"),'PSS'], 5))
# 0.44871 0.98905
pssWGDran[which(pssWGDran$Type == "WGD3"),][which(round(pssWGDran[which(pssWGDran$Type == "WGD3"),'PSS'], 5) == round(0.44871, 5)),]
pssWGDran[which(pssWGDran$Type == "WGD3"),][which(round(pssWGDran[which(pssWGDran$Type == "WGD3"),'PSS'], 5) == round(0.98905, 5)),]

# Write CSV: Top 5 Most Divergent of Each Type with SVM Preds
svmProps = read.csv("/ptetPCP/ptetPCP_csv/ptetPCP_Properties.csv", header= T)
svmProps$svm.pred = gsub(pattern = " ", replacement = "\n", svmProps$svm.pred)
svmProps$svm.pred = gsub(pattern = "-", replacement = "\n", svmProps$svm.pred)

mostDivgernt = merge(rbind(
  pssWGDran[which(pssWGDran$PSS %in% sort(pssWGDran[which(pssWGDran$Type == "WGD1"),'PSS'])[1:5]),],
  pssWGDran[which(pssWGDran$PSS %in% sort(pssWGDran[which(pssWGDran$Type == "WGD2"),'PSS'])[1:5]),],
  pssWGDran[which(pssWGDran$PSS %in% sort(pssWGDran[which(pssWGDran$Type == "WGD3"),'PSS'])[1:5]),]), svmProps, by.x="GeneA", by.y = "Accession")
mostDivgernt = merge(mostDivgernt, svmProps, by.x="GeneB", by.y="Accession")

write.csv(mostDivgernt[,c(1:4, 14:16, 52:54)], 
  file = "/ptetPCP/ptetPCP_csv/ptetPCP_WGD_mostDivergent.csv", quote = F, row.names = F)

# Get sequences
library(seqinr)
load("./ptetPCP_rdata/ptetPCP_SVM.Rdata")

ptetProts = read.fasta("/Paramecium_Proteomes-Predicted/ptetraurelia_mac_51_annotation_v2.0.protein.fa", seqtype = "AA", as.string = T)
setwd("/ptetPCP/ptetPCP_fasta/WGD/")
for(q in 1:nrow(pssWGDran)){
  seqinr::write.fasta(sequences = getSequence(ptetProts[c(pssWGDran$GeneA[q], pssWGDran$GeneB[q])], as.string = T), 
                      names = getName(ptetProts[c(pssWGDran$GeneA[q], pssWGDran$GeneB[q])]), 
                      file.out = paste(pssWGDran$GeneA[q], "_", pssWGDran$GeneB[q], ".fasta", sep="")
                      )
}

# Make alignments + measure evolutionary distance
  # align in bash with muscle
    # muscle -in $fasta -out $outfile
  # distance in bash with megacc
    # megacc -a ../distance_estimation_pairwise_amino_acid.mao -d $fasta -o $outfile

megFiles = list.files("/ptetPCP/ptetPCP_mega/ptetPCP_WGD-dist/", full.names = T)
megFiles = megFiles[grep(pattern = "*[0-9].meg", x = megFiles)]

FmegTab = data.frame(matrix(ncol = 3, nrow = 1))
colnames(FmegTab) = c("GeneA", "GeneB", "Sequence_Distance")

for(meg in megFiles){
  megTab = read.table(meg, skip = 37)[2]  # .meg files with 2 seqs have 37 lines of unimportance
  colnames(megTab) = "Sequence_Distance"
  megTab$GeneA = unlist(strsplit(x = substr(x = basename(meg), start = 1, 37), "_"))[1]
  megTab$GeneB = unlist(strsplit(x = substr(x = basename(meg), start = 1, 37), "_"))[2]
  FmegTab = rbind(FmegTab, megTab)
}
FmegTab = FmegTab[2:nrow(FmegTab),]
tmp = merge(FmegTab, wgdDF, by="GeneA")  # merge mega distances and PSS
megaEucDF = tmp[-which(tmp$GeneB.x != tmp$GeneB.y),]  # artefact of merging is removed and cp'ed into new var
megaEucDF = as.data.frame(cbind(megaEucDF$Type, megaEucDF$GeneA, megaEucDF$GeneB.x, megaEucDF$Sequence_Distance, megaEucDF$PSS))  # rearrange
colnames(megaEucDF) = c("Type", "GeneA", "GeneB", "Sequence_Distance", "PSS")
megaEucDF$Sequence_Distance = as.numeric(megaEucDF$Sequence_Distance)
megaEucDF$PSS = as.numeric(megaEucDF$PSS)

write.csv(megaEucDF, "/ptetPCP/ptetPCP_csv/ptetPCP_WGD_seqPSS.csv", quote = F, row.names = F)

# Plots
library(tidyverse)
ggplot(megaEucDF, aes(x=Sequence_Distance, y=PSS)) + 
  geom_point() + geom_smooth(method = "lm", col = "black") + 
  xlim(0,1.5) + ylim(0,1) + xlab("P-Distance") + ylab("Protein Similarity Score") + facet_wrap(~ Type)

# Stats
WGD1.lm = lm(PSS ~ Sequence_Distance, data=megaEucDF[which(megaEucDF$Type == "WGD1"),]) 
summary(WGD1.lm) # p-value: 0.1001
WGD2.lm = lm(PSS ~ Sequence_Distance, data=megaEucDF[which(megaEucDF$Type == "WGD2"),]) 
summary(WGD2.lm) # p-value: p-value: 1.887e-06
WGD3.lm = lm(PSS ~ Sequence_Distance, data=megaEucDF[which(megaEucDF$Type == "WGD3"),]) 
summary(WGD3.lm) # p-value: p-value: 0.006464


plot(megaEucDF[which(megaEucDF$Type == "WGD1"),'Sequence_Distance'], resid(WGD1.lm), 
     ylab="Residuals", xlab="WGD1: P-Distance")
abline(0, 0)  
plot(megaEucDF[which(megaEucDF$Type == "WGD2"),'Sequence_Distance'], resid(WGD2.lm), 
     ylab="Residuals", xlab="WGD2: P-Distance") 
abline(0, 0) 
plot(megaEucDF[which(megaEucDF$Type == "WGD3"),'Sequence_Distance'], resid(WGD3.lm), 
     ylab="Residuals", xlab="WGD3: P-Distance") 
abline(0, 0) 

# Compartment enrichment?
setwd("/ptetPCP/")
svmProps = read.csv("./ptetPCP_csv/ptetPCP_Properties.csv", header= T)
svmProps$svm.pred = gsub(pattern = " ", replacement = "\n", svmProps$svm.pred)
svmProps$svm.pred = gsub(pattern = "-", replacement = "\n", svmProps$svm.pred)

ggplot(svmProps) + geom_count(aes(x=svm.pred, y=Has.WGD1.paralog..T.F.)) +   scale_y_continuous(trans = "log10") + 
  xlab("Predicted Compartment") + ylab("Log mRNA Expression (VEG Growth)")

write.csv(cbind(data.frame(svmProps[which(complete.cases(svmProps[,'Has.WGD1.paralog..T.F.'])),] %>% 
             group_by(svm.pred) %>% 
             summarise(WGD1 = sum(Has.WGD1.paralog..T.F.==TRUE)/(sum(Has.WGD1.paralog..T.F.==TRUE) + sum(Has.WGD1.paralog..T.F.==FALSE)))),
data.frame(svmProps[which(complete.cases(svmProps[,'Has.WGD2.paralog..T.F.'])),] %>% 
             group_by(svm.pred) %>% 
             summarise(WGD2 = sum(Has.WGD2.paralog..T.F.==TRUE)/(sum(Has.WGD2.paralog..T.F.==TRUE) + sum(Has.WGD2.paralog..T.F.==FALSE)))),
data.frame(svmProps[which(complete.cases(svmProps[,'Has.WGD3.paralog..T.F.'])),] %>% 
             group_by(svm.pred) %>% 
             summarise(WGD3 = sum(Has.WGD3.paralog..T.F.==TRUE)/(sum(Has.WGD3.paralog..T.F.==TRUE) + sum(Has.WGD3.paralog..T.F.==FALSE))))),
file = "/ptetPCP/ptetPCP_csv/ptetPCP_WGG-compartment_enrichment.csv", quote = F, row.names = F)

# Evo Translocations
mostDivgernt_all = merge(rbind(
  pssWGDran[which(pssWGDran$PSS %in% sort(pssWGDran[which(pssWGDran$Type == "WGD1"),'PSS'])),],
  pssWGDran[which(pssWGDran$PSS %in% sort(pssWGDran[which(pssWGDran$Type == "WGD2"),'PSS'])),],
  pssWGDran[which(pssWGDran$PSS %in% sort(pssWGDran[which(pssWGDran$Type == "WGD3"),'PSS'])),]), svmProps, by.x="GeneA", by.y = "Accession")
mostDivgernt_all = merge(mostDivgernt_all, svmProps, by.x="GeneB", by.y="Accession")

ggplot(mostDivgernt_all) + geom_count(aes(x=svm.pred.x, y=svm.pred.y, color = ..n..)) + scale_size(range=c(0,10)) + 
  xlab("Predicted Compartment (Gene A)") + ylab("Predicted Compartment (Gene B)")
ggplot(mostDivgernt_all) + geom_count(aes(x=svm.pred.x, y=svm.pred.y, color = ..n..)) + scale_size(range=c(0,10)) + facet_wrap(~ Type)
  # meh


