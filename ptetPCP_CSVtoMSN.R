# ptetPCP_CSVMtoMSN.R

library(dplyr)
library(pRoloc)

setwd("/ptetPCP/")

ptetPCP = read.csv(file = "ptetPCP_csv/ptetPCP.csv", header=T) 
dim(ptetPCP)  # 12579 proteins X 248 columns

### Remove cRaP Proteins
ptetPCP = rbind(ptetPCP[grep(pattern = "PTET*", x = ptetPCP$Accession),],   
                ptetPCP[grep(pattern = "YP_00*", x = ptetPCP$Accession),],  
                ptetPCP[grep(pattern = "*lcl", x = ptetPCP$Accession),])   
dim(ptetPCP)  # 12545 proteins X 248 columns (34 proteins removed)
# 11564 ptet proteins, 617 Kleb proteins, 31 mito ORFs

### Change Mitochondrial ORF names
ptetPCP[grep(pattern = "*lcl", x = ptetPCP$Accession),'Accession'] = ptetPCP[grep(pattern = "*lcl", x = ptetPCP$Accession),'Description']
ptetPCP[grep(pattern = "ORF*", x = ptetPCP$Accession),'Accession'] = unlist(lapply(ptetPCP[grep(pattern = "ORF*", x = ptetPCP$Accession),'Accession'], strsplit, " "))[grep(pattern = "ORF*", x = unlist(lapply(ptetPCP[grep(pattern = "ORF*", x = ptetPCP$Accession),'Accession'], strsplit, " ")))]
ptetPCP[grep(pattern = "ORF*", x = ptetPCP$Accession),'Accession'] = unlist(lapply(strsplit(ptetPCP[grep(pattern = "ORF*", x = ptetPCP$Accession),'Accession'], ":"), "[", 1))

### Remove Low Confidence
ptetPCP = ptetPCP[-which(ptetPCP$Protein.FDR.Confidence..Combined == "Low"),]
dim(ptetPCP)  # 12212 proteins X 248 columns (333 proteins removed)

### Remove Proteins with < 3 PSMs
ptetPCP = ptetPCP[-which(ptetPCP$X..PSMs < 3),]
dim(ptetPCP)  # 11437 proteins X 248 columns (775 proteins removed)

### Remove Proteins without Unique Peptides
ptetPCP = ptetPCP[-which(ptetPCP$X..Unique.Peptides == 0),]
dim(ptetPCP)  # 11436 proteins X 248 columns (1 protein removed)
# 10887 ptet proteins, 518 Kleb proteins, 31 mito ORFs

### Raw abundance file
# Notes: 1 5K is missing (corrupt), 1 12K was rerun due to machine error (F146 --> F148)
ptetPCP_raw = ptetPCP %>% dplyr::select(colnames(ptetPCP)[c(2,1,4:9)], colnames(ptetPCP)[34:140])
dim(ptetPCP_raw)  # 11436 proteins by 115 columns
write.csv(ptetPCP_raw, "./ptetPCP_csv/ptetPCP_raw.csv", row.names = F, quote = F)

### Split Replicate Experiments
ptetPCP_raw_MAC = ptetPCP_raw[, grep(pattern = "*Mac", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_MAC = as.data.frame(rowMeans(ptetPCP_raw_MAC[,c(1:3)]))
ptetPCP_raw_rep2_MAC = as.data.frame(rowMeans(ptetPCP_raw_MAC[,c(4:6)]))
ptetPCP_raw_rep3_MAC = as.data.frame(rowMeans(ptetPCP_raw_MAC[,c(7:9)]))

ptetPCP_raw_300g = ptetPCP_raw[, grep(pattern = "*300g", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_300g = as.data.frame(rowMeans(ptetPCP_raw_300g[,c(1:3)]))
ptetPCP_raw_rep2_300g = as.data.frame(rowMeans(ptetPCP_raw_300g[,c(4:6)]))
ptetPCP_raw_rep3_300g = as.data.frame(rowMeans(ptetPCP_raw_300g[,c(7:9)]))

ptetPCP_raw_1K = ptetPCP_raw[, grep(pattern = "*1K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_1K = as.data.frame(rowMeans(ptetPCP_raw_1K[,c(1:3)]))
ptetPCP_raw_rep2_1K = as.data.frame(rowMeans(ptetPCP_raw_1K[,c(4:6)]))
ptetPCP_raw_rep3_1K = as.data.frame(rowMeans(ptetPCP_raw_1K[,c(7:9)]))

ptetPCP_raw_3K = ptetPCP_raw[, grep(pattern = "*3K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_3K = as.data.frame(rowMeans(ptetPCP_raw_3K[,c(1:3)]))
ptetPCP_raw_rep2_3K = as.data.frame(rowMeans(ptetPCP_raw_3K[,c(4:6)]))
ptetPCP_raw_rep3_3K = as.data.frame(rowMeans(ptetPCP_raw_3K[,c(7:9)]))

ptetPCP_raw_5K = ptetPCP_raw[, grep(pattern = "*..Sample..5K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_5K = as.data.frame(rowMeans(ptetPCP_raw_5K[,c(1:3)]))
ptetPCP_raw_rep2_5K = as.data.frame(rowMeans(ptetPCP_raw_5K[,c(4:6)]))
ptetPCP_raw_rep3_5K = as.data.frame(rowMeans(ptetPCP_raw_5K[,c(7:8)]))  # Missing one file

ptetPCP_raw_9K = ptetPCP_raw[, grep(pattern = "*..Sample..9K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_9K = as.data.frame(rowMeans(ptetPCP_raw_9K[,c(1:3)]))
ptetPCP_raw_rep2_9K = as.data.frame(rowMeans(ptetPCP_raw_9K[,c(4:6)]))
ptetPCP_raw_rep3_9K = as.data.frame(rowMeans(ptetPCP_raw_9K[,c(7:9)]))

ptetPCP_raw_12K = ptetPCP_raw[, grep(pattern = "*12K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_12K = as.data.frame(rowMeans(ptetPCP_raw_12K[,c(1:3)]))
ptetPCP_raw_rep2_12K = as.data.frame(rowMeans(ptetPCP_raw_12K[,c(4:6)]))
ptetPCP_raw_rep3_12K = as.data.frame(rowMeans(ptetPCP_raw_12K[,c(7:9)]))

ptetPCP_raw_15K = ptetPCP_raw[, grep(pattern = "*15K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_15K = as.data.frame(rowMeans(ptetPCP_raw_15K[,c(1:3)]))
ptetPCP_raw_rep2_15K = as.data.frame(rowMeans(ptetPCP_raw_15K[,c(4:6)]))
ptetPCP_raw_rep3_15K = as.data.frame(rowMeans(ptetPCP_raw_15K[,c(7:9)]))

ptetPCP_raw_30K = ptetPCP_raw[, grep(pattern = "*30K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_30K = as.data.frame(rowMeans(ptetPCP_raw_30K[,c(1:3)]))
ptetPCP_raw_rep2_30K = as.data.frame(rowMeans(ptetPCP_raw_30K[,c(4:6)]))
ptetPCP_raw_rep3_30K = as.data.frame(rowMeans(ptetPCP_raw_30K[,c(7:9)]))

ptetPCP_raw_79K = ptetPCP_raw[, grep(pattern = "*79K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_79K = as.data.frame(rowMeans(ptetPCP_raw_79K[,c(1:3)]))
ptetPCP_raw_rep2_79K = as.data.frame(rowMeans(ptetPCP_raw_79K[,c(4:6)]))
ptetPCP_raw_rep3_79K = as.data.frame(rowMeans(ptetPCP_raw_79K[,c(7:9)]))

ptetPCP_raw_120K = ptetPCP_raw[, grep(pattern = "*120K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_120K = as.data.frame(rowMeans(ptetPCP_raw_120K[,c(1:3)]))
ptetPCP_raw_rep2_120K = as.data.frame(rowMeans(ptetPCP_raw_120K[,c(4:6)]))
ptetPCP_raw_rep3_120K = as.data.frame(rowMeans(ptetPCP_raw_120K[,c(7:9)]))

ptetPCP_raw_SUP = ptetPCP_raw[, grep(pattern = "*Sup", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_SUP = as.data.frame(rowMeans(ptetPCP_raw_SUP[,c(1:3)]))
ptetPCP_raw_rep2_SUP = as.data.frame(rowMeans(ptetPCP_raw_SUP[,c(4:6)]))
ptetPCP_raw_rep3_SUP = as.data.frame(rowMeans(ptetPCP_raw_SUP[,c(7:9)]))

ptetPCP_raw_rep1 = cbind(ptetPCP_raw[,c(1:8)], ptetPCP_raw_rep1_MAC, ptetPCP_raw_rep1_300g, ptetPCP_raw_rep1_1K, ptetPCP_raw_rep1_3K, ptetPCP_raw_rep1_5K, ptetPCP_raw_rep1_9K, ptetPCP_raw_rep1_12K, ptetPCP_raw_rep1_15K, ptetPCP_raw_rep1_30K, ptetPCP_raw_rep1_79K, ptetPCP_raw_rep1_120K, ptetPCP_raw_rep1_SUP)
ptetPCP_raw_rep2 = cbind(ptetPCP_raw[,c(1:8)], ptetPCP_raw_rep2_MAC, ptetPCP_raw_rep2_300g, ptetPCP_raw_rep2_1K, ptetPCP_raw_rep2_3K, ptetPCP_raw_rep2_5K, ptetPCP_raw_rep2_9K, ptetPCP_raw_rep2_12K, ptetPCP_raw_rep2_15K, ptetPCP_raw_rep2_30K, ptetPCP_raw_rep2_79K, ptetPCP_raw_rep2_120K, ptetPCP_raw_rep2_SUP)
ptetPCP_raw_rep3 = cbind(ptetPCP_raw[,c(1:8)], ptetPCP_raw_rep3_MAC, ptetPCP_raw_rep3_300g, ptetPCP_raw_rep3_1K, ptetPCP_raw_rep3_3K, ptetPCP_raw_rep3_5K, ptetPCP_raw_rep3_9K, ptetPCP_raw_rep3_12K, ptetPCP_raw_rep3_15K, ptetPCP_raw_rep3_30K, ptetPCP_raw_rep3_79K, ptetPCP_raw_rep3_120K, ptetPCP_raw_rep3_SUP)
colnames(ptetPCP_raw_rep1) = c("Accession", "Confidence-FDR", "Coverage", "nPeptides", "nPSM", "nUniquePeptides","nAA", "nRazorPeptides", 
                               "MAC-1", "300g-1", "1K-1", "3K-1", "5K-1", "9K-1", "12K-1", "15K-1", "30K-1", "79K-1", "120K-1", "SUP-1")
colnames(ptetPCP_raw_rep2) = c("Accession", "Confidence-FDR", "Coverage", "nPeptides", "nPSM", "nUniquePeptides","nAA", "nRazorPeptides", 
                               "MAC-2", "300g-2", "1K-2", "3K-2", "5K-2", "9K-2", "12K-2", "15K-2", "30K-2", "79K-2", "120K-2", "SUP-2")
colnames(ptetPCP_raw_rep3) = c("Accession", "Confidence-FDR", "Coverage", "nPeptides", "nPSM", "nUniquePeptides","nAA", "nRazorPeptides", 
                               "MAC-3", "300g-3", "1K-3", "3K-3", "5K-3", "9K-3", "12K-3", "15K-3", "30K-3", "79K-3", "120K-3", "SUP-3")

### Remove + Output Sup only proteins
source("/ptetPCP/ptetPCP_scripts/NAprocessPCP_SUP.R")
ptetPCP_raw_rep1_SUPonly = NAprocessPCP_SUP(ptetPCP_raw_rep1)  # 333 proteins
ptetPCP_raw_rep2_SUPonly = NAprocessPCP_SUP(ptetPCP_raw_rep2)  # 235 proteins
ptetPCP_raw_rep3_SUPonly = NAprocessPCP_SUP(ptetPCP_raw_rep3)  # 277 proteins
ptetPCP_raw_SUPonly = intersect(intersect(ptetPCP_raw_rep1_SUPonly$Accession, ptetPCP_raw_rep2_SUPonly$Accession), ptetPCP_raw_rep3_SUPonly$Accession)
# 107 proteins shared
write.csv(ptetPCP_raw_SUPonly, file = "/ptetPCP/ptetPCP_csv/ptetPCP_SupernatantONLY.csv", quote = F, row.names = F)

ptetPCP_raw_rep1 = ptetPCP_raw_rep1[-which(ptetPCP_raw_rep1$Accession %in% ptetPCP_raw_rep1_SUPonly$Accession),]
ptetPCP_raw_rep2 = ptetPCP_raw_rep2[-which(ptetPCP_raw_rep2$Accession %in% ptetPCP_raw_rep2_SUPonly$Accession),]
ptetPCP_raw_rep3 = ptetPCP_raw_rep3[-which(ptetPCP_raw_rep3$Accession %in% ptetPCP_raw_rep3_SUPonly$Accession),]

### Remove MAC only proteins
source("/ptetPCP/ptetPCP_scripts/NAprocessPCP_MAC.R")
ptetPCP_raw_rep1_MAConly = NAprocessPCP_MAC(ptetPCP_raw_rep1)  # 55 proteins
ptetPCP_raw_rep2_MAConly = NAprocessPCP_MAC(ptetPCP_raw_rep2)  # 64 proteins
ptetPCP_raw_rep3_MAConly = NAprocessPCP_MAC(ptetPCP_raw_rep3)  # 69 proteins
ptetPCP_raw_MAConly = intersect(intersect(ptetPCP_raw_rep1_MAConly$Accession, ptetPCP_raw_rep2_MAConly$Accession), ptetPCP_raw_rep3_MAConly$Accession)

ptetPCP_raw_rep1 = ptetPCP_raw_rep1[-which(ptetPCP_raw_rep1$Accession %in% NAprocessPCP_MAC(ptetPCP_raw_rep1)[,1]),]
ptetPCP_raw_rep2 = ptetPCP_raw_rep2[-which(ptetPCP_raw_rep2$Accession %in% NAprocessPCP_MAC(ptetPCP_raw_rep2)[,1]),]
ptetPCP_raw_rep3 = ptetPCP_raw_rep3[-which(ptetPCP_raw_rep3$Accession %in% NAprocessPCP_MAC(ptetPCP_raw_rep3)[,1]),]

### Write to CSV
write.csv(ptetPCP_raw_rep1, file = "/ptetPCP/ptetPCP_csv/ptetPCP_raw_rep1.csv", quote = F, row.names = F)
write.csv(ptetPCP_raw_rep2, file = "/ptetPCP/ptetPCP_csv/ptetPCP_raw_rep2.csv", quote = F, row.names = F)
write.csv(ptetPCP_raw_rep3, file = "/ptetPCP/ptetPCP_csv/ptetPCP_raw_rep3.csv", quote = F, row.names = F)

### Convert to MSN + Remove Unimportant R Objects
ptetPCP_raw_rep1_msn = readMSnSet2(file = "./ptetPCP_CSV/ptetPCP_raw_rep1.csv",
                                   ecol = c("MAC.1", "X300g.1", "X1K.1", "X3K.1", "X5K.1", "X9K.1", 
                                            "X12K.1", "X15K.1", "X30K.1", "X79K.1", "X120K.1", "SUP.1"), 
                                   fnames = "Accession")
ptetPCP_raw_rep2_msn = readMSnSet2(file = "./ptetPCP_CSV/ptetPCP_raw_rep2.csv",
                                   ecol = c("MAC.2", "X300g.2", "X1K.2", "X3K.2", "X5K.2", "X9K.2", 
                                            "X12K.2", "X15K.2", "X30K.2", "X79K.2", "X120K.2", "SUP.2"), 
                                   fnames = "Accession")
ptetPCP_raw_rep3_msn = readMSnSet2(file = "./ptetPCP_CSV/ptetPCP_raw_rep3.csv",
                                   ecol = c("MAC.3", "X300g.3", "X1K.3", "X3K.3", "X5K.3", "X9K.3", 
                                            "X12K.3", "X15K.3", "X30K.3", "X79K.3", "X120K.3", "SUP.3"), 
                                   fnames = "Accession")
rm(list=setdiff(ls(), c("ptetPCP_raw_rep1_msn", "ptetPCP_raw_rep2_msn", "ptetPCP_raw_rep3_msn")))


###### Data Sets

### Remove NAs and Normalize
ptetPCP_raw_rep1_norm = normalize(ptetPCP_raw_rep1_msn, method = "sum")
ptetPCP_raw_rep2_norm = normalize(ptetPCP_raw_rep2_msn, method = "sum")
ptetPCP_raw_rep3_norm = normalize(ptetPCP_raw_rep3_msn, method = "sum")

### Impute NBAVG and Normalize
ptetPCP_raw_rep1_nbavg_norm = normalize(impute(ptetPCP_raw_rep1_msn, method = "nbavg"), method = "sum")
ptetPCP_raw_rep2_nbavg_norm = normalize(impute(ptetPCP_raw_rep2_msn, method = "nbavg"), method = "sum")
ptetPCP_raw_rep3_nbavg_norm = normalize(impute(ptetPCP_raw_rep3_msn, method = "nbavg"), method = "sum")

### Impute Zeros
ptetPCP_raw_rep1_nbavg_norm_zero = impute(ptetPCP_raw_rep1_nbavg_norm, method = "zero")
ptetPCP_raw_rep2_nbavg_norm_zero = impute(ptetPCP_raw_rep2_nbavg_norm, method = "zero")
ptetPCP_raw_rep3_bavg_norm_zero = impute(ptetPCP_raw_rep3_nbavg_norm, method = "zero")

### Concatenate above datasets
# Remove NAs and Normalize
tmp1a =  merge(exprs(ptetPCP_raw_rep1_norm), exprs(ptetPCP_raw_rep2_norm), by=0)
rownames(tmp1a) = tmp1a$Row.names
tmp1a$Row.names = NULL 
tmp1b = merge(tmp1a, exprs(ptetPCP_raw_rep3_norm), by=0)
rownames(tmp1b) = tmp1b$Row.names
tmp1c = fData(ptetPCP_raw_rep1_norm)[which((fData(ptetPCP_raw_rep1_norm)$'Accession') %in% (rownames(tmp1b))),]
tmp1d = fData(ptetPCP_raw_rep2_norm)[which((fData(ptetPCP_raw_rep2_norm)$'Accession') %in% (rownames(tmp1b))),]
tmp1e = fData(ptetPCP_raw_rep3_norm)[which((fData(ptetPCP_raw_rep3_norm)$'Accession') %in% (rownames(tmp1b))),]
tmp1f = intersect(intersect(tmp1c, tmp1d), tmp1e)
ptetPCP_raw_all_norm = merge(tmp1f, tmp1b, by=0)
ptetPCP_raw_all_norm$Row.names = NULL

# Impute NBAVG and Normalize
tmp2a =  merge(exprs(ptetPCP_raw_rep1_nbavg_norm), exprs(ptetPCP_raw_rep2_nbavg_norm), by=0)
rownames(tmp2a) = tmp2a$Row.names
tmp2a$Row.names = NULL 
tmp2b = merge(tmp2a, exprs(ptetPCP_raw_rep3_nbavg_norm), by=0)
rownames(tmp2b) = tmp2b$Row.names
tmp2c = fData(ptetPCP_raw_rep1_nbavg_norm)[which((fData(ptetPCP_raw_rep1_nbavg_norm)$'Accession') %in% (rownames(tmp2b))),]
tmp2d = fData(ptetPCP_raw_rep2_nbavg_norm)[which((fData(ptetPCP_raw_rep2_nbavg_norm)$'Accession') %in% (rownames(tmp2b))),]
tmp2e = fData(ptetPCP_raw_rep3_nbavg_norm)[which((fData(ptetPCP_raw_rep3_nbavg_norm)$'Accession') %in% (rownames(tmp2b))),]
tmp2f = intersect(intersect(tmp2c, tmp2d), tmp2e)
ptetPCP_raw_all_nbavg_norm = merge(tmp2f, tmp2b, by=0)
ptetPCP_raw_all_nbavg_norm$Row.names = NULL

# Impute Zeros
tmp3a =  merge(exprs(ptetPCP_raw_rep1_nbavg_norm_zero), exprs(ptetPCP_raw_rep2_nbavg_norm_zero), by=0)
rownames(tmp3a) = tmp3a$Row.names
tmp3a$Row.names = NULL 
tmp3b = merge(tmp3a, exprs(ptetPCP_raw_rep3_bavg_norm_zero), by=0)
rownames(tmp3b) = tmp3b$Row.names
tmp3c = fData(ptetPCP_raw_rep1_nbavg_norm_zero)[which((fData(ptetPCP_raw_rep1_nbavg_norm_zero)$'Accession') %in% (rownames(tmp3b))),]
tmp3d = fData(ptetPCP_raw_rep2_nbavg_norm_zero)[which((fData(ptetPCP_raw_rep2_nbavg_norm_zero)$'Accession') %in% (rownames(tmp3b))),]
tmp3e = fData(ptetPCP_raw_rep3_bavg_norm_zero)[which((fData(ptetPCP_raw_rep3_bavg_norm_zero)$'Accession') %in% (rownames(tmp3b))),]
tmp3f = intersect(intersect(tmp3c, tmp3d), tmp3e)
ptetPCP_raw_all_bavg_norm_zero = merge(tmp3f, tmp3b, by=0)
ptetPCP_raw_all_bavg_norm_zero$Row.names = NULL

### Remove Aberrations (Flat Lines caused by Imputation Error; should have all values == 0.08333...)
ptetPCP_raw_all_bavg_norm_zero = ptetPCP_raw_all_bavg_norm_zero %>% filter(round(MAC.1, 4) != round(0.0833, 4)) %>%
  filter(round(X15K.1, 4) != round(0.0833, 4)) %>%  filter(round(SUP.1, 4) != round(0.0833, 4)) %>%
  filter(round(MAC.2, 4) != round(0.0833, 4)) %>% filter(round(X15K.2, 4) != round(0.0833, 4)) %>%
  filter(round(SUP.2, 4) != round(0.0833, 4)) %>% filter(round(MAC.3, 4) != round(0.0833, 4)) %>%
  filter(round(X15K.3, 4) != round(0.0833, 4)) %>% filter(round(SUP.3, 4) != round(0.0833, 4))

ptetPCP_raw_all_nbavg_norm = ptetPCP_raw_all_nbavg_norm %>% filter(round(MAC.1, 4) != round(0.0833, 4)) %>%
  filter(round(X15K.1, 4) != round(0.0833, 4)) %>%  filter(round(SUP.1, 4) != round(0.0833, 4)) %>%
  filter(round(MAC.2, 4) != round(0.0833, 4)) %>% filter(round(X15K.2, 4) != round(0.0833, 4)) %>%
  filter(round(SUP.2, 4) != round(0.0833, 4)) %>% filter(round(MAC.3, 4) != round(0.0833, 4)) %>%
  filter(round(X15K.3, 4) != round(0.0833, 4)) %>% filter(round(SUP.3, 4) != round(0.0833, 4))


###### Data Sets:
# Organellar Classification:    ptetPCP_raw_all_bavg_norm_zero, ptetPCP_raw_rep1_nbavg_norm_zero, ptetPCP_raw_rep2_nbavg_norm_zero, ptetPCP_raw_rep3_nbavg_norm_zero
# Nearest Neighbors:            ptetPCP_raw_all_norm, ptetPCP_raw_rep1_msn, ptetPCP_raw_rep2_msn, ptetPCP_raw_rep3_msn

### Write to new CSV
write.csv(ptetPCP_raw_all_bavg_norm_zero, file = "./ptetPCP_csv/ptetPCP_nbavg_norm_zero.csv", quote = F, row.names = F)
write.csv(ptetPCP_raw_all_nbavg_norm, file = "./ptetPCP_csv/ptetPCP_nbavg_norm.csv", quote = F, row.names = F)
write.csv(ptetPCP_raw_all_norm, file = "./ptetPCP_csv/ptetPCP_norm.csv", quote = F, row.names = F)

### Read in new MSNs
ptetPCP_raw_all_nbavg_norm_zero_NA = filterNA(readMSnSet2("./ptetPCP_csv/ptetPCP_nbavg_norm_zero.csv",
                                                          ecol = c("MAC.1", "X300g.1", "X1K.1", "X3K.1", "X5K.1", "X9K.1", 
                                                                   "X12K.1", "X15K.1", "X30K.1", "X79K.1", "X120K.1", "SUP.1", 
                                                                   "MAC.2", "X300g.2", "X1K.2", "X3K.2", "X5K.2", "X9K.2", 
                                                                   "X12K.2", "X15K.2", "X30K.2", "X79K.2", "X120K.2", "SUP.2",
                                                                   "MAC.3", "X300g.3", "X1K.3", "X3K.3", "X5K.3", "X9K.3", 
                                                                   "X12K.3", "X15K.3", "X30K.3", "X79K.3", "X120K.3", "SUP.3"),                                           
                                                          fnames = "Accession"))
ptetPCP_raw_all_norm_NA = filterNA(readMSnSet2("./ptetPCP_csv/ptetPCP_norm.csv",
                                               ecol = c("MAC.1", "X300g.1", "X1K.1", "X3K.1", "X5K.1", "X9K.1", 
                                                        "X12K.1", "X15K.1", "X30K.1", "X79K.1", "X120K.1", "SUP.1", 
                                                        "MAC.2", "X300g.2", "X1K.2", "X3K.2", "X5K.2", "X9K.2", 
                                                        "X12K.2", "X15K.2", "X30K.2", "X79K.2", "X120K.2", "SUP.2",
                                                        "MAC.3", "X300g.3", "X1K.3", "X3K.3", "X5K.3", "X9K.3", 
                                                        "X12K.3", "X15K.3", "X30K.3", "X79K.3", "X120K.3", "SUP.3"),                                           
                                               fnames = "Accession"))
ptetPCP_raw_all_nbavg_norm_NA = filterNA(readMSnSet2("./ptetPCP_csv/ptetPCP_nbavg_norm.csv",
                                                     ecol = c("MAC.1", "X300g.1", "X1K.1", "X3K.1", "X5K.1", "X9K.1", 
                                                              "X12K.1", "X15K.1", "X30K.1", "X79K.1", "X120K.1", "SUP.1", 
                                                              "MAC.2", "X300g.2", "X1K.2", "X3K.2", "X5K.2", "X9K.2", 
                                                              "X12K.2", "X15K.2", "X30K.2", "X79K.2", "X120K.2", "SUP.2",
                                                              "MAC.3", "X300g.3", "X1K.3", "X3K.3", "X5K.3", "X9K.3", 
                                                              "X12K.3", "X15K.3", "X30K.3", "X79K.3", "X120K.3", "SUP.3"),                                           
                                                     fnames = "Accession"))

### Summary
plot2D(ptetPCP_raw_all_nbavg_norm_zero_NA, fcol = NULL)  # 10349 proteins
plot2D(ptetPCP_raw_all_nbavg_norm_NA, fcol = NULL)       # 4480 proteins
plot2D(ptetPCP_raw_all_norm_NA, fcol = NULL)             # 2377 proteins

### Clean up
rm(list=setdiff(ls(), 
                c("ptetPCP_raw_rep1_msn", "ptetPCP_raw_rep2_msn", "ptetPCP_raw_rep3_msn", 
                  "ptetPCP_raw_rep1_nbavg_norm", "ptetPCP_raw_rep2_nbavg_norm", "ptetPCP_raw_rep3_nbavg_norm",
                  "ptetPCP_raw_rep1_nbavg_norm_zero", "ptetPCP_raw_rep2_nbavg_norm_zero", "ptetPCP_raw_rep3_nbavg_norm_zero",
                  "ptetPCP_raw_all_nbavg_norm_zero_NA", 
                  "ptetPCP_raw_all_norm_NA", 
                  "ptetPCP_raw_all_nbavg_norm_NA")))

#save.image("./ptetPCP_rdata/ptetPCP_rawMSN.RData")
load("./ptetPCP_rdata/ptetPCP_rawMSN.RData")


# tSNE
set.seed(42)
plot2D(ptetPCP_raw_all_norm_NA, fcol = NULL, method = "t-SNE", col = "black")
set.seed(42)
plot2D(ptetPCP_raw_all_nbavg_norm_NA, fcol = NULL, method = "t-SNE", col = "black")
set.seed(42)
plot2D(ptetPCP_raw_all_nbavg_norm_zero_NA, fcol = NULL, method = "t-SNE", col = "black")

# PCs
plot2D(ptetPCP_raw_all_norm_NA, fcol = NULL, col = "black")
plot2D(ptetPCP_raw_all_norm_NA, fcol = NULL, col = "black", dims = c(1,3))
plot2D(ptetPCP_raw_all_norm_NA, fcol = NULL, col = "black", dims = c(1,4))
plot2D(ptetPCP_raw_all_norm_NA, method = "scree")

# PCA
plot2D(ptetPCP_raw_all_norm_NA, fcol = NULL, col = "black")
plot2D(ptetPCP_raw_all_nbavg_norm_NA, fcol = NULL, col = "black")
plot2D(ptetPCP_raw_all_nbavg_norm_zero_NA, fcol = NULL, col = "black")

# venn Diagram
library(VennDiagram)
t1 = fData(ptetPCP_raw_rep1_msn)[complete.cases(exprs(ptetPCP_raw_rep1_msn)),]
t2 = fData(ptetPCP_raw_rep2_msn)[complete.cases(exprs(ptetPCP_raw_rep2_msn)),]
t3 = fData(ptetPCP_raw_rep3_msn)[complete.cases(exprs(ptetPCP_raw_rep3_msn)),]
venn.diagram(
  x = list(t1$Accession, t2$Accession, t3$Accession),
  category.names = c("Experiment 1" , "Experiment 2 " , "Experiment 3"), print.mode = c("raw", "percent"),
  filename = './ptetPCP_plots/vennDiagram_noImputation.png',
  output=TRUE, main = "Proteins Identfied in Each Experiment (Missing Values Removed)"
)
t1 = fData(ptetPCP_raw_rep1_nbavg_norm)[complete.cases(exprs(ptetPCP_raw_rep1_nbavg_norm)),]
t2 = fData(ptetPCP_raw_rep2_nbavg_norm)[complete.cases(exprs(ptetPCP_raw_rep2_nbavg_norm)),]
t3 = fData(ptetPCP_raw_rep3_nbavg_norm)[complete.cases(exprs(ptetPCP_raw_rep3_nbavg_norm)),]
venn.diagram(
  x = list(t1$Accession, t2$Accession, t3$Accession),
  category.names = c("Experiment 1" , "Experiment 2 " , "Experiment 3"), print.mode = c("raw", "percent"),
  filename = './ptetPCP_plots/vennDiagram_nbavg.png',
  output=TRUE, main = "Proteins Identfied in Each Experiment (NBAVG Imputation)"
)
t1 = fData(ptetPCP_raw_rep1_nbavg_norm_zero)[complete.cases(exprs(ptetPCP_raw_rep1_nbavg_norm_zero)),]
t2 = fData(ptetPCP_raw_rep2_nbavg_norm_zero)[complete.cases(exprs(ptetPCP_raw_rep2_nbavg_norm_zero)),]
venn.diagram(
  x = list(t1$Accession, t2$Accession),
  category.names = c("Experiment 1" , "Experiment 2 "), print.mode = c("raw", "percent"),
  filename = './ptetPCP_plots/vennDiagram_nbavg.png',
  output=TRUE, main = "Proteins Identfied in Each Experiment (ZERO Imputation)"
)
