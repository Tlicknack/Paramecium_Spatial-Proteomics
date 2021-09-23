# This will walk through the entire LOPIT Data Import, Visualization, and Classification
library("MSnbase")
library("pRoloc")
library("pRolocGUI")
#library("dplyr")
library("ggplot2")
library("vsn")

##### Read in Original csv
ptCsv = read.csv("/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/Protein-Level-csv/Ptet_LOPIT_Licknack.csv", header=T)
ptCsv$Checked = NULL
ptCsv$Master = NULL

# Make new .csv's for Downstream Analysis
#####
##### 12 Fractions
ptetLOPIT12_csv = cbind(ptCsv[,which(names(ptCsv) == "Accession"):which(names(ptCsv) == "X..Razor.Peptides")], 
                        ptCsv[,which(names(ptCsv) == "Abundances..Grouped...120K"):which(names(ptCsv) == "Abundances..Grouped...Sup")])
ptetLOPIT12_csv = ptetLOPIT12_csv[c("Accession", "Description", "Exp..q.value..Combined", "Sum.PEP.Score", "Coverage....", "X..Peptides",
                                    "X..PSMs", "X..Unique.Peptides", "X..AAs", "MW..kDa.", "calc..pI", "Score.Sequest.HT..Sequest.HT", "X..Peptides..by.Search.Engine...Sequest.HT", "X..Razor.Peptides",
                                    "Abundances..Grouped...Mac", "Abundances..Grouped...300g", "Abundances..Grouped...1K", 
                                    "Abundances..Grouped...3K", "Abundances..Grouped...5K", "Abundances..Grouped...9K", "Abundances..Grouped...12K",
                                    "Abundances..Grouped...15K", "Abundances..Grouped...30K", "Abundances..Grouped...79K", 
                                    "Abundances..Grouped...120K", "Abundances..Grouped...Sup")]
      # extract columns of interest
colnames(ptetLOPIT12_csv) = c("Accession", "Description", "q.value", "Sum-PEP-Score", "Coverage", "Number-Peptides", "Number-PSMs", 
                              "Number-Unique-Peptides", "Number-AA", "MW", "pI", "Score-Sequest", "Number-Peptides-Sequest-HT", 
                              "Number-Razor-Peptides", "MAC", "300g", "1K", "3K", "5K", "9K", "12K","15K", "30K", "79K", "120K", "Sup")
      # Put Fractions in order
write.csv(x = ptetLOPIT_final, file = "/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/Protein-Level-csv/ptet_LOPIT12.csv", row.names = F)
#####
#####
#####
##### 36 Fractions 
ptetLOPIT36_csv = cbind(ptCsv[,which(names(ptCsv) == "Accession"):which(names(ptCsv) == "X..Razor.Peptides")], 
                        ptCsv[,which(names(ptCsv) == "Abundances..Normalized...F74..Sample..120K"):which(names(ptCsv) == "Abundances..Normalized...F109..Sample..Sup")])
# Average each replicate
ptetLOPIT36_csv_abundances = ptetLOPIT36_csv[,15:length(names(ptetLOPIT36_csv))]  # get the columns with sum normalized data
ptetLOPIT36_csv_abundances120_one = rowMeans(ptetLOPIT36_csv_abundances[,1:3])    # 120K
ptetLOPIT36_csv_abundances120_two = rowMeans(ptetLOPIT36_csv_abundances[,4:6])
ptetLOPIT36_csv_abundances120_three = rowMeans(ptetLOPIT36_csv_abundances[,7:9])
ptetLOPIT36_csv_abundances12_one = rowMeans(ptetLOPIT36_csv_abundances[,10:12])   # 12K
ptetLOPIT36_csv_abundances12_two = rowMeans(ptetLOPIT36_csv_abundances[,13:15])
ptetLOPIT36_csv_abundances12_three = rowMeans(ptetLOPIT36_csv_abundances[,16:18])
ptetLOPIT36_csv_abundances15_one = rowMeans(ptetLOPIT36_csv_abundances[,19:21])   # 15K
ptetLOPIT36_csv_abundances15_two = rowMeans(ptetLOPIT36_csv_abundances[,22:24])   
ptetLOPIT36_csv_abundances15_three = rowMeans(ptetLOPIT36_csv_abundances[,25:27])
ptetLOPIT36_csv_abundances1_one = rowMeans(ptetLOPIT36_csv_abundances[,28:30])    # 1K
ptetLOPIT36_csv_abundances1_two = rowMeans(ptetLOPIT36_csv_abundances[,31:33])
ptetLOPIT36_csv_abundances1_three = rowMeans(ptetLOPIT36_csv_abundances[,34:36])
ptetLOPIT36_csv_abundances300_one = rowMeans(ptetLOPIT36_csv_abundances[,37:39])  # 300g
ptetLOPIT36_csv_abundances300_two = rowMeans(ptetLOPIT36_csv_abundances[,40:42])
ptetLOPIT36_csv_abundances300_three = rowMeans(ptetLOPIT36_csv_abundances[,43:45])
ptetLOPIT36_csv_abundances30_one = rowMeans(ptetLOPIT36_csv_abundances[,46:48])   # 30K
ptetLOPIT36_csv_abundances30_two = rowMeans(ptetLOPIT36_csv_abundances[,49:51])
ptetLOPIT36_csv_abundances30_three = rowMeans(ptetLOPIT36_csv_abundances[,52:54])
ptetLOPIT36_csv_abundances3_one = rowMeans(ptetLOPIT36_csv_abundances[,55:57])    # 3K
ptetLOPIT36_csv_abundances3_two = rowMeans(ptetLOPIT36_csv_abundances[,58:60])
ptetLOPIT36_csv_abundances3_three = rowMeans(ptetLOPIT36_csv_abundances[,61:63])
ptetLOPIT36_csv_abundances5_one = rowMeans(ptetLOPIT36_csv_abundances[,64:66])    # 5K
ptetLOPIT36_csv_abundances5_two = rowMeans(ptetLOPIT36_csv_abundances[,67:69])
ptetLOPIT36_csv_abundances5_three = rowMeans(ptetLOPIT36_csv_abundances[,70:72])
ptetLOPIT36_csv_abundances79_one = rowMeans(ptetLOPIT36_csv_abundances[,73:75])   # 79K
ptetLOPIT36_csv_abundances79_two = rowMeans(ptetLOPIT36_csv_abundances[,76:78])
ptetLOPIT36_csv_abundances79_three = rowMeans(ptetLOPIT36_csv_abundances[,79:81])
ptetLOPIT36_csv_abundances9_one = rowMeans(ptetLOPIT36_csv_abundances[,82:84])    # 9K
ptetLOPIT36_csv_abundances9_two = rowMeans(ptetLOPIT36_csv_abundances[,85:87])
ptetLOPIT36_csv_abundances9_three = rowMeans(ptetLOPIT36_csv_abundances[,88:90])
ptetLOPIT36_csv_abundancesMAC_one = rowMeans(ptetLOPIT36_csv_abundances[,91:93])  # MAC
ptetLOPIT36_csv_abundancesMAC_two = rowMeans(ptetLOPIT36_csv_abundances[,94:96])
ptetLOPIT36_csv_abundancesMAC_three = rowMeans(ptetLOPIT36_csv_abundances[,97:99])
ptetLOPIT36_csv_abundancesSUP_one = rowMeans(ptetLOPIT36_csv_abundances[,100:102])# SUP
ptetLOPIT36_csv_abundancesSUP_two = rowMeans(ptetLOPIT36_csv_abundances[,103:105])
ptetLOPIT36_csv_abundancesSUP_three = rowMeans(ptetLOPIT36_csv_abundances[,106:108])

ptetLOPIT36_csv_means = data.frame(  # combine averages
  ptetLOPIT36_csv_abundancesMAC_one, ptetLOPIT36_csv_abundances300_one, ptetLOPIT36_csv_abundances1_one, 
  ptetLOPIT36_csv_abundances3_one, ptetLOPIT36_csv_abundances5_one, ptetLOPIT36_csv_abundances9_one, 
  ptetLOPIT36_csv_abundances12_one, ptetLOPIT36_csv_abundances15_one, ptetLOPIT36_csv_abundances30_one, 
  ptetLOPIT36_csv_abundances79_one, ptetLOPIT36_csv_abundances120_one,ptetLOPIT36_csv_abundancesSUP_one,
  ptetLOPIT36_csv_abundancesMAC_two, ptetLOPIT36_csv_abundances300_two, ptetLOPIT36_csv_abundances1_two, 
  ptetLOPIT36_csv_abundances3_two, ptetLOPIT36_csv_abundances5_two, ptetLOPIT36_csv_abundances9_two, 
  ptetLOPIT36_csv_abundances12_two, ptetLOPIT36_csv_abundances15_two, ptetLOPIT36_csv_abundances30_two, 
  ptetLOPIT36_csv_abundances79_two, ptetLOPIT36_csv_abundances120_two,ptetLOPIT36_csv_abundancesSUP_two,
  ptetLOPIT36_csv_abundancesMAC_three, ptetLOPIT36_csv_abundances300_three, ptetLOPIT36_csv_abundances1_three, 
  ptetLOPIT36_csv_abundances3_three, ptetLOPIT36_csv_abundances5_three, ptetLOPIT36_csv_abundances9_three, 
  ptetLOPIT36_csv_abundances12_three, ptetLOPIT36_csv_abundances15_three, ptetLOPIT36_csv_abundances30_three, 
  ptetLOPIT36_csv_abundances79_three, ptetLOPIT36_csv_abundances120_three,ptetLOPIT36_csv_abundancesSUP_three)
ptetLOPIT36_csv_means = cbind(ptCsv[,which(names(ptCsv) == "Accession")], ptetLOPIT36_csv_means)

colnames(ptetLOPIT36_csv_means) = c("Accession", "MAC-1", "300g-1", "1kg-1", "3kg-1", "5kg-1", "9kg-1", "12kg-1", "15kg-1", "30kg-1", "79kg-1", "120kg-1", "Sup-1",
                                     "MAC-2", "300g-2", "1kg-2", "3kg-2", "5kg-2", "9kg-2", "12kg-2", "15kg-2", "30kg-2", "79kg-2", "120kg-2", "Sup-2",
                                     "MAC-3", "300g-3", "1kg-3", "3kg-3", "5kg-3", "9kg-3", "12kg-3", "15kg-3", "30kg-3", "79kg-3", "120kg-3", "Sup-3")

write.csv(ptetLOPIT36_csv_means, "/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/Protein-Level-csv/ptet_LOPIT_36.csv", row.names = F)

#####
#####
##### Read in MSnSet Object
ptetLOPIT12_msn = readMSnSet2(file = "/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/Protein-Level-csv/ptet_LOPIT12.csv", 
                        ecol = c("MAC", "X300g", "X1K", "X3K", "X5K", "X9K", "X12K","X15K", "X30K", "X79K", "X120K", "Sup"), 
                        fnames = "Accession")

ptetLOPIT36_msn = readMSnSet2(file = "/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/Protein-Level-csv/ptet_LOPIT_36.csv", 
                                         ecol = c("MAC.1", "X300g.1", "X1kg.1", "X3kg.1", "X5kg.1", "X9kg.1", "X12kg.1","X15kg.1", "X30kg.1", "X79kg.1", "X120kg.1", "Sup.1",
                                                  "MAC.2", "X300g.2", "X1kg.2", "X3kg.2", "X5kg.2", "X9kg.2", "X12kg.2","X15kg.2", "X30kg.2", "X79kg.2", "X120kg.2", "Sup.2",
                                                  "MAC.3", "X300g.3", "X1kg.3", "X3kg.3", "X5kg.3", "X9kg.3", "X12kg.3","X15kg.3", "X30kg.3", "X79kg.3", "X120kg.3", "Sup.3"), 
                                         fnames = "Accession")

save.image("/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/RData-LOPIT/ptetLOPIT_12-36.RData")

#####
#####
##### Structure of Data
plot2D(object = ptetLOPIT12_msn, fcol = NULL)  
plot2D(object = ptetLOPIT36_msn, fcol = NULL)  

#####
#####
##### Examine/Handle Missing Data
# 12 Fraction
# Impute Missing Data
ptetLOPIT12_knnImpute = impute(object = ptetLOPIT12_msn, method = "knn")
ptetLOPIT12_zeroImpute = impute(object = ptetLOPIT12_msn, method = "zero")   
head(exprs(ptetLOPIT12_msn))
head(exprs(ptetLOPIT12_knnImpute))
head(exprs(ptetLOPIT12_zeroImpute))
plot2D(object = ptetLOPIT12_msn, fcol = NULL)           # PC1 = 33.89%
plot2D(object = ptetLOPIT12_knnImpute, fcol = NULL)     # PC1 = 29.1%
plot2D(object = ptetLOPIT12_zeroImpute, fcol = NULL)    # PC1 = 29.32%
# Remove N/As
ptetLOPIT12_NA = filterNA(ptetLOPIT12_msn)
head(exprs(ptetLOPIT12_NA))
# Noramlize- Proceeding with NA's Removed
ptetLOPIT12_NA_sumNormalized = normalise(object = ptetLOPIT12_NA, method = "sum")  # normalize to sum
ptetLOPIT12_NA_vscNormalized = normalise(object = ptetLOPIT12_NA, method = "vsn")  # normalize to variance
meanSdPlot(ptetLOPIT12_NA)                # plot of standard deviation vs mean... no relationship
meanSdPlot(ptetLOPIT12_NA_sumNormalized)  # plot of standard deviation vs mean... no relationship, but much lower
meanSdPlot(ptetLOPIT12_NA_vscNormalized)  # plot of standard deviation vs mean... sd drops as mean increases, but still higher than sum
    # Sum normalization reduces SD far more
    # ptetLOPIT12_NA_sumNormalized is the go-to dataset     

# For 36 Fraction Dataset
ptetLOPIT36_knnImpute = impute(object = ptetLOPIT36_msn, method = "knn")
ptetLOPIT36_zeroImpute = impute(object = ptetLOPIT36_msn, method = "zero")   
plot2D(object = ptetLOPIT36_msn, fcol = NULL)           # PC1 = 60.63% ... PCA is a terrible shape
plot2D(object = ptetLOPIT36_knnImpute, fcol = NULL)     # PC1 = 61.28%
plot2D(object = ptetLOPIT36_zeroImpute, fcol = NULL)    # PC1 = 60.66%
# Remove N/As
ptetLOPIT36_NA = filterNA(ptetLOPIT36_msn)
head(exprs(ptetLOPIT36_NA))
# Noramlize- Proceeding with NA's Removed
ptetLOPIT36_NA_sumNormalized = normalise(object = ptetLOPIT36_NA, method = "sum")  # normalize to sum
ptetLOPIT36_NA_vscNormalized = normalise(object = ptetLOPIT36_NA, method = "vsn")  # normalize to variance
meanSdPlot(ptetLOPIT36_NA)                # plot of standard deviation vs mean... exponential increase in sd as mean increases
meanSdPlot(ptetLOPIT36_NA_sumNormalized)  # plot of standard deviation vs mean... no relationship, but much lower
meanSdPlot(ptetLOPIT36_NA_vscNormalized)  # plot of standard deviation vs mean... weak positive relationship, but higher than sum
    # Sum normalization reduces SD far more
    # ptetLOPIT36_NA_sumNormalized is the go-to dataset     
save.image("/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/RData-LOPIT/ptetLOPIT_12-36_impute_normalize.RData")

#####
#####
##### Loading Markers 
########## See LOPIT-MarkerProteins.R for Details










