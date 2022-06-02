  # ptetPCP_makeEcuDistDB.R

setwd("/ptetPCP/")
load("./ptetPCP_rdata/ptetPCP_rawMSN.RData")
euclidean <- function(a, b) sqrt(sum((a - b)^2))

filterNA(ptetPCP_raw_all_norm_NA)
filterNA(ptetPCP_raw_all_nbavg_norm_NA)
filterNA(ptetPCP_raw_all_nbavg_norm_zero_NA)

noImpute = as.data.frame(exprs(ptetPCP_raw_all_norm_NA))
noImpute = noImpute[grep("PTET*", rownames(noImpute)),]

write.csv(noImpute, "./ptetPCP_csv/ptetPCP_NA_EXPRS.csv", quote = F, row.names = T)

nbavgImpute = as.data.frame(exprs(ptetPCP_raw_all_nbavg_norm_NA))
nbavgImpute = nbavgImpute[grep("PTET*", rownames(nbavgImpute)),]

write.csv(nbavgImpute, "./ptetPCP_csv/ptetPCP_NBAVG_EXPRS.csv", quote = F, row.names = T)

lDist_NA = list()
for(i in 1:nrow(noImpute)){
  for(j in 1:nrow(noImpute)){
    chi_ij = euclidean(noImpute[i,2:ncol(noImpute)], noImpute[j,2:ncol(noImpute)])
    
    if(chi_ij > 0){
      lDist_NA[[noImpute[i,1]]][noImpute[j,1]] = chi_ij
    }
  }
  save.image("./ptetPCP_eucDist.RData")
}

lDist_NBAVG = list()
for(i in 1:nrow(nbavgImpute)){
  for(j in 1:nrow(nbavgImpute)){
    chi_ij = euclidean(nbavgImpute[i,2:ncol(nbavgImpute)], nbavgImpute[j,2:ncol(nbavgImpute)])
    
    if(chi_ij > 0){
      lDist_NBAVG[[nbavgImpute[i,1]]][nbavgImpute[j,1]] = chi_ij
    }
  }
  save.image("./ptetPCP_eucDist.RData")
}


# Yeast
library(pRolocdata)
data("yeast2018")
yeast_exprs = 2^(exprs(yeast2018))
write.table(rownames(fData(yeast2018)), file = "/ptetPCP/ptetPCP_csv/yeast2018_genes.tab", 
            sep = ";", row.names = F, quote = F, col.names = F)
    # uploaded this to YeastMine to get Descriptions. wrote that result to yeast2018_properties.csv
yeastProps = read.csv("/ptetPCP/ptetPCP_csv/yeast2018_properties.csv", header=T)
write.csv(merge(yeastProps, fData(yeast2018), by.x="Accession", by.y=0), 
          file = "/ptetPCP/ptetPCP_csv/yeast2018_properties.csv", quote = F, row.names = F)

for(i in 1:nrow(yeast_exprs)){   # make DF
  for(j in 1:nrow(yeast_exprs)){
    chi_ij = euclidean(yeast_exprs[i,1:ncol(yeast_exprs)], yeast_exprs[j,1:ncol(yeast_exprs)])
    
    if(chi_ij > 0){
      lDist_yeast[[rownames(yeast_exprs)[i]]] [rownames(yeast_exprs)[j]] = chi_ij
    }
  }
  save.image("./yeast2018_eucDist.RData")
}
