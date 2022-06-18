  # ptetPCP_pcpulldown.R

setwd('/ptetPCP/')

library(pRoloc)

load("./ptetPCP_rdata/ProteinSimilarityScores.RData")

ptetProps = read.csv("./ptetPCP_csv/ptetPCP_Properties.csv", header= T)
yeastProps = read.csv("./ptetPCP_csv/yeast2018_properties.csv", header=T)
toxoProps = read.csv("./ptetPCP_csv/toxo2020_properties.csv", header=T)

source("/ptetPCP/ptetPCP_scripts/PCPulldown.R")

#
# Ptet
  # CAM
PCPulldown(msn = preds2, protein = "PTET.51.1.P0460139", 
            lDist = lDist_NA, n_top_hits = 20, svmProps = ptetProps, 
           destination_dir = "/ptetPCP/ptetPCP_eucDist/PCPulldowns/")
  # rab_A69
PCPulldown(msn = preds2, protein = "PTET.51.1.P0020101", 
           lDist = lDist_NA, n_top_hits = 20, svmProps = ptetProps, 
           destination_dir = "/ptetPCP/ptetPCP_eucDist/PCPulldowns/")
  # Tom40
PCPulldown(msn = preds2, protein = "PTET.51.1.P0280026", 
           lDist = lDist_NA, n_top_hits = 20, svmProps = ptetProps, 
           destination_dir = "/ptetPCP/ptetPCP_eucDist/PCPulldowns/")

#
# Yeast- ~43% of the proteome covered; ~29 proteins would be 99th percentile
  # Tom40
PCPulldown(msn = yeast2018, protein = "P23644", 
           lDist = lDist_yeast, n_top_hits = 50, svmProps = yeastProps, 
           destination_dir = "/ptetPCP/ptetPCP_eucDist/PCPulldowns/") 
  # ASC1: https://www.mcponline.org/article/S1535-9476(20)32315-X/fulltext 
PCPulldown(msn = yeast2018, protein = "P38011", 
           lDist = lDist_yeast, n_top_hits = 29, svmProps = yeastProps, 
           destination_dir = "/ptetPCP/ptetPCP_eucDist/PCPulldowns/") 
  # ATP1: 
PCPulldown(msn = yeast2018, protein = "P07251", 
           lDist = lDist_yeast, n_top_hits = 5, svmProps = yeastProps, 
           destination_dir = "/ptetPCP/ptetPCP_eucDist/PCPulldowns/") 
  # RPA49: 
PCPulldown(msn = yeast2018, protein = "Q01080", 
           lDist = lDist_yeast, n_top_hits = 5, svmProps = yeastProps, 
           destination_dir = "/ptetPCP/ptetPCP_eucDist/PCPulldowns/") 
# FBP26: 
PCPulldown(msn = yeast2018, protein = "P32604", 
           lDist = lDist_yeast, n_top_hits = 5, svmProps = yeastProps, 
           destination_dir = "/ptetPCP/ptetPCP_eucDist/PCPulldowns/") 

#
# Toxo
  # MORN3 https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001081 
PCPulldown(msn = Barylyuk2020ToxoLopit, protein = "TGME49_245710", 
           lDist = lDist_toxo, n_top_hits = 20, svmProps = toxoProps, 
           destination_dir = "/ptetPCP/ptetPCP_eucDist/PCPulldowns/") 
  # TOM40
PCPulldown(msn = Barylyuk2020ToxoLopit, protein = "TGME49_218280", 
           lDist = lDist_toxo, n_top_hits = 20, svmProps = toxoProps, 
           destination_dir = "/ptetPCP/ptetPCP_eucDist/PCPulldowns/") 


