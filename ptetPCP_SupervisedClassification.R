  # ptetPCP_SupervisedClassification.R

setwd("/ptetPCP/")
load("./ptetPCP_rdata/ptetPCP_knn17.RData")  # KNN includes markers

# Create class weights
w = table(getMarkers(marked_knn17, verbose = TRUE))
w <- 1/w[names(w) != "unknown"]
svmMarked  = svmOptimisation(filterNA(marked_knn17), fcol = "markers", 
                             times = 100, xval = 5, classWeights = w, verbose = FALSE)                        
svmresMarked = svmClassification(filterNA(marked_knn17), fcol = "markers", assessRes = svmMarked)  
processingData(svmresMarked)  
preds1 = getPredictions(svmresMarked, fcol = "svm", mcol = "markers")   
medSVM = median(fData(svmresMarked)$svm.scores)  #median svm score: .53
preds2 = getPredictions(svmresMarked, fcol = "svm", t = medSVM, mcol = "markers") 
ts1 = tapply(fData(preds1)$'svm.scores', fData(preds1)$'svm', median)
preds3 = getPredictions(svmresMarked, fcol = "svm", t = ts1, mcol = "markers") 
  
#save.image("./ptetPCP_rdata/ptetPCP_SVM.Rdata")
load("./ptetPCP_rdata/ptetPCP_SVM.Rdata")
write.csv(fData(preds2), "./ptetPCP_csv/SVMpredictions.csv", row.names = F, quote = F)
  # Plugged the gene list into Biomart for properties. Manually assembled into PCP_properties.csv
write.csv(data.frame(ts1), file = "/ptetPCP/ptetPCP_csv/ptetPCP_SVMscores_median.csv", quote = F, row.names = F)

#
# tSNE
set.seed(42)
plot2D(preds2, fcol = "svm", method = "t-SNE")
addLegend(preds2, fcol = "svm", where = "topright", bty = "n", cex = .6)
set.seed(42)
plot2D(preds2, fcol = "svm.pred", method = "t-SNE")
addLegend(preds2, fcol = "svm.pred", where = "topright", bty = "n", cex = .6)
set.seed(42)
plot2D(preds3, fcol = "svm.pred", method = "t-SNE")
addLegend(preds3, fcol = "svm.pred", where = "topright", bty = "n", cex = .6)

# Summary Table
library(gt)
library(dplyr)
tmp1 = merge(data.frame(table(fData(preds2)$'markers')), data.frame(table(fData(preds2)$'svm')), by = "Var1")
preds2_compare = merge(tmp1, data.frame(table(fData(preds2)$'svm.pred')))
colnames(preds2_compare) = c("Organellar Compartment", "# Markers", "# Predicted", "# Predicted (post-filtering)")
preds2_compare$Description = c("Core Basal Body", "Ciliary- Basal Body", "Cytosolic Enzymes; Not Membrane-bound Organelles", 
                               "Phagosomal/Lysosomal Enzymes", "ER Chaperones", "Contractile Vacuole; Vesicular Transport", 
                               "Cytoskeletal; Vesicular Transport", "Mitochondrial Inner Membrane/Matrix", 
                               "Mitochondrial Inner Membrane/Matrix", "Mitochondrial Matrix only", "Mitochondrial Outer Membrane", 
                               "Macro/Micronuclear", "Surface Antigens", "Peroxisomal Enzymes/Import", 
                               "20S and 26S Proteasome", "40S and 60S Ribosome", "Trichocyst Matrix proteins (TMPs)")
preds2_compare = preds2_compare %>% 
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))
preds2_compare$Description[18] = "."
preds2_compare$SVMmed = as.numeric(ts1)
write.csv(preds2_compare, "./ptetPCP_csv/SVMpredictions_Summary.csv", row.names = F, quote = F)

preds2_compare %>% gt() %>% cols_align(align = "center") 


# Knn vs SVM
library(ggplot2)
knnsvm = fData(preds2)
knnsvm$Cluster = as.character(knnsvm$Cluster)
knnsvm$Cluster = factor(knnsvm$Cluster, levels=c("6", "5", "1", "10", "2", "13", "16", "8", "17", "4", "11", "7", "15", "3", "9", "14", "12"))

ggplot(knnsvm) + geom_count(aes(x=Cluster, y=svm.pred, color = ..n.., size = ..n..)) + 
  guides(color="legend") +
  xlab("De Novo Cluster") + ylab("SVM Predicted Compartment")


# Klebsiella?
table(fData(preds2[grep(pattern = "YP_*", x = fData(preds2)$'Accession')])$'svm.pred')
plotDist(preds2[grep(pattern = "YP_*", x = fData(preds2)$'Accession')], fcol = "svm.pred")


