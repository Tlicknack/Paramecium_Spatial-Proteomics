  # ptetPCP_DeNovoClustering.R

setwd("/ptetPCP/")
save.image("./ptetPCP_rdata/ptetPCP_marked.RData")

library(pRoloc)

# K-nearest Neighbors
marked_knn17 = MLearn( ~ ., filterNA(marked), kmeansI, centers=17)

marked_knn17_df = as.data.frame(unlist(marked_knn17@RObject$cluster))

colnames(marked_knn17_df) = "Cluster"
write.csv(marked_knn17_df, "./ptetPCP_csv/ptetPCP_knn17.csv", quote = F, row.names = T) 

marked_knn17 = addMarkers(marked, markers = "./ptetPCP_csv/ptetPCP_knn17.csv", 
                          mcol = "Cluster", fcol = "Accession", verbose = T)

#save.image("./ptetPCP_rdata/ptetPCP_knn17.RData")
load("./ptetPCP_rdata/ptetPCP_knn17.RData")


# t-SNE with KNN Clusters
library("randomcoloR")
set.seed(42)
knnColors = distinctColorPalette(17)
set.seed(42)
plot2D(marked_knn17, fcol = "Cluster", method = "t-SNE", col = knnColors)
addLegend(marked_knn17, fcol = "Cluster", where = "topright", bty = "n", cex = .7)
