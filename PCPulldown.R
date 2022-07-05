PCPulldown <- function(msn, protein, lDist, n_top_hits, svmProps, destination_dir){
  # Condition 1
  if(dir.exists(destination_dir) == F){
    stop(print("The directory isn't writeable (doesn't exist)"))
  }
  
  setwd(destination_dir)
  n_names = c(protein, names(sort(lDist[[protein]], decreasing = T)[1:n_top_hits]))
  n_names = n_names[which(n_names %in% rownames(exprs(msn)))]   # Check
  table(factor(n_names, levels=unique(n_names)))
  
  # Condition 2
  if(is.null(lDist[[protein]])){
    stop(print("This protein had one or more missing values. Try another lDist?"))
  }
  
  # PSS Plot
  tiff_name1 = paste(protein, "_", as.character(n_top_hits), "-nearest-neighbors_PSS.tiff", sep="")
  tiff(filename = tiff_name1)
  hist(lDist[[protein]], breaks = 500, main = protein, xlab = "Protein Similarity Score")
  dev.off()
  
  # Dist Plot
  tiff_name2 = paste(protein, "_", as.character(n_top_hits), "-nearest-neighbors_dist.tiff", sep="")
  tiff(filename = tiff_name2)
  plotDist(object = msn[n_names], fcol = NULL, pcol = c("red", rep("black", n_top_hits)))  # protein is red, rest are black
  dev.off()
  
  # output table
  outDF = svmProps[which(svmProps$Accession %in% n_names),]
  outDF = outDF[order(match(outDF$Accession,n_names)),]   # ensure order is correct
  write.csv(outDF, file = paste(protein, "_", as.character(n_top_hits), 
                                "-nearest-neighbors.csv", sep=""), row.names = F)
}
