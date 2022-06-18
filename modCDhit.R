  # modCDhit.R
library(seqinr)
library(stringr)

modCDhitClust = function(x){
  x = unlist(x)
  xx = paste(unlist(str_extract_all(x, "PTET.51.1.P[0-9][0-9][0-9][0-9][0-9][0-9][0-9]")), collapse = "; ")
  return(xx)
}

cdHit = read.table("/ptetPCP/ptetPCP_csv/ptet-reduced.protein_99percent_identicalLen.fa.clstr", 
                   header=F, as.is = T, fill = T)

cdHit = read.fasta("/ptetPCP/ptetPCP_csv/ptet-reduced.protein_99percent_identicalLen.fa.clstr", as.string = T, seqtype = "AA") 
newDF = data.frame(matrix(ncol=3, nrow=length(cdHit)))
colnames(newDF) = c("Cluster", "NumberOfGenes", "Accession(s)")
newDF$Cluster = names(cdHit)  
test = unlist(lapply(getSequence(cdHit, as.string = T), modCDhitClust))
newDF$`Accession(s)` = test
newDF$NumberOfGenes = as.numeric(sapply(newDF$`Accession(s)`, str_count, ";"))+1

write.csv(newDF, file = "/ptetPCP/ptetPCP_csv/CDHIT_clusters.csv", quote=F, row.names = F)
