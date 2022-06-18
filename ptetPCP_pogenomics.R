  # ptetPCP_pogenomics.R

library(data.table)

ptetGFF = read.table("/Paramecium_Annotations/ptet-gene.tab", header=F)
ptetgncGFF = read.table("/ptetPCP/ptetPCP_popgenomics/tetraurelia51_EuGene_annotation_edit_intron_intergenic_oriented.gff", header=F)
ptetgncGFF = ptetgncGFF[grep("gene", ptetgncGFF$V3),]
ptetgncGFF$V1 = gsub("scaffold_", "scaffold51_", ptetgncGFF$V1)
ptetgncGFF$V2 = NULL
ptetgncGFF$V3 = NULL
ptetgncGFF$V6 = NULL
ptetgncGFF$V8 = NULL
ptetgncGFF$V9 = as.character(sapply(X = unlist(lapply(strsplit(ptetgncGFF$V9, split = ";"), "[", 1)), FUN = substr, start=4, stop=16))

colnames(ptetGFF) = c("scf", "start", "end", "strand", "Accession")
colnames(ptetgncGFF) = c("scf", "start", "end", "strand", "Accession")


ptetConversion = data.frame(scaffold = NA, GeneA = NA, GeneAcoords1 = NA, GeneAcoords2 = NA, 
                            GeneB = NA , GeneBcoords1 = NA, GeneBcoords2 = NA)

for(scaf in unique(ptetGFF$scf)){
  ptGFF1 = data.table(
    scaffold = scaf,
    GeneA = ptetGFF[which(ptetGFF$scf == scaf),5],
    GeneAcoords1 =     ptetGFF[which(ptetGFF$scf == scaf),2],
    GeneAcoords2 =     ptetGFF[which(ptetGFF$scf == scaf),3])
  setkey(ptGFF1, scaffold, GeneA, GeneAcoords1, GeneAcoords2)
  
  ptGFF2 = data.table(
    scaffold = scaf,
    GeneB = ptetgncGFF[which(ptetgncGFF$scf == scaf),5],
    GeneBcoords1 =     ptetgncGFF[which(ptetgncGFF$scf == scaf),2],
    GeneBcoords2 =     ptetgncGFF[which(ptetgncGFF$scf == scaf),3])  
  setkey(ptGFF2, scaffold, GeneB, GeneBcoords1, GeneBcoords2)
  
  for(q in 1:nrow(ptGFF1)){
    ptetConversion = rbind(ptetConversion, merge(ptGFF1[q], ptGFF2[which((ptGFF2$GeneBcoords1 <= ptGFF1$GeneAcoords1[q]) & (ptGFF2$GeneBcoords2 >= ptGFF1$GeneAcoords1[q])),], 
          by="scaffold"))
  }
}
ptetConversion$GeneB = sub("PTETEGNG", "PTETEGNC", ptetConversion$GeneB)  #match Parul's format
write.csv(ptetConversion, "/ptetPCP/ptetPCP_popgenomics/ptetIDconversion.csv", quote = F, row.names = F)

# merge with pi/ps table from Parul 
load("./ptetPCP_rdata/ptetPCP_SVM.Rdata")
library(ggplot2)

ptPPI = read.table("/ptetPCP/ptetPCP_popgenomics/all_genes_filtered_tetraurelia.char", header=T)

ptPopGen = merge(merge(ptetConversion, ptPPI, by.x="GeneB", by.y="gene"), 
                 fData(preds2), by.x="GeneA", by.y="Accession")
ptPopGen$svm.pred = gsub(pattern = " ", replacement = "\n", ptPopGen$svm.pred)
ptPopGen$svm.pred = gsub(pattern = "-", replacement = "\n", ptPopGen$svm.pred)

ggplot(ptPopGen) + geom_boxplot(aes(x=svm.pred, y=piNpiS)) + ylim(c(0,1)) + 
  xlab("Predicted Compartment") + ylab("piN / piS")
TukeyHSD(aov(piNpiS ~ svm.pred, data = ptPopGen))  # none

#save.image("./ptetPCP_rdata/ptetPCP_popGenome.Rdata")
load("./ptetPCP_rdata/ptetPCP_popGenome.Rdata")

