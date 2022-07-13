  # ptetPCP_eucDist_plots.R
setwd("/ptetPCP/")
load("./ptetPCP_rdata/ptetPCP_SVM.Rdata")
load("./ptetPCP_rdata/ptetPCP_eucDist.RData")

library(tidyverse)

preds2_fd = fData(preds2)[grep("PTET*", fData(preds2)$'Accession'),]

# Range of Euclidean Distances
range(unlist(lDist_NA))

minEuc = range(unlist(lDist_NA))[1]
maxEuc = range(unlist(lDist_NA))[2]

# Convert to Similarity Score
  # Normalize Distance, subtract by 1 to make Protein Similarity Score
normalize <- function(x) {return (1 - ((x - minEuc) / (maxEuc - minEuc))) }
for(q in 1:length(lDist_NA)){
  lDist_NA[[names(lDist_NA)[q]]] = sapply(lDist_NA[[q]], FUN = normalize)
}
range(unlist(lDist_NA))

# Make Plots
setwd("/ptetPCP/ptetPCP_eucDist/NA-norm//")
for(compartment in unique(preds2_fd$svm.pred)){
  dir.create(compartment)
  setwd(compartment)
  genesC = preds2_fd[which(preds2_fd$svm.pred == compartment),'Accession']
  
  for(gene in genesC){
    if(length(lDist_NA[[gene]]) > 0){
      png(paste(gene, "_hist.png", sep = ""))  # saves histogram with this name and next plot
      hist(lDist_NA[[gene]], breaks = 500, main = as.character(gene), xlab = ("Protein Similarity Score"))
      dev.off()  # shut off plot saving
    }
  }
  setwd("../")  # go back to parent dir
}

# Create an "average" distribution for each compartment
  # Sort so protein distributions are in order
lDist_NA = lapply(lDist_NA, sort)

setwd("/ptetPCP/ptetPCP_eucDist/NA-norm/")
for(compartment in unique(preds2_fd$svm.pred)){  # make plots
  genesC = preds2_fd[which(preds2_fd$svm.pred == compartment),'Accession']
  png(paste(compartment, "_hist.png", sep = ""))
  hist(rowMeans(data.frame(lDist_NA[which(names(lDist_NA) %in% genesC)])), breaks = 100, 
      main = as.character(compartment), xlab = "Protein Similarity Score", col = group.colors[compartment])
  # normalize each part of list containing each organelle's proteins. rowmeans() does all the work after binding together a DF
  dev.off()
}

  # Statistics-- normality? bi/trimodality? Different from unknown? mixed guassian?
library("LaplacesDemon")
library("mixtools")
lCompartmentDFs = list()
lCompartmentShapiroWilk = list()
lCompartmentBiModal = list()
lCompartmentTriModal = list()
lCompartmentMixedG = list()

for(compartment in unique(preds2_fd$svm.pred)){  # compare statistically
  genesC = preds2_fd[which(preds2_fd$svm.pred == compartment),'Accession']
  dfC = data.frame(rowMeans(data.frame(lDist_NA[which(names(lDist_NA) %in% genesC)])))
  lCompartmentDFs[[compartment]] = dfC
  lCompartmentShapiroWilk[[compartment]] = shapiro.test(dfC$rowMeans.data.frame.lDist_NA.which.names.lDist_NA...in..genesC....)
  lCompartmentBiModal[[compartment]] = is.bimodal(dfC$rowMeans.data.frame.lDist_NA.which.names.lDist_NA...in..genesC....)
  lCompartmentTriModal[[compartment]] = is.trimodal(dfC$rowMeans.data.frame.lDist_NA.which.names.lDist_NA...in..genesC....)
  lCompartmentMixedG[[compartment]] = normalmixEM(dfC$rowMeans.data.frame.lDist_NA.which.names.lDist_NA...in..genesC...., k=2)
}
lCompartmentShapiroWilk  # only pellicle is normal
lCompartmentBiModal  # Mitochondria-1, Basal Body 1 are bimodal
lCompartmentTriModal  # none
lCompartmentMixedG

  # make plots: comparing all to unknown
setwd("/ptetPCP/ptetPCP_plots/")
for(compartment in unique(preds2_fd$svm.pred)){ 
  png(paste("scatterplot_lDist_", compartment, "vsUnknown.png", sep = "")) 
  plot(as.numeric(unlist(get(compartment, lCompartmentDFs))), as.numeric(unlist(lCompartmentDFs$unknown)), 
       xlab = "Uknown Distribution of PSSs", ylab = compartment)
  dev.off()  # shut off plot saving
}

  # make plots: Mixed guassian
setwd("/ptetPCP/ptetPCP_eucDist/NA-norm/")
for(compartment in unique(preds2_fd$svm.pred)){ 
  png(paste("scatterplot_lDist_", compartment, "vsUnknown.png", sep = "")) 
  plot(lCompartmentMixedG[[compartment]], breaks=100, which=2, xlab2="Protein Similarity Score",
       main2=compartment)
  dev.off()  # shut off plot saving
}

# measure non-parametric skewness
npSkew = function(x) { return ((mean(x)-median(x))/sd(x))}

npSkew_NA = unlist(lapply(lDist_NA, npSkew))
hist(npSkew_NA, breaks = 500)
npSkew_NA = merge(fData(preds2), data.frame(NPSKEW = npSkew_NA), by=0)[,c("svm.pred", "NPSKEW")]

ggplot(npSkew_NA) + geom_histogram(aes(x=NPSKEW, col=svm.pred))

group.colors = c("Basal Body 1" = "red", "Basal Body 2" = "deepskyblue4", 
                 "Pellicle" = "coral", "Trichocyst Matrix" = "blueviolet",
                 "Lysosome" = "orange", "Peroxisome" = "salmon", "Nuclei" = "blue4", 
                 "Membrane Trafficking 1" = "darkgoldenrod1", "Membrane Trafficking 2" = "cyan2", 
                 "Membrane Trafficking 3" = "brown",  
                 "Mitochondria-1" = "deeppink", "Mitochondria-2" = "darkorchid", 
                 "Mitochondria-3" = "limegreen", "Mitochondria-4" = "cornflowerblue",
                 "Cytosol" = "darkgreen", "Proteasome" = "aquamarine4", "Ribosome" = "darkgoldenrod",
                 "unknown" = "white")
ggplot(npSkew_NA) + geom_histogram(aes(x=NPSKEW, fill=svm.pred), bins=50) + 
  scale_fill_manual(values = group.colors) + xlab("Non-Parametric Skewness of Protein Similarity Score")

#save.image("./ptetPCP_rdata/ptetPCP_PSS.RData")
load("/ptetPCP/ptetPCP_rdata/ptetPCP_PSS.RData")

#####
#####
#####
##### NBAVG Data
range(unlist(lDist_NBAVG))

minEuc_nb = range(unlist(lDist_NBAVG))[1]
maxEuc_nb = range(unlist(lDist_NBAVG))[2]

  # Convert to Similarity Score
    # Normalize Distance, subtract by 1 to make Protein Similarity Score
normalize <- function(x) {return (1 - ((x - minEuc_nb) / (maxEuc_nb - minEuc_nb))) }
for(q in 1:length(lDist_NBAVG)){
  lDist_NBAVG[[names(lDist_NBAVG)[q]]] = sapply(lDist_NBAVG[[q]], FUN = normalize)
}
range(unlist(lDist_NBAVG))

  # nonparametric skew
npSkew = function(x) { return ((mean(x)-median(x))/sd(x))}

npSkew_NBAVG = unlist(lapply(lDist_NBAVG, npSkew))
hist(npSkew_NBAVG, breaks = 500)
npSkew_NBAVG = merge(fData(preds2), data.frame(NPSKEW = npSkew_NBAVG), by=0)[,c("svm.pred", "NPSKEW")]

ggplot(npSkew_NBAVG) + geom_histogram(aes(x=NPSKEW, col=svm.pred))

group.colors = c("Basal Body 1" = "red", "Basal Body 2" = "deepskyblue4", 
                 "Pellicle" = "coral", "Trichocyst Matrix" = "blueviolet",
                 "Lysosome" = "orange", "Peroxisome" = "salmon", "Nuclei" = "blue4", 
                 "Membrane Trafficking 1" = "darkgoldenrod1", "Membrane Trafficking 2" = "cyan2", 
                 "Membrane Trafficking 3" = "brown",  
                 "Mitochondria-1" = "deeppink", "Mitochondria-2" = "darkorchid", 
                 "Mitochondria-3" = "limegreen", "Mitochondria-4" = "cornflowerblue",
                 "Cytosol" = "darkgreen", "Proteasome" = "aquamarine4", "Ribosome" = "darkgoldenrod",
                 "unknown" = "white")
ggplot(npSkew_NBAVG) + geom_histogram(aes(x=NPSKEW, fill=svm.pred), bins=50) + 
  scale_fill_manual(values = group.colors) + xlab("Non-Parametric Skewness of Protein Similarity Score")

save.image("./ptetPCP_rdata/ptetPCP_PSS.RData")

#####
#####
#####
##### Yeast Data
load("/ptetPCP/ptetPCP_rdata/yeast2018_eucDist.RData")
  # Range of Euclidean Distances
range(unlist(lDist_yeast))
  # 0.02220629 1.08988061
minEucY = range(unlist(lDist_yeast))[1]
maxEucY = range(unlist(lDist_yeast))[2]

  # Convert to Similarity Score
   # Normalize Distance, subtract by 1 to make Protein Similarity Score
normalizeY <- function(x) {return (1 - ((x - minEucY) / (maxEucY - minEucY))) }
for(q in 1:length(lDist_yeast)){
  lDist_yeast[[names(lDist_yeast)[q]]] = sapply(lDist_yeast[[q]], FUN = normalizeY)
}
range(unlist(lDist_yeast))
  # 0 1

  # Make Plots
library(pRolocdata)
data("yeast2018")

setwd("/ptetPCP/ptetPCP_eucDist/yeast//")
for(compartment in unique(fData(yeast2018)$'predicted.location')){
  dir.create(compartment)
  setwd(compartment)
  genesC = fData(yeast2018)[which(fData(yeast2018)$'predicted.location' == compartment),'Accession']
  
  for(gene in genesC){
    if(length(lDist_yeast[[gene]]) > 0){
      png(paste(gene, "_hist.png", sep = ""))  # saves histogram with this name and next plot
      hist(lDist_yeast[[gene]], breaks = 500, main = as.character(gene), xlab = ("Protein Similarity Score"))
      dev.off()  # shut off plot saving
    }
  }
  setwd("../")  # go back to parent dir
}

  # Create an "average" distribution for each compartment
    # Sort so protein distributions are in order
lDist_yeast = lapply(lDist_yeast, sort)

setwd("/ptetPCP/ptetPCP_eucDist/yeast/")
for(compartment in unique(fData(yeast2018)$'predicted.location')){
  genesC = fData(yeast2018)[which(fData(yeast2018)$'predicted.location' == compartment),'Accession']
  png(paste(compartment, "_hist.png", sep = ""))
  hist(rowMeans(data.frame(lDist_yeast[which(names(lDist_yeast) %in% genesC)])), breaks = 500, 
       main = as.character(compartment), xlab = "Protein Similarity Score")  
  # normalize each part of list containing each organelle's proteins. rowmeans() does all the work after binding together a DF
  dev.off()
}

  # Measure nonparametric skew
npSkew = function(x) { return ((mean(x)-median(x))/sd(x))}

npSkew_yeast = unlist(lapply(lDist_yeast, npSkew))
hist(npSkew_yeast, breaks = 500)
npSkew_yeast = merge(fData(yeast2018), data.frame(NPSKEW = npSkew_yeast), by=0)[,c("predicted.location", "NPSKEW")]

group.colors_y = c("Mitochondrion" = "deeppink", "Proteasome regulatory subunit" = "aquamarine3", 
                   "Ribosome" = "darkgoldenrod","unknown" = "white", "Nucleus" = "blue4",
                   "Cytosol" = "darkgreen", "Histone" = "yellow", "Plasma membrane" = "cyan2",
                   "ER" = "darkgoldenrod1", "Golgi" = "darkgoldenrod3", "Vacuole" = "brown", 
                   "V ATPase" = "gray", "Proteasome catalytic subunit" = "aquamarine4")
ggplot(npSkew_yeast) + geom_histogram(aes(x=NPSKEW, fill=predicted.location), bins=50) + 
  scale_fill_manual(values = group.colors_y) + xlab("Non-Parametric Skewness of Protein Similarity Score")

rm(list=setdiff(ls(), c("lDist_NA", "lDist_NBAVG", "lDist_yeast", "npSkew_NA", "npSkew_NBAVG", "npSkew_yeast", "preds2", "yeast2018", "euclidean", "normalize")))

#save.image("./ptetPCP_rdata/ProteinSimilarityScores.RData")
load("./ptetPCP_rdata/ProteinSimilarityScores.RData")

#####
#####
#####
##### Toxo
load("./ptetPCP_rdata/toxo2020_eucDist.RData")

range(unlist(lDist_toxo))

minEuc_t = range(unlist(lDist_toxo))[1]
maxEuc_t = range(unlist(lDist_toxo))[2]

  # Convert to Similarity Score
    # Normalize Distance, subtract by 1 to make Protein Similarity Score
normalize <- function(x) {return (1 - ((x - minEuc_t) / (maxEuc_t - minEuc_t))) }
for(q in 1:length(lDist_toxo)){
  lDist_toxo[[names(lDist_toxo)[q]]] = sapply(lDist_toxo[[q]], FUN = normalize)
}
range(unlist(lDist_toxo))

  # Plots
setwd("/ptetPCP/ptetPCP_eucDist/toxo/")
for(compartment in unique(fData(Barylyuk2020ToxoLopit)$'tagm.mcmc.allocation.pred')){
  dir.create(as.character(compartment))
  setwd(as.character(compartment))
  genesC = fData(Barylyuk2020ToxoLopit)[which(fData(Barylyuk2020ToxoLopit)$'tagm.mcmc.allocation.pred' == as.character(compartment)),'Accession']
  
  for(gene in genesC){
    if(length(lDist_toxo[[gene]]) > 0){
      png(paste(gene, "_hist.png", sep = ""))  # saves histogram with this name and next plot
      hist(lDist_toxo[[gene]], breaks = 500, main = as.character(gene), xlab = ("Protein Similarity Score"))
      dev.off()  # shut off plot saving
    }
  }
  setwd("../")  # go back to parent dir
}

  # Create an "average" distribution for each compartment
    # Sort so protein distributions are in order
lDist_toxo = lapply(lDist_toxo, sort)

setwd("/ptetPCP/ptetPCP_eucDist/toxo/")
for(compartment in unique(fData(Barylyuk2020ToxoLopit)$'tagm.mcmc.allocation.pred')){
  genesC = fData(Barylyuk2020ToxoLopit)[which(fData(Barylyuk2020ToxoLopit)$'tagm.mcmc.allocation.pred' == as.character(compartment)),'Accession']
  png(paste(compartment, "_hist.png", sep = ""))
  hist(rowMeans(data.frame(lDist_toxo[which(names(lDist_toxo) %in% genesC)])), breaks = 500, 
       main = as.character(compartment), xlab = "Protein Similarity Score")  
  # normalize each part of list containing each organelle's proteins. rowmeans() does all the work after binding together a DF
  dev.off()
}

  # Measure nonparametric skew
npSkew = function(x) { return ((mean(x)-median(x))/sd(x))}

npSkew_toxo = unlist(lapply(lDist_toxo, npSkew))
hist(npSkew_toxo, breaks = 500)
npSkew_toxo = merge(fData(Barylyuk2020ToxoLopit), data.frame(NPSKEW = npSkew_toxo), by=0)[,c("tagm.mcmc.allocation.pred", "NPSKEW")]

ggplot(npSkew_toxo) + geom_histogram(aes(x=NPSKEW, fill=tagm.mcmc.allocation.pred), bins=50) + 
  scale_fill_manual(values = group.colors_y) + xlab("Non-Parametric Skewness of Protein Similarity Score")

rm(list=setdiff(ls(), c("lDist_toxo", "Barylyuk2020ToxoLopit", "npSkew_toxo")))

#save.image("/ptetPCP/ptetPCP_rdata/ProteinSimilarityScores.RData")
load("/ptetPCP/ptetPCP_rdata/ProteinSimilarityScores.RData")


