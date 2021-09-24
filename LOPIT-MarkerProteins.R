# This is a comprehensive analysis of marker protein distributions
library("MSnbase")
library("pRoloc")
library("pRolocGUI")

load("/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/RData-LOPIT/ptetLOPIT_12-36_impute_normalize.RData")

##### 12 Fractions
# Add Markers
ptetLOPIT12_NA_sumNormalized_markedB = addMarkers(ptetLOPIT12_NA_sumNormalized, markers = "/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/markerProts/markersBroad.csv", 
                                                  mcol = "markersBroad", fcol = "Accession", verbose = T)
ptetLOPIT12_NA_sumNormalized_markedS = addMarkers(ptetLOPIT12_NA_sumNormalized, markers = "/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/markerProts/markersSpecific.csv", 
                                                  mcol = "markersSpecific", fcol = "Accession", verbose = T)
ptetLOPIT12_NA_sumNormalized_markedY = addMarkers(ptetLOPIT12_NA_sumNormalized, markers = "/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/markerProts/markersYeast.csv", 
                                                  mcol = "markersYeast", fcol = "Accession", verbose = T)
ptetLOPIT12_NA_sumNormalized_markedTMP = addMarkers(ptetLOPIT12_NA_sumNormalized, markers = "/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/markerProts/markersTMP.csv", 
                                                  mcol = "markersTMP", fcol = "Accession", verbose = T)

# PCA
plot2D(ptetLOPIT12_NA_sumNormalized_markedB, fcol = "markersBroad")
addLegend(ptetLOPIT12_NA_sumNormalized_markedB, fcol = "markersBroad", where = "bottomright", ncol=10, cex = 0.75)
plot2D(ptetLOPIT12_NA_sumNormalized_markedS, fcol = "markersSpecific")
addLegend(ptetLOPIT12_NA_sumNormalized_markedS, fcol = "markersSpecific", where = "bottomright", ncol=10, cex = 0.75)
plot2D(ptetLOPIT12_NA_sumNormalized_markedY, fcol = "markersYeast")
addLegend(ptetLOPIT12_NA_sumNormalized_markedY, fcol = "markersYeast", where = "bottomright", ncol=2, cex = 0.75)

plot2D(ptetLOPIT12_NA_sumNormalized_markedTMP, fcol = "markersTMP")
kleb12 = featureNames(ptetLOPIT12_NA_sumNormalized_markedB[grep(pattern = "YP_*", featureNames(ptetLOPIT12_NA_sumNormalized_markedB))])
plot2D(object = ptetLOPIT12_NA_sumNormalized_markedB, fcol = NULL)     # Klebsiella Proteins
highlightOnPlot(object = ptetLOPIT12_NA_sumNormalized_markedB, foi = kleb12, col = "black", lwd = 2)  

# Distributions
plotDist(ptetLOPIT12_NA_sumNormalized_markedB)
#Markers in data: 69 out of 2388
#organelleMarkers
#BB          BB    Centrosome        Cilia           CV           ER          MAC          MIC Mitochondria   Peroxisome 
#5            1            1            6            1            1           14            1           20            3 
#PM     Ribosome   Trichocyst      unknown      Vesicle 
#4            7            2         2319            3 
bb = which(fData(ptetLOPIT12_NA_sumNormalized_markedB)["markersBroad"] == "BB")  # High in MAC, peaks in 3K. 
plotDist(ptetLOPIT12_NA_sumNormalized_markedB, markers = bb)
mit = which(fData(ptetLOPIT12_NA_sumNormalized_markedB)["markersBroad"] == "Mitochondria")
plotDist(ptetLOPIT12_NA_sumNormalized_markedB, markers = mit)                    # Peaks in 1K. Some peak in Sup
rib = which(fData(ptetLOPIT12_NA_sumNormalized_markedB)["markersBroad"] == "Ribosome")
plotDist(ptetLOPIT12_NA_sumNormalized_markedB, markers = rib)                   # Some peak in 1K. Some peak in 79/120K
ves = which(fData(ptetLOPIT12_NA_sumNormalized_markedB)["markersBroad"] == "Vesicle")
plotDist(ptetLOPIT12_NA_sumNormalized_markedB, markers = ves)                   # Moderate in 3-30K... 
er = which(fData(ptetLOPIT12_NA_sumNormalized_markedB)["markersBroad"] == "ER")
plotDist(ptetLOPIT12_NA_sumNormalized_markedB, markers = er)                    # only 1 ER protein... moderate in 3-30K
cv = which(fData(ptetLOPIT12_NA_sumNormalized_markedB)["markersBroad"] == "CV")
plotDist(ptetLOPIT12_NA_sumNormalized_markedB, markers = cv)                    # a few peaked in 30K.. generally high from 3-30K
perox= which(fData(ptetLOPIT12_NA_sumNormalized_markedB)["markersBroad"] == "Peroxisome")
plotDist(ptetLOPIT12_NA_sumNormalized_markedB, markers = perox)                 # gradual increased from 300 --> 1K --> 3K
golgi= which(fData(ptetLOPIT12_NA_sumNormalized_markedB)["markersBroad"] == "Golgi")
plotDist(ptetLOPIT12_NA_sumNormalized_markedB, markers = golgi)                 # peaking from 9-30K
mic= which(fData(ptetLOPIT12_NA_sumNormalized_markedB)["markersBroad"] == "MIC") 
plotDist(ptetLOPIT12_NA_sumNormalized_markedB, markers = mic)                   # only 1 MIC.... peaks in Sup
mac= which(fData(ptetLOPIT12_NA_sumNormalized_markedB)["markersBroad"] == "MAC")
plotDist(ptetLOPIT12_NA_sumNormalized_markedB, markers = mac)                   # peaks in MAC but is high throughout
act= which(fData(ptetLOPIT12_NA_sumNormalized_markedB)["markersBroad"] == "actin")
plotDist(ptetLOPIT12_NA_sumNormalized_markedB, markers = act)                   # none
trich= which(fData(ptetLOPIT12_NA_sumNormalized_markedB)["markersBroad"] == "Trichocyst")
plotDist(ptetLOPIT12_NA_sumNormalized_markedB, markers = trich)                   # Big peak in 300g and nowhere else
plotDist(ptetLOPIT12_NA_sumNormalized_markedB[kleb12], pcol = "blue")             # KLEB- Hard to tell... 

plotDist(ptetLOPIT12_NA_sumNormalized_markedY)
#Cytoplasm            ER Mitochondrial Mitochondrion       Nucleus        Plasma    Proteasome      Ribosome       unknown       Vacuole 
#2             7             2             8             2             1             8             7          2348             3 
erY = which(fData(ptetLOPIT12_NA_sumNormalized_markedY)["markersYeast"] == "ER")
plotDist(ptetLOPIT12_NA_sumNormalized_markedY, markers = erY)                    # peaks around 9-30K
mitoY = c(which(fData(ptetLOPIT12_NA_sumNormalized_markedY)["markersYeast"] == "Mitochondrial"),
            which(fData(ptetLOPIT12_NA_sumNormalized_markedY)["markersYeast"] == "Mitochondrion"))
plotDist(ptetLOPIT12_NA_sumNormalized_markedY, markers = mitoY)                    # looks identical to manual Mitos
riboY = which(fData(ptetLOPIT12_NA_sumNormalized_markedY)["markersYeast"] == "Ribosome")
plotDist(ptetLOPIT12_NA_sumNormalized_markedY, markers = riboY)                    # looks identical to manual Ribos
protY = which(fData(ptetLOPIT12_NA_sumNormalized_markedY)["markersYeast"] == "Proteasome")
plotDist(ptetLOPIT12_NA_sumNormalized_markedY, markers = protY)                    # increasing in 120K and peaks at Sup
cytoY = which(fData(ptetLOPIT12_NA_sumNormalized_markedY)["markersYeast"] == "Cytoplasm")
plotDist(ptetLOPIT12_NA_sumNormalized_markedY, markers = cytoY)                    # only 1/3 peaks in Sup. another peaks in 30K?

#####
#####
##### Let's get into the GUI
pRolocGUI::pRolocVis(ptetLOPIT12_NA_sumNormalized_markedB, app = "pca", fcol = "markersBroad")
  # One by one, going through each group and adding paralogs in the dataset that weren't originally in the marker set

#### Basal Bodies
  # After adding all paralogs and associated genes, I can upload my new list
  #It appears as though there are 3 clusters for the BBs. I'll separate those below
    # 1st: Largest- Peaks in 300g
    # 2nd: Intermediate- Peaks in 3K
    # 3rd: Smallest- Peaks in Sup
  # Are these 3 clusters biologically meaningful??
ptetLOPIT12_NA_sumNormalized_markedBB = addMarkers(ptetLOPIT12_NA_sumNormalized, markers = "/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/markerProts/markersBB.csv", 
                                                  mcol = "markersBB", fcol = "Accession", verbose = T)  # Add our new markers. 64 total
tmpN = which(fData(ptetLOPIT12_NA_sumNormalized_markedBB)$"markersBB" == "BB")    # Separate marker proteins for BB
length(tmpN)                                                                      # 45 of 64 BB proteins in dataset
tmpDF = exprs(ptetLOPIT12_NA_sumNormalized_markedBB[tmpN])                        # rows for BB markers
tmpMax = colnames(tmpDF)[apply(tmpDF,1,which.max)]                                # column with largest value
table(tmpMax)                                                   
#   MAC   Sup  X15K X300g   X3K 
#   24     4     1     8     8
earlyBB = rownames(tmpDF[which(tmpMax == "MAC" | tmpMax == "X300g"),])            # 1st group: 300g and 1K
midBB = rownames(tmpDF[which(tmpMax == "X3K"),])                                  # 2nd group: 3K
lateBB = rownames(tmpDF[which(tmpMax == "X15K" | tmpMax == "Sup"),])              # Thirs group: Sup (and 15K, because where else would it go?)
# PCA with each group highlighted to show good overlap between peaked columns and clusters
plot2D(object = ptetLOPIT12_NA_sumNormalized_markedBB, fcol = NULL)             
highlightOnPlot(object = ptetLOPIT12_NA_sumNormalized_markedBB, foi = earlyBB, col = "black", lwd = 2) 
highlightOnPlot(object = ptetLOPIT12_NA_sumNormalized_markedBB, foi = midBB, col = "red", lwd = 2) 
highlightOnPlot(object = ptetLOPIT12_NA_sumNormalized_markedBB, foi = lateBB, col = "green", lwd = 2) 

#####
##### Back to the GUI
pRolocGUI::pRolocVis(ptetLOPIT12_NA_sumNormalized_markedB, app = "pca", fcol = "markersBroad")

###### Actin
  # Very few actine genes in data, so this is a place where using GO terms will be helpful
  # Loading in all genes for the GOterm: 0015629 (Cytoskeleton)
ptetLOPIT12_NA_sumNormalized_markedCytoskel = addMarkers(ptetLOPIT12_NA_sumNormalized, markers = "/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/markerProts/markersCytoskeleton.csv", 
                                                   mcol = "markersCytosk", fcol = "Accession", verbose = T)  # Add our new markers. 64 total
plot2D(object = ptetLOPIT12_NA_sumNormalized_markedCytoskel, fcol = "markersCytosk")             
pRolocGUI::pRolocVis(ptetLOPIT12_NA_sumNormalized_markedCytoskel, app = "pca", fcol = "markersCytosk")
  # Really interesting observation that PTET.51.1.P0710058 is the only protein spiking in the 3K fraction, while the rest peak in the Sup
  # There's another protein peaking in 3K with a far lower relative intensity: PTET.51.1.P0600063
    # THis is a paralog of a protein studied in Sehring et al. 2007 
  # Another interesting protein is PTET.51.1.P0900031, which is higher in 9K, 12K, and 15K, but still peaks in the Sup. Maybe associated with heavier things? e.g. food vacuoles
    # It's duplicate is PTET.51.1.P0770188, which only peaks in the SUp and doesnt increase elsewhere
      # PTET.51.1.P0770188 clusters with the rest of the cytoskeleton
  # Last one is PTET.51.1.P0320240, which peaks at 1K. Only protein of this group to do so
    # It's duplicate is PTET.51.1.P0200282, which only peaks in Sup and clusters with most other Cytoskeleton
tmpN = which(fData(ptetLOPIT12_NA_sumNormalized_markedCytoskel)$"markersCytosk" == "Cytoskeleton")    # Separate marker proteins for cytoskeleton
length(tmpN)                                                                      # 31 of 61 proteins in dataset
tmpDF = exprs(ptetLOPIT12_NA_sumNormalized_markedCytoskel[tmpN])                        # rows for cyto markers
tmpMax = colnames(tmpDF)[apply(tmpDF,1,which.max)]                                # column with largest value
table(tmpMax)                # let's break this into heavy and light clusters                                   
#   MAC Sup X1K X3K X5K 
#   3   24   1   2   1 
heavyCytosk = rownames(tmpDF[which(tmpMax == "MAC" | tmpMax == "X300g" | tmpMax == "X1K" | tmpMax == "X3K" | tmpMax == "X5K"),])
lightCytosk = rownames(tmpDF[which(tmpMax == "Sup"),])                      
plot2D(object = ptetLOPIT12_NA_sumNormalized_markedCytoskel, fcol = NULL)             
highlightOnPlot(object = ptetLOPIT12_NA_sumNormalized_markedCytoskel, foi = heavyCytosk, col = "black", lwd = 2) 
highlightOnPlot(object = ptetLOPIT12_NA_sumNormalized_markedCytoskel, foi = lightCytosk, col = "red", lwd = 2) 

#####
##### Back to the GUI
pRolocGUI::pRolocVis(ptetLOPIT12_NA_sumNormalized_markedB, app = "pca", fcol = "markersBroad")













ptetLOPIT36_NA_sumNormalized_markedB = addMarkers(ptetLOPIT36_NA_sumNormalized, markers = "/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/markerProts/markersBroad.csv", 
                     mcol = "markersBroad", fcol = "Accession", verbose = T)
ptetLOPIT36_NA_sumNormalized_markedS = addMarkers(ptetLOPIT36_NA_sumNormalized, markers = "/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/markerProts/markersSpecific.csv", 
                     mcol = "markersSpecific", fcol = "Accession", verbose = T)
ptetLOPIT36_NA_sumNormalized_markedY = addMarkers(ptetLOPIT36_NA_sumNormalized, markers = "/data/LynchLabCME/Paramecium/Tim/Paramecium_LOPIT/markerProts/markersYeast.csv", 
                     mcol = "markersYeast", fcol = "Accession", verbose = T)

