# ptetPCP_Markers.R

setwd("/ptetPCP/")
load("./ptetPCP_rdata/ptetPCP_rawMSN.RData")

library(pRoloc)

filterNA(ptetPCP_raw_all_norm_NA)
filterNA(ptetPCP_raw_all_nbavg_norm_NA)
filterNA(ptetPCP_raw_all_nbavg_norm_zero_NA)

marked = addMarkers(ptetPCP_raw_all_nbavg_norm_zero_NA, markers = "./ptetPCP_csv/markersFinal.csv", 
                    mcol = "markers", fcol = "Accession", verbose = T)

#save.image("./ptetPCP_rdata/ptetPCP_marked.RData")
load("./ptetPCP_rdata/ptetPCP_marked.RData")

# tSNE with Markers
set.seed(42)
plot2D(marked, fcol = "markers", method = "t-SNE")
addLegend(marked, fcol = "markers", where = "topright", bty = "n", cex = .6)

# Marker Protein Resolution
mrkHClust(marked, fcol = "markers")
qs = QSep(marked)
plot(qs)
levelPlot(qs)

# Distributions
  # Nuclei
plotDist(marked[which(fData(marked)$'markers' == "Basal Body 1")], pcol = "red", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Basal Body 2")], pcol = "deepskyblue4", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Cytosol")], pcol = "darkgreen", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Lysosome")], pcol = "orange", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Membrane Trafficking 1")], pcol = "darkgoldenrod1", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Membrane Trafficking 2")], pcol = "cyan2", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Membrane Trafficking 3")], pcol = "brown", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Mitochondria-1")], pcol = "deeppink", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Mitochondria-2")], pcol = "darkorchid", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Mitochondria-3")], pcol = "limegreen", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Mitochondria-4")], pcol = "cornflowerblue", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Nuclei")], pcol = "blue4", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Pellicle")], pcol = "coral", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Peroxisome")], pcol = "salmon", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Proteasome")], pcol = "aquamarine4", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Ribosome")], pcol = "darkgoldenrod", ylim = c(0,1))
plotDist(marked[which(fData(marked)$'markers' == "Trichocyst Matrix")], pcol = "blueviolet", ylim = c(0,1))



#####
##### 
#####
  # Nuclear: Nuclei 
    # Dist: MAC/300g + 9K-30K
exprs(ptetPCP_raw_all_nbavg_norm_NA["PTET.51.1.P1370127"])  # RPB1
plotDist(ptetPCP_raw_all_nbavg_norm_NA["PTET.51.1.P1370127"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0480005"])  # RPB2 (missing 2 columns)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0480005"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1410114"])  # PtRPC
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1410114"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0610112"])  # TFIIS3 (missing many columns)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0610112"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0540024"])  # LIG4a (missing a few columns)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0540024"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0090213"])  # PTETG900022001-Poly(A) Processing
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0090213"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0890115"])  # LIGK03 (high in MAC and 30K, missing many)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0890115"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1740061"])  # UNKNOWN- Topoisomerase I
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P1740061"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0380097"])  # UNKNOWN- DNA Helicase
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0380097"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0040268"])  # UNKNOWN- Pre-mRNA splicing factor
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0040268"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1360088"])  # UNKNOWN- DNA-directed RNA polymerase subunit 2
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P1360088"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1370064"])  # H3P3 (missing a few)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1370064"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0550265"])  # H3P2
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0550265"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010480"])  # PTMB.80c(H1) (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010480"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1470101"])  # act4-1 (missing a few)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1470101"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0190343"])  # alp1-1 (missing a few)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0190343"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0010578"])  # PTMB.03(Nup)
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0010578"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0750016"])  # UNKNOWN- Peptidase S59, nucleoporin
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0750016"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0450164"])  # UNKNOWN- Utp11
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0450164"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0010267"])  # PTMB.259(nucleoli)
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0010267"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0010095"])  # PTMB.400c(Helicase)
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0010095"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0280244"])  # UNKOWN- Nucleolar, Nop52
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0280244"], fcol = NULL, pcol = "red", ylim = c(0,1))

  # Cortical: Basal Body Core
    # Dist: MAC/300g only
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0240208"])  # BLD10
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0240208"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0990148"])  # EPI49
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0990148"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1180141"])  # PtCen2a (missing 1 column)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1180141"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0060088"])  # EPI131 (missing a few columns)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0060088"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0520250"])  # Epi34
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0520250"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1400098"])  # Epi51-5a (missing many columns)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1400098"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0250081"])  # Epi38-2b (missing many columns)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0250081"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0030411"])  # Epi30-3a (missing many columns)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0030411"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0760096"])  # Epi2-2a (missing many columns)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0760096"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0900049"])  # Epi41-1a  (missing many columns)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0900049"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1290151"])  # Epi11-4a  (missing many columns)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1290151"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0270102"])  # Epi46-3b  (missing many columns)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0270102"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0340252"])  # RPGRIP1L/MKS5/NPHP8
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0340252"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0760096"])  # Epi2-2a (missing many columns)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0760096"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0520183"])  # Epi43 (missing many columns)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0520183"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0170217"])  # FOR20a (only abundant in MAC and 300g)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0170217"], fcol = NULL, pcol = "red", ylim = c(0,1))

  # Cortical: Basal Body Rootlet
    # MAC/300g + 3K
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1110112"])  # PCM1 (1 missing value) cortical
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1110112"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1640062"])  # PCM4 (1 missing value)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1640062"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0170094"])  # Bug22p cilia
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0170094"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0900028"])  # SFA1a (1 missing value)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0900028"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0670013"])  # SFA2 (many missing values)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0670013"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0220060"])  # SFA4 (many missing values)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0220060"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0940007"])  # SFA5a
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0940007"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0290029"])  # SFA6a
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0290029"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1300010"])  # SFA8a
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P1300010"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1480020"])  # SFA10a
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P1480020"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0040427"])  # SFA11a (missing many columns)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0040427"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0110226"])  # SFA12d
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0110226"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0890188"])  # SFA13a
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0890188"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1030030"])  # KdD5 (a few missing columns)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1030030"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0270089"])  # Ptcen_icl5a
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0270089"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0700039"])  # Ptcen_icl7a
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0700039"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1280088"])  # Ptcen_icl8a
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1280088"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0180150"])  # Ptcen27
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0180150"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1440071"])  # PtCenBP1
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1440071"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0890151"])  # tub_alphaPT2
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0890151"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1500039"])  # tub_betaPT3
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1500039"], fcol = NULL, pcol = "red", ylim = c(0,1))

  # Cortical: Pellicle
    # Dist: No MAC + 3K-5K + 12K-15K
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0010348"])  # mag3 (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010348"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0410067"])  # mag5 (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0410067"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0870001"])  # mag6 (few NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0870001"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1410053"])  # mag7
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P1410053"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1060180"])  # sAG_51A
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P1060180"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1590123"])  # sAG_alpha51D
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P1590123"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1280004"])  # sAG_epsilon51D (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1280004"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1420113"])  # UNKNOWN- Paramecium surface antigen
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P1420113"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0120257"])  #  UNKNOWN- Paramecium surface antigen (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0120257"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1530070"])  #  UNKNOWN- Paramecium surface antigen
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P1530070"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0230216"])  # UNKNOWN- G surface protein, allelic form 168 (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0230216"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1430110"])  # PTETG14300005001(surfaceAntigen)
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P1430110"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0430141"])  # PTETG4300022001(surfaceAntigen) (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0430141"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0440143"])  # PTETG4400015001(surfaceAntigen) (few NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0440143"], fcol = NULL, pcol = "red", ylim = c(0,1))

  # Cortical: Trichocyst Matrix
    # Dist: 300g only
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0200318"])  # TMP2a (few NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0200318"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0260014"])  # TMP1h
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0260014"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1540054"])  # TMP4e (few NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1540054"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0020298"])  # TMP52c (few NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0020298"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010552"])  # TMP33a (one NA)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010552"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1470130"])  # TMP21e (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1470130"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1490007"])  # TMP24b (few NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1490007"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0150197"])  # TMP51e (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0150197"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0180247"])  # TMP31d
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0180247"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0180301"])  # TMP23a (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0180301"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1840025"])  # TMP11g (few NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1840025"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1590103"])  # TMP25b (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1590103"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0160048"])  # TMP12f (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0160048"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0260111"])  # TMP35b (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0260111"], fcol = NULL, pcol = "red", ylim = c(0,1))

  # Membrane Trafficking: Vacuole-1 (Phagolysosomal)
    # Dist: 3K-5K sharp decline + 3:9K dip
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0670042"])  # TGL6
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0670042"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0760066"])  # PLA1 (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0760066"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1150002"])  # UNKNOWN- Peptidase C1A
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P1150002"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1620036"])  # UNKNOWN- Peptidase C1A (few NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1620036"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1840043"])  # UNKNOWN- Peptidase C1A 
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P1840043"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0210028"])  # UNKNOWN- Peptidase C1A (few NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0210028"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0220279"])  # UNKNOWN- Peptidase C1A 
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0220279"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0360163"])  # UNKNOWN- Peptidase C1A (few NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0360163"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0450177"])  # UNKNOWN- Peptidase C1A 
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0450177"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0460227"])  # UNKNOWN- Peptidase C1A (few NAs) 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0460227"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0520159"])  # UNKNOWN- Peptidase C1A (few NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0520159"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0560051"])  # UNKNOWN- Peptidase C1A 
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0560051"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0720197"])  # UNKNOWN- Peptidase C1A (few NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0720197"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0750188"])  # UNKNOWN- Peptidase C1A 
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0750188"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0090011"])  # UNKNOWN- Peptidase C1A 
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0090011"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0980127"])  # UNKNOWN- Peptidase C1A 
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0980127"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0730225"])  # UNKNOWN- Ribonuclease T2-like
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0730225"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0150174"])  # UNKNOWN- Lysosomal cystine transporter (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0150174"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0060101"])  # UNKNOWN- Gamma interferon inducible lysosomal thiol reductase GILT
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0060101"], fcol = NULL, pcol = "red", ylim = c(0,1))

  # Membrane Trafficking: Vacuole-2  (ER Chaperones?)
    # Dist: 9K-30K + 1:12K dip
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0020398"])  # PtSyx1-1 (many NAs)  exocytosis
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0020398"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0640022"])  # ptSERCA1 ER+Alveola
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0640022"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1020045"])  # Hsp70Pt08
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P1020045"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0980088"])  # pdi1_1 (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0980088"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0100050"])  # UNKNOWN- Endoplasmic reticulum oxidoreductin 1
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0100050"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0560116"])  # UNKOWN- DnaJ domain (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0560116"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0850007"])  # UNKNOWN- Heat shock protein 70 family
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0850007"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0250198"])  # UNKNOWN- Heat shock protein 70 family
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0250198"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0250311"])  # UNKNOWN- HSP40/DnaJ peptide-binding
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0250311"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0340057"])  # UNKNOWN- Stress-associated endoplasmic reticulum protein
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0340057"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0460141"])  # UNKNOWN- Retrieval of early ER protein Rer1
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0460141"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0910101"])  # UNKNOWN- Translocation protein Sec63
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0910101"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1220012"])  # UNKOWN- Glycosyl transferase, family 1 (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1220012"], fcol = NULL, pcol = "red", ylim = c(0,1))


  # Membrane Trafficking: Vacuole-3 (Trafficking to CV; Golgi)
    # Dist: 300g + 9K + 30K + 3:peaked 15/30K
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0760174"])  # PtSyb2-1 (many NAs) CV
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0760174"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P1670084"])  # STO1c CV
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P1670084"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0410185"])  # NSF2 CV
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0410185"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0340190"])  # IP3Rn-2 CV
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0340190"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0810071"])  # PtSyx3-2 (many NAs) endosome
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0810071"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0750126"])  # PtSec1
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0750126"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0070093"])  # PtSyx14-1 (many NAs) CV
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0070093"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0240212"])  # PtSec22 (few NAs)  golgi
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0240212"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0740007"])  # rab_C61
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0740007"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1000099"])  # rab_B23
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1000099"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0090039"])  # rab_B27
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0090039"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0150054"])  # rab_B11
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0150054"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0580179"])  # PtSto3a
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0580179"], fcol = NULL, pcol = "red", ylim = c(0,1))

# Membrane Trafficking: Vacuole-4 (Trafficking to PM)
  # Dist: 300g + 15K + Sup
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010071"])  # CK2alpha1-1 cortical/cytoplasmic
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010071"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0040411"])  # DHC-8 cytoplasmic dynein
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0040411"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0850133"])  # Act1-5 diffuse cyto
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0850133"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0400134"])  # Act2_1 diffuse cyto
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0400134"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0140342"])  # arp2-1
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0140342"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0620211"])  # arp3-1
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0620211"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0070002"])  # alphaPT8
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0070002"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0840194"])  # IFT46 Cortex/cilia
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0840194"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1400023"])  # CaNA4a Cortex
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1400023"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1380073"])  # syb8-1 (many NAs) cortical vesicles
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1380073"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0150355"])  # DRPD endosomes
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0150355"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1130109"])  # PtPP2C ER/Cilia
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1130109"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0700129"])  # PKAc2-2
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0700129"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0300272"])  # UNKNOWN- Armadillo-type fold
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0300272"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1820010"])  # UNKNOWN- Clathrin, heavy chain/VPS, 7-fold repeat
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1820010"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0140233"])  # UNKNOWN- Intraflagellar transport protein 81
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0140233"], fcol = NULL, pcol = "red", ylim = c(0,1))

  # Peroxisomal: Peroxisome
    # Dist: 3K peak with lots of variability
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0810085"])  # THIKa
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0810085"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0770146"])  # UNKNOWN- Thiolase-like, subgroup
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0770146"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_norm_NA["PTET.51.1.P0600217"])  # UNKNOWN- Isocitrate Lyase
plotDist(ptetPCP_raw_all_norm_NA["PTET.51.1.P0600217"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0590096"])  # UNKNOWN- Peroxisome membrane protein, Pex16 (many NAs)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0590096"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0160203"])  # UNKNOWN- Peroxisome membrane protein, Pex16
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0160203"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0050397"])  # UNKNOWN- Peroxisomal biogenesis factor 11
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0090034"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0140220"])  # UNKNOWN- Acyl-CoA dehydrogenase/oxidase C-terminal
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0140220"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0200182"])  # UNKNOWN- Acyl-CoA dehydrogenase/oxidase, N-terminal and middle domain
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0200182"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0700161"])  # UNKNOWN- Acyl-CoA dehydrogenase/oxidase C-terminal
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0700161"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0350075"])  # UNKNOWN- Peroxisomal fatty acyl CoA transporter, transmembrane domain
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0350075"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0880188"])  # UNKNOWN- Peroxisomal biogenesis factor 11
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0880188"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1180146"])  # UNKNOWN- Acyl-CoA dehydrogenase/oxidase, N-terminal and middle domain
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1180146"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0650053"])  # UNKNOWN- Isocitrate dehydrogenase NADP-dependent
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0650053"], fcol = NULL, pcol = "red", ylim = c(0,1))

  # Mitochondrial: Mitochondria-1
    # 300g-1K (sometimes 3-5K)
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF2"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF2"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF71"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF71"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF102"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF102"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF122"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF122"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF123"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF123"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF222"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF222"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF223"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF223"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF224"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF224"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF275"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF275"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0340280"])  # RPL3a
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0340280"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0350267"])  # RPL28a
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0350267"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0020194"])  # UNKNOWN- GCA2a
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0020194"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0120293"])  # NAD11a
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0120293"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1440138"])  # NDUFAB2a
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1440138"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1560071"])  # SDH1a
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1560071"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0040135"])  # QCR1a
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0040135"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0330242"])  # QCR2a
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0330242"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010254"])  # SURF1a
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010254"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0970034"])  # ATPTT2a
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0970034"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0170096"])  # UNKNOWN- Tim23
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0170096"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010208"])  # PTMB.302c(TIM10)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010208"], fcol = NULL, pcol = "red", ylim = c(0,1))


# Mitochondrial: Mitochondria-2
  # 300g-1K + 1:5K + 3:12K
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF6"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF6"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF17"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF17"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF18"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF18"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF28"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF28"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF32"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF32"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF38"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF38"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF72"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF72"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF82"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF82"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF86"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF86"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF99"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF99"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF100"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF100"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF101"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF101"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF104"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF104"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF138"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF138"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF140"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF140"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF144"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF144"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0160089"])  # CBP3
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0160089"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0020022"])  # SCO1
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0020022"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0020234"])  # UNKNOWN- succinate dehydrogenase
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0020234"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1570128"])  # UNKNOWN- Mitochondrial carrier domain
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1570128"], fcol = NULL, pcol = "red", ylim = c(0,1))

# Mitochondrial: Mitochondria-3
  # 300g-1K + Sup
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF60"])  # 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["ORF60"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1250149"])  # UNKNOWN- Phosphoenolpyruvate carboxykinase, ATP-utilising
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1250149"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0120175"])  # UNKNOWN- isocitrate dehydrogenase
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0120175"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010339"])  # PTMB.197c(Aldo/Keto Reductase)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010339"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0590214"])  # Enolase
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0590214"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0360202"])  # UNKNOWN- Phosphoenolpyruvate carboxykinase, ATP-utilising
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0360202"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1240009"])  # UNKNOWN- Succinyl-CoA synthetase, beta subunit
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1240009"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0430154"])  # UNKNOWN- Citrate synthase-like
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0430154"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0690157"])  # UNKNOWN- NAD(P)-binding domain
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0690157"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0040212"])  # UNKNOWN- Succinyl-CoA synthetase, beta subunit
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0040212"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0030311"])  # UNKNOWN- Succinyl-CoA ligase, alpha subunit
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0030311"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0040159"])  # UNKNOWN- Mitochondrial carrier domain
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0040159"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0720074"])  # UNKNOWN- NAD(P)-binding domain
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0720074"], fcol = NULL, pcol = "red", ylim = c(0,1))

# Mitochondrial: Mitochondria-4
  # 300g-1K + 1:30K + 2/3:30K+120K
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0280026"])  # Tom40a
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0280026"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0480261"])  # Tom40b
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0480261"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0520173"])  # KdC1(Porin)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0520173"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0600104"])  # KdC2(Porin) 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0600104"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0950122"])  # KdC3(Porin)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0950122"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P3030001"])  # UNKNOWN- Molybdopterin synthase catalytic subunit 2 (Porin)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P3030001"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0520211"])  # PtSep4
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0520211"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0260012"])  # UNKNOWN- Glutathione S-transferase, C-terminal-like
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0260012"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0370192"])  # UNKNOWN- Amidase
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0370192"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0910128"])  # UNKNOWN- Aldehyde/histidinol dehydrogenase
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0910128"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0170080"])  # UNKNOWN- NAD(P)-binding domain
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0170080"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0480063"])  # UNKNOWN- Phosphofructokinase domain
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0480063"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0920146"])  # UNKNOWN- Aldehyde/histidinol dehydrogenase
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0920146"], fcol = NULL, pcol = "red", ylim = c(0,1))

  # Cytoplasmic: Proteasome
    # Dist: 120K rise to Sup + 30K + sometimes 300g small bump
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010416"])  # PTMB.133c(20S Proteasome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010416"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0800115"])  # UNKNOWN- Nucleophile aminohydrolases, N-terminal(20S Proteasome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0800115"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0030132"])  # UNKNOWN- Nucleophile aminohydrolases, N-terminal(20S Proteasome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0030132"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0870157"])  # UNKNOWN- Nucleophile aminohydrolases, N-terminal(20S Proteasome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0870157"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0100257"])  # UNKNOWN- Nucleophile aminohydrolases, N-terminal(20S Proteasome) 
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0100257"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1030078"])  # UNKNOWN- Nucleophile aminohydrolases, N-terminal(20S Proteasome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1030078"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1120142"])  # UNKNOWN- Nucleophile aminohydrolases, N-terminal(20S Proteasome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1120142"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0120126"])  # UNKNOWN- Nucleophile aminohydrolases, N-terminal(20S Proteasome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0120126"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0080262"])  # UNKNOWN- Nucleophile aminohydrolases, N-terminal(20S Proteasome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0080262"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0200281"])  # UNKNOWN- Nucleophile aminohydrolases, N-terminal(20S Proteasome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0200281"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0560080"])  # UNKNOWN- P-loop containing nucleoside triphosphate hydrolase(26S Proteasome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0560080"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0680085"])  # UNKNOWN- Proteasome component (PCI) domain(26S Proteasome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0680085"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0820051"])  # UNKNOWN- PCI/PINT associated module(26S Proteasome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0820051"], fcol = NULL, pcol = "red", ylim = c(0,1))

# Cytoplasmic: Ribosome
  # Dist: 79K-120K bump
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0590250"])  # PTETG5900030001(40S Ribosome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0590250"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0050386"])  # UNKNOWN- Plectin/S10, N-terminal(40S Ribosome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0050386"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0180136"])  # UNKNOWN- Ribosomal protein S7e(40S Ribosome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0180136"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0630160"])  # UNKNOWN- K homology domain, prokaryotic type(40S Ribosome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0630160"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0520192"])  # UNKNOWN- Ribosomal protein S5 domain 2-type fold(40S Ribosome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0520192"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0800141"])  # Ribosomal protein L4 domain(60S Ribosome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0800141"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0510170"])  # UNKNOWN- Nucleotide-binding, alpha-beta plait(60S Ribosome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0510170"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0110039"])  # UNKNOWN- Ribosomal protein L10e/L16(60S Ribosome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0110039"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0190090"])  # UNKNOWN- Ribosomal protein L13(60S Ribosome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0190090"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0210110"])  # UNKNOWN- Ribosomal protein L35A(60S RIbosome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0210110"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0050332"])  # PTETG500001001(60S Ribosome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0050332"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1370018"])  # PTETG13700004001(60S Ribosome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1370018"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0320288"])  # UNKNOWN- Ribosomal protein L5 eukaryotic/L18 archaeal(60S Ribosome)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0320288"], fcol = NULL, pcol = "red", ylim = c(0,1))


# Cytoplasmic: Cytosol
  # Dist: Sup Peak
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0330286"])  # UNKNOWN- Alanine--tRNA ligase, chloroplastic/mitochondrial
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0330286"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0080058"])  # UNKNOWN- Phenylalanyl-tRNA synthetase
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0080058"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0150089"])  # UNKNOWN- Serine-tRNA ligase, type1
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0150089"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0360019"])  # UNKNOWN- aminoacyl-tRNA ligase
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0360019"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1440110"])  # UNKNWN- Valine-tRNA ligase
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1440110"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0770173"])  # UNKNOWN- Fructose-bisphosphate aldolase class-I, eukaryotic-type
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0770173"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0210289"])  # UNKNOWN- Rossmann-like alpha/beta/alpha sandwich fold
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0210289"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0020142"])  # UNKNOWN- Alcohol dehydrogenase superfamily, zinc-type
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0020142"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0040389"])  # UNKNOWN- Fructose-1,6-bisphosphatase class 1/Sedoheputulose-1,7-bisphosphatase
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0040389"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0210130"])  # UNKNOWN- Aminoacyl-tRNA synthetase, class II (D/K/N)-like
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0210130"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0410147"])  # UNKNOWN- Malate dehydrogenase, type 2
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0410147"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0470183"])  # UNKNOWN- Serine/threonine protein phosphatase 5
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0470183"], fcol = NULL, pcol = "red", ylim = c(0,1))
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1430120"])  # UNKNOWN- Alanyl-tRNA synthetase, class IIc, core domain
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1430120"], fcol = NULL, pcol = "red", ylim = c(0,1))




  # Lipid Droplet?
    # Dist: 300g + 9K + 15K + Sup
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0660104"])  # UNKNOWN- Glycolipid transfer protein domain
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0660104"], fcol = NULL, pcol = "red")
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1160095"])  # UNKNOWN- Glycolipid transfer protein domain
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1160095"], fcol = NULL, pcol = "red")
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0170086"])  # UNKNOWN- Cation-transporting  P-type ATPase, subfamily IV
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0170086"], fcol = NULL, pcol = "red")
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0580022"])  # UNKNOWN- Phospholipid-transporting P-type ATPase
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0580022"], fcol = NULL, pcol = "red")
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0160163"])  # UNKNOWN- Phospholipid-transporting ATPase IG
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0160163"], fcol = NULL, pcol = "red")
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010203"])  # PTMB.306(Lipid-droplet associated hydrolase)
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0010203"], fcol = NULL, pcol = "red")
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1050178"])  # UNKNOWN- Phospholipid/glycerol acyltransferase
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P1050178"], fcol = NULL, pcol = "red")
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0080048"])  # UNKNOWN- Phospholipid/glycerol acyltransferase
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0080048"], fcol = NULL, pcol = "red")
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0260075"])  # UNKNOWN- Tensin phosphatase, lipid phosphatase domain
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0260075"], fcol = NULL, pcol = "red")
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0790198"])  # UNKNOWN- Phospholipid:diacylglycerol acyltransferase 1
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0790198"], fcol = NULL, pcol = "red")
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0300072"])  # vhac2
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0300072"], fcol = NULL, pcol = "red")
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0840033"])  # UNKNOWN- Adenylyl cyclase class-3/4/guanylyl cyclase
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0840033"], fcol = NULL, pcol = "red")
exprs(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0660104"])  # UNKNOWN- Glycolipid transfer protein domain
plotDist(ptetPCP_raw_all_nbavg_norm_zero_NA["PTET.51.1.P0660104"], fcol = NULL, pcol = "red")




