# Work flow for Licknack et al. (Spatial Proteomics of Paramecium tetraurelia)

Raw Data Location: ____

Scripts in order, external steps, plots/figures/data:

1) ptetPCP_CSVMtoMSN.R
  - Goal: Converting from protein-level csv to MSnBase S4 object
  - Packages: pRoloc, tidyverse
  - Custom R Functions: NAprocessPCP_SUP.R; NAprocessPCP_MAC.R
  - RData: ptetPCP_rawMSN.RData
  - Plots: PCA, t-SNE (w/ + w/o imputation), venn diagram
  - Output Files: ptetPCP_norm.csv; ptetPCP_nbavg_norm.csv; ptetPCP_nbavg_norm_zero.csv
  
2) ptetPCP_Markers.R
  - Goal: Overview of marker proteins for all 17 organellar compartments
  - Packages: pRoloc
  - External: These were manually compiled and annotated in Supp Table 3 (Licknack et al. 2022)
  - RData: ptetPCP_marked.RData
  - Plots: t-SNE (imputed+markers); Dendrogram (markers), Protein Abundance Distributions (markers)
  
3) ptetPCP_DeNovoClustering.R
  - Goal: KNN Clustering (k=17; same number as marker classes)
  - Packages: pRoloc
  - RData: ptetPCP_knn17.RData
  - Plots: t-SNE (knn clusters)
  - Output Files: ptetPCP_knn17.csv
  
4) ptetPCP_SupervisedClassification.R
  - Goal: SVM Classification for 17 organellar compartments
  - Packages: pRoloc, tidyverse
  - RData: ptetPCP_SVM.Rdata
  - Plots: t-SNE (preds, preds-filtered), Dot Plots (SVM vs KNN; sorted)
  - Output Files: SVMpredictions.csv, SVMpredictions_Summary.csv
 
No Particular Order to the following:

 A) ptetPCP_Properties.R
  - Goal: Make plots for gene properties of predicted organellar proteins
  - Package: tidyverse, tm, wordcloud
  - Custom R Functions: makeWordCloud.R
  - External: proteins were plugged into Biomart/Sherlock (https://paramecium.i2bc.paris-saclay.fr/page/tools) and merged in Excel
            
            - additionally, protein sequences were plugged into TargetP and merged in the same way
              
              - File: ptetPCP_Properties.csv
  - Plots: Violin Plot (mRNA; PSM); Regression (mRNA vs PSM); Wordclouds; Dot Plots (Signal Peptides; Target Peptides; TM Domains); Box Plots (Amino Acids, pI; piNpiS)
    - piNpiS plots require: ptetPCP_pogenomics.R
       - This simply modifies a file obtained from Johri et al. 2018
       - RData: ptetPCP_popGenome.Rdata
   
   
  B1) ptetPCP_makeEucDistDB.R
    - Goal: Generate all pairwise Euclidean Distances 
    - External: Recomend using a computing cluster or something strong
    - RData: ptetPCP_eucDist.RData
  
  
  B2) ptetPCP_eucDist_plots.R
    - Goal: Convert euclidean distances to Protein Similarity Scores (PSSs), measure skew of PSS distributions, plots
    - Package: tidyverse
    - Plots: Histograms (all PSS Distributions; average PSS Distribution by Compartment; non-parametric skew)
    - RData: ptetPCP_PSS.RData
  
  
  C) ptetPCP_WGD.R
    - Goal: Compare PSSs between all WGD1/2/3 ohnologs
    - Package: tidyverse; seqinr
    - Custom R Functions: WGDpairs.R; WGDtestMakeBind.R
    - RData: ptetPCP_WGD-ProteinSimilarity.RData
    - Output Files: ptetPCP_WGD.R; ptetPCP_WGD_mostDivergent.csv; GeneA_GeneB.fasta (protein Fasta for all WGD ohnologs); ptetPCP_WGD_seqPSS.csv
    - External: align fasta files with muscle, measure P-distance with megacc
    - Plots: Regression (PSS vs P-Distance); Residual Plots 
