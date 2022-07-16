# alorenzetti 202008

# this is a scaffold script
# here you can find the adequate
# script order

# creating a few directories necessary to store results and plots
if(!dir.exists("plots")){dir.create("plots")}
if(!dir.exists("results")){dir.create("results")}

# loading libs
source("scripts/01_loadingLibs.R")

# parsers of misc features
source("scripts/02_parseHalfLives.R")
source("scripts/03_codonUsage.R") # check why there are a few proteins (<10) out of frame; possibly pseudogenes
source("scripts/04_parseSmAP1.R")
source("scripts/05_parseAntisenseRNAs.R")
source("scripts/06_parseLocation.R")
source("scripts/07_parseTPS.R")
source("scripts/08_parseAndProcess2099.R")
source("scripts/09_parseMembraneAndOtherProteins.R")
source("scripts/10_getAndParseMiscFeatures.R")
source("scripts/11_combineAndWrangleFunCat.R")

# proteome parser and analysis
source("scripts/12_parseProteomicsSpectronaut.R")
source("scripts/13_DEanalysisProtein.R")

# transcriptome parser and analysis
source("scripts/14_parseRNAtpms.R")
source("scripts/15_DEanalysisRNA.R")

# exploratory analysis plots
source("scripts/16_exploratoryAnalysis.R")

# models and functions to work with fold change data
source("scripts/17_models.R")

# working with abundance 
source("scripts/18_unifyAbundance.R")

# heatmap generation
source("scripts/19_heatmapsAllAbsolute_v2_english.R")

# enrichment functions and analyses
source("scripts/20_enrichment.R")

# other figures and analyses
# including abundance correlation plots
source("scripts/21_analysesAndFigures_v3.R") 

# supplementary material
source("scripts/22_outputSupMaterial.R")
