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
source("scripts/03_codonUsage.R") # check why there are a few proteins (<10) out of frame
source("scripts/04_parseMicrobesOnlinePredOperons.R")
source("scripts/05_parseSmAP1.R")
source("scripts/06_parseAntisenseRNAs.R")
source("scripts/07_parseTFChipSeq.R")
source("scripts/08_parseTFChipChip.R")
source("scripts/09_parseLocation.R")
source("scripts/10_parseTPS.R")
source("scripts/11_parseiTSS.R")
source("scripts/12_genAndParseIRs.R")
# source("scripts/parseAndProcess2647.R")
# source("scripts/parseAndProcess1503.R")
source("scripts/13_parseAndProcess2099.R")
source("scripts/parseMembraneProteins.R")
source("scripts/14_getAndParseMiscFeatures.R")
source("scripts/15_combineAndWrangleFunCat.R")

# proteome parser and analysis
source("scripts/16_parseProteomicsSpectronaut.R")
source("scripts/17_DEanalysisProtein.R")

# transcriptome parser and analysis
source("scripts/18_parseRNAtpms.R")
source("scripts/19_DEanalysisRNA.R")

# exploratory analysis plots
source("scripts/20_exploratoryAnalysis.R")

# models and functions to work with fold change data
source("scripts/21_models.R")

# working with abundance 
source("scripts/22_unifyAbundance.R")

# heatmap generation
source("scripts/23_heatmapsAllAbsolute_v2_english.R")

# enrichment functions and analyses
source("scripts/24_enrichment.R")

# other figures and analyses
# including abundance correlation plots
source("scripts/25_analysesAndFigures_v3.R") 

# supplementary material
source("scripts/26_outputSupMaterial.R")
