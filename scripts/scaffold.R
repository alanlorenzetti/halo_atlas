# alorenzetti 202008

# this is a scaffold script
# here you can find the adequate
# script order

# creating a few directories necessary to store results and plots
if(!dir.exists("plots")){dir.create("plots")}
if(!dir.exists("results")){dir.create("results")}

# loading libs
source("scripts/loadingLibs.R")

# parsers of misc features
source("scripts/parseHalfLives.R")
source("scripts/codonUsage.R") # check why there are a few proteins (<10) out of frame
source("scripts/parseMicrobesOnlinePredOperons.R")
source("scripts/parseAntisenseRNAs.R")
source("scripts/parseTFChipSeq.R")
source("scripts/parseTFChipChip.R")
source("scripts/parseLocation.R")
source("scripts/parseTPS.R")
source("scripts/parseiTSS.R")
source("scripts/genAndParseIRs.R")
# source("scripts/parseAndProcess2647.R")
# source("scripts/parseAndProcess1503.R")
source("scripts/parseAndProcess2099.R")
source("scripts/getAndParseMiscFeatures.R")
source("scripts/combineAndWrangleFunCat.R")

# proteome parser and analysis
source("scripts/parseProteomicsSpectronaut.R")
source("scripts/DEanalysisProtein.R")

# transcriptome parser and analysis
source("scripts/parseRNAtpms.R")
source("scripts/DEanalysisRNA.R")

# exploratory analysis plots
source("scripts/exploratoryAnalysis.R")

# models and functions to work with fold change data
source("scripts/models.R")

# working with abundance 
source("scripts/unifyAbundance.R")

# heatmap generation
source("scripts/heatmapsAllAbsolute_v2_english.R")

# enrichment functions and analyses
source("scripts/enrichment.R")

# other figures and analyses
# including abundance correlation plots
source("scripts/analysesAndFigures_v3.R") 

# supplementary material
source("scripts/outputSupMaterial.R")
