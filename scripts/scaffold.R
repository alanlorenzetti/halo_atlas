# alorenzetti 202008

# this is a scaffold script
# here you can find the adequate
# script order

# creating a few directories necessary to store results and plots
if(!dir.exists("plots")){dir.create("plots")}
if(!dir.exists("results")){dir.create("results")}

# loading libs
source("scripts/loadingLibs.R")

# ~/gdrive/runKallisto/plotKmerTests.R
#source("scripts/locusTagDictGbk2Tpa.R") # about to get deprecated; replaced by https://alanlorenzetti.github.io/halo_nr_tx/
source("scripts/parseHalfLives.R")
source("scripts/codonUsage.R") # check why there are a few proteins (<10) out of o frame
source("scripts/parseMicrobesOnlinePredOperons.R")
source("scripts/parseAntisenseRNAs.R")
#source("scripts/getAndParseKEGG.R") # about to get deprecated
#source("scripts/getAndParseUniprotInfo.R") # about to get deprecated
#source("scripts/getAndParseArCOGs.R") # about to get deprecated
source("scripts/parseTFChipSeq.R")
source("scripts/parseTFChipChip.R")
source("scripts/parseLocation.R")
source("scripts/parseTPS.R")
source("scripts/genAndParseIRs.R")
source("scripts/parseAndProcess2647.R")
source("scripts/parseAndProcess1503.R")
source("scripts/parseAndProcess2099.R")
source("scripts/getAndParseMiscFeatures.R")
source("scripts/combineAndWrangleFunCat.R")

# this version is using unnormalized new version of spectronaut data
# I think we should get back to the previous version
source("scripts/parseProteomicsSpectronaut.R")
#source("scripts/parseProteomicsOneOmics.R") # about to get deprecated
source("scripts/DEanalysisProtein.R")

source("scripts/parseRNAtpms.R")
source("scripts/DEanalysisRNA.R")

source("scripts/exploratoryAnalysis.R") # exploratory analysis plots

source("scripts/models.R")
# source("scripts/proteinRegRules.R") # about to get deprecated
# source("scripts/timecourseAnalysis.R") # about to get deprecated

# source("scripts/heatmaps.R") about to get deprecated
source("scripts/unifyAbundance.R")
# source("scripts/heatmapsAllRelativeChanges.R") # about to get deprecated
# source("scripts/heatmapsAllAbsolute.R") # about to get deprecated
#source("scripts/heatmapsAllAbsolute_v2_tese.R")
source("scripts/heatmapsAllAbsolute_v2_english.R")
source("scripts/enrichment.R")
#source("scripts/clusterEnrichmentAnalysis.R") # about to get deprecated
# source("scripts/computeAnglesAndMag.R") # about to get deprecated
#source("scripts/globalFCanalysis.R")
#source("scripts/analysesAndFigures.R")
#source("scripts/analysisAndFiguresRequested20201009.R") # correlation plots
source("scripts/analysesAndFigures_v3.R") # correlation plots
source("scripts/outputSupMaterial.R") # correlation plots
#source("scripts/findPutativeRegulatedTx.R")
# source("scripts/norm_issues_tmp.R") # accessory and not essential
# source("scripts/bepeReportFigures.R")