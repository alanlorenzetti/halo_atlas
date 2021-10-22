# alorenzetti 202008

# description ####
# this script will get and parse KEGG data
# for Halobacterium salinarum
# using an API

# loading libraries ####
source("scripts/loadingLibs.R")

# processing starts ####
# getting kegg info
if(file.exists("data/keggSetHalDf.RData")){
  load("data/keggSetHalDf.RData")
  
} else {
  keggSetHal = kegg.gsets(species = "hal", id.type = "kegg")
  keggSetHal = keggSetHal$kg.sets
  
  keggSetHalDf = keggSetHal %>%
    unlist %>%
    tibble::enframe()
  
  keggSetHalDf$name = sub(keggSetHalDf$name, pattern = "[0-9]{1,}$", replacement = "")
  colnames(keggSetHalDf) = c("KEGGpathway", "locus_tag")
  
  keggSetHalDf %<>% 
    dplyr::group_by(locus_tag) %>%
    dplyr::summarise(KEGGpathway = base::paste(KEGGpathway, collapse = "; "))
  
  # this entry had gene name instead of locus_tag
  # so I am manually replacing it according to UniProt VNG
  keggSetHalDf[keggSetHalDf$locus_tag == "tbpD","locus_tag"] = "VNG_5163G"
  
  # removing all those entries not starting with VNG
  # except for the above case, all of those are represented
  # in other VNG
  keggSetHalDf = keggSetHalDf[grepl(pattern = "^VNG", keggSetHalDf$locus_tag),]
  
  # also removing underscore
  keggSetHalDf = keggSetHalDf %>% mutate(locus_tag = sub("VNG_", "VNG", locus_tag))
  
  save(keggSetHalDf, file="data/keggSetHalDf.RData")
}