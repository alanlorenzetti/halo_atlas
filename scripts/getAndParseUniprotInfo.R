# alorenzetti 202008

# description ####
# this script will get and parse Uniprot data
# for Halobacterium salinarum
# using an API

# loading libraries ####
source("scripts/loadingLibs.R")

# processing starts ####
################# WARNING ##################
# somehow eggnog data for halo has disappeared
# from Uniprot; We will have to use the previsously
# downloaded version instead of downloading it again
# other than having the file read to use, this will need 
# a MAJOR fix
############################################
# getting Uniprot info
# uniprot is not letting us get all the queried terms
# we have to split the queries by up to 99 each time
if(file.exists("data/uniprotHalDf.RData")){
  load("data/uniprotHalDf.RData")
} else {
  uniprotHal = UniProt.ws(taxId=64091)
  columns = c("UNIPROTKB", "GO", "EGGNOG")
  
  keys2get = keys(uniprotHal, keytype = "KEGG")
  
  m = keys2get %>% length()
  v = seq(from=1, to=m, by=98)
  if(v[v %>% length()]){v[v %>% length() + 1] = m}
  table=NULL
  
  for(i in seq(1, length(v)-1)){
    if(i == 1){
      chunk = UniProt.ws::select(uniprotHal,
                                 keys=keys2get[v[i]:v[i+1]],
                                 columns = columns,
                                 keytype = "KEGG")
      
      table = bind_rows(table,chunk)
      
    }else{
      chunk = UniProt.ws::select(uniprotHal,
                                 keys=keys2get[v[i]:v[i+1]],
                                 columns = columns,
                                 keytype = "KEGG")
      
      table = bind_rows(table,chunk)
    }
  }
  
  table = table %>%
    filter(str_detect(KEGG, "^hal:VNG")) %>% 
    filter(str_detect(EGGNOG, "^arCOG")) %>% 
    distinct() %>% 
    rename(locus_tag = "KEGG",
           arCOG_ID = "EGGNOG") %>% 
    mutate(locus_tag = sub("hal:VNG_", "VNG", locus_tag))
  
  uniprotHalDf = table
  
  save(uniprotHalDf, file="data/uniprotHalDf.RData")
}
