# alorenzetti 202008

# description ####
# this script will get and parse ArCOG data
# for Halobacterium salinarum

# loading libraries ####
source("scripts/loadingLibs.R")

# processing starts ####
# getting arcog info
if(file.exists("data/arcogHal.RData")){
  load("data/arcogHal.RData")
} else {
  # reading arcog
  # I disabled quotes, otherwise it
  # would truncate the reading process
  # and now it seems to work well
  arcogNames = read_delim("ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/ar14.arCOGdef19.tab",
                          delim = "\t",
                          quote = "",
                          col_names = c("arCOG_ID",
                                        "func",
                                        "geneName",
                                        "name",
                                        "COG_ID",
                                        "pfam_ID",
                                        "cdd_ID",
                                        "TIGR_ID"))
  
  arcogClasses = read_delim("ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/funclass.tab",
                            delim = "\t",
                            quote = "",
                            col_names = c("class_ID",
                                          "class",
                                          "description"))
  # it is going to throw a warning
  # since file has heterogeneous number
  # of fields
  arcogData = read_csv("ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/ar14.arCOG.csv",
                       quote = "",
                       col_names = c("domain_ID", "genome_name", "protein_ID",
                                     "protein_length", "domain_start", "domain_end",
                                     "arCOG_ID", "membership_class",
                                     "deprecated_field1", "deprecated_field2"))
  
  # cogData could have more than one COG for a GI
  # sometimes there is even a repetitive COG for the same GI
  # due to distinct domain regions
  arcogData = arcogData %>% 
    filter(genome_name == "Halobacterium_NRC_1_uid57769") %>% 
    dplyr::select(domain_ID, arCOG_ID) %>% 
    dplyr::distinct()
  
  # selecting cols of interest
  arcogNames = arcogNames %>%
    dplyr::select(arCOG_ID, func, name)
  
  arcogFunct = left_join(arcogNames, arcogClasses, by = c("func" = "class_ID")) %>%
    dplyr::select("arCOG_ID", "func", "description", "name") %>% 
    dplyr::group_by(arCOG_ID) %>% 
    summarise(funcCode = paste0(func, collapse = "|"),
              arCOG = paste0(description, collapse = "|"),
              arCOGproduct = paste0(name, collapse = "|"))
  
  arcogFinal = left_join(arcogData, arcogFunct, by = "arCOG_ID") %>% 
    mutate(arCOG = arCOG,
           arCOGcode = funcCode) %>% 
    dplyr::select(gi = domain_ID,
                  arCOG_ID, 
                  arCOGcode,
                  arCOG,
                  arCOGproduct) %>%
    mutate(gi = as.character(gi)) %>% 
    dplyr::select(-gi) %>% 
    filter(!str_detect(arCOG, "^NA$")) %>%
    distinct()
  
  arcogHal = arcogFinal
  
  save(arcogHal, file="data/arcogHal.RData")
}
