# alorenzetti 202006

# description #####
# this script will create qualitative regulation
# clusters based on status of proteins, mRNAs, and RPFs

# loading libs #####
source("scripts/loadingLibs.R")

##########################finding protein regulation status using mRNA and prot
findNumbersReg = function(df, timepointContrast){
  A = list() ; B = list()
  time = timepointContrast
  
  # setting higher and lower time point contrasts
  ihigh = sub("^(.*)_vs_.*$","\\1",time)
  ilow = sub("^.*_vs_(.*$)","\\1",time)
  
  # merging dfs to be compared
  dfhigh = df[[ihigh]] %>% rename_at(vars(contains("TP")), list(~sub("TP.$","",.)))
  dflow = df[[ilow]] %>% rename_at(vars(contains("TP")), list(~sub("TP.$","",.)))
  combdf = left_join(dfhigh, dflow, by = "locus_tag", suffix = c(paste0("_", ihigh, "_high"),
                                                                 paste0("_", ilow, "_low")))
  # setting up varnames according to desired contrasts
  vartotrna = paste0("mRNA_", c(ihigh,ilow), "_", c("high","low"))
  varribo = paste0("RPF_", c(ihigh,ilow), "_", c("high","low"))
  varprot = paste0("protein_", c(ihigh,ilow), "_", c("high","low"))
  
  varsigtotrna = paste0("sigtotrna_", c(ihigh,ilow), "_", c("high","low"))
  varsigribo = paste0("sigribo_", c(ihigh,ilow), "_", c("high","low"))
  varsigprot = paste0("sigProt_", c(ihigh,ilow), "_", c("high","low"))
  
  varborderlinetotrna = paste0("borderlineZerototrna_", c(ihigh,ilow), "_", c("high","low"))
  varborderlineribo = paste0("borderlineZeroribo_", c(ihigh,ilow), "_", c("high","low"))
  varborderlineprot = paste0("borderlineZeroProt_", c(ihigh,ilow), "_", c("high","low"))
  
  ### using only PROTEIN and mRNA
  ## within the set of UPREGULATED PROTEINS
  # upregulated mRNAs
  A[["ProteinUp_mRNAUp"]] = combdf %>%
    filter(!!as.symbol(varprot[1]) > 0 &
             !!as.symbol(varsigprot[1]) == "yes" &
             !!as.symbol(vartotrna[2]) > 0 &
             !!as.symbol(varsigtotrna[2]) == "yes") %>% 
    mutate(regRule = "ProteinUp_mRNAUp")
  
  # non changing mRNAs
  A[["ProteinUp_mRNAZero"]] = combdf %>%
    filter(!!as.symbol(varprot[1]) > 0 &
             !!as.symbol(varsigprot[1]) == "yes" &
             !!as.symbol(varborderlinetotrna[2]) == "yes") %>% 
    mutate(regRule = "ProteinUp_mRNAZero")
  
  # downregulated mRNAs
  A[["ProteinUp_mRNADown"]] = combdf %>%
    filter(!!as.symbol(varprot[1]) > 0 &
             !!as.symbol(varsigprot[1]) == "yes" &
             !!as.symbol(vartotrna[2]) < 0 &
             !!as.symbol(varsigtotrna[2]) == "yes") %>% 
    mutate(regRule = "ProteinUp_mRNADown")
  
  ## within the set of NONCHANGING PROTEINS
  # upregulated mRNAs
  A[["ProteinZero_mRNAUp"]] = combdf %>%
    filter(!!as.symbol(varborderlineprot[1]) == "yes" &
             !!as.symbol(vartotrna[2]) > 0 &
             !!as.symbol(varsigtotrna[2]) == "yes") %>% 
    mutate(regRule = "ProteinZero_mRNAUp")
  
  # non changing mRNAs
  A[["ProteinZero_mRNAZero"]] = combdf %>%
    filter(!!as.symbol(varborderlineprot[1]) == "yes" &
             !!as.symbol(varborderlinetotrna[2]) == "yes") %>% 
    mutate(regRule = "ProteinZero_mRNAZero")
  
  # downregulated mRNAs
  A[["ProteinZero_mRNADown"]] = combdf %>%
    filter(!!as.symbol(varborderlineprot[1]) == "yes" &
             !!as.symbol(vartotrna[2]) < 0 &
             !!as.symbol(varsigtotrna[2]) == "yes") %>% 
    mutate(regRule = "ProteinZero_mRNADown")
  
  ## within the set of DOWNREGULATED PROTEINS
  # upregulated mRNAs
  A[["ProteinDown_mRNAUp"]] = combdf %>%
    filter(!!as.symbol(varprot[1]) < 0 &
             !!as.symbol(varsigprot[1]) == "yes" &
             !!as.symbol(vartotrna[2]) > 0 &
             !!as.symbol(varsigtotrna[2]) == "yes") %>% 
    mutate(regRule = "ProteinDown_mRNAUp")
  
  # non changing mRNAs
  A[["ProteinDown_mRNAZero"]] = combdf %>%
    filter(!!as.symbol(varprot[1]) < 0 &
             !!as.symbol(varsigprot[1]) == "yes" &
             !!as.symbol(varborderlinetotrna[2]) == "yes") %>% 
    mutate(regRule = "ProteinDown_mRNAZero")
  
  # downregulated mRNAs
  A[["ProteinDown_mRNADown"]] = combdf %>%
    filter(!!as.symbol(varprot[1]) < 0 &
             !!as.symbol(varsigprot[1]) == "yes" &
             !!as.symbol(vartotrna[2]) < 0 &
             !!as.symbol(varsigtotrna[2]) == "yes") %>% 
    mutate(regRule = "ProteinDown_mRNADown")
  
  ### using protein, mRNA, and RPFs
  ### to find TLR genes depending on ribosome to explain it
  ## within the set of UPREGULATED proteins
  # non changing mRNAs and upregulated RPFs
  B[["ProteinUp_mRNAZero_RPFUp"]] = combdf %>% filter(!!as.symbol(varprot[1]) > 0 &
                                                        !!as.symbol(varsigprot[1]) == "yes" &
                                                        !!as.symbol(varborderlinetotrna[2]) == "yes" &
                                                        !!as.symbol(varribo[2]) > 0 &
                                                        !!as.symbol(varsigribo[2]) == "yes")
  
  # downregulated mRNAs and upregulated RPFs
  B[["ProteinUp_mRNADown_RPFUp"]] = combdf %>% filter(!!as.symbol(varprot[1]) > 0 &
                                                        !!as.symbol(varsigprot[1]) == "yes" &
                                                        !!as.symbol(vartotrna[2]) < 0 &
                                                        !!as.symbol(varsigtotrna[2]) == "yes" &
                                                        !!as.symbol(varribo[2]) > 0 &
                                                        !!as.symbol(varsigribo[2]) == "yes")
  
  ## within the set of NON CHANGING proteins
  ## whose RPFs compensates for mRNA changes
  # upregulated RPFs with a downregulated mRNA
  B[["ProteinZero_mRNADown_RPFUp"]] = combdf %>% filter(!!as.symbol(varborderlineprot[1]) == "yes" &
                                                          !!as.symbol(vartotrna[2]) < 0 &
                                                          !!as.symbol(varsigtotrna[2]) == "yes" &
                                                          !!as.symbol(varribo[2]) > 0 &
                                                          !!as.symbol(varsigribo[2]) == "yes")
  
  # downregulated RPFs with a upregulated mRNA
  B[["ProteinZero_mRNAUp_RPFDown"]] = combdf %>% filter(!!as.symbol(varborderlineprot[1]) == "yes" &
                                                          !!as.symbol(vartotrna[2]) > 0 &
                                                          !!as.symbol(varsigtotrna[2]) == "yes" &
                                                          !!as.symbol(varribo[2]) < 0 &
                                                          !!as.symbol(varsigribo[2]) == "yes")
  
  ## within the set of DOWNREGULATED proteins
  # non-changing mRNAs with downregulated RPFs
  B[["ProteinDown_mRNAZero_RPFDown"]] = combdf %>% filter(!!as.symbol(varprot[1]) < 0 &
                                                            !!as.symbol(varsigprot[1]) == "yes" &
                                                            !!as.symbol(varborderlinetotrna[2]) == "yes" &
                                                            !!as.symbol(varribo[2]) < 0 &
                                                            !!as.symbol(varsigribo[2]) == "yes")
  
  # correspond to upregulated mRNAs with downregulated RPFs
  B[["ProteinDown_mRNAUp_RPFDown"]] = combdf %>% filter(!!as.symbol(varprot[1]) < 0 &
                                                          !!as.symbol(varsigprot[1]) == "yes" &
                                                          !!as.symbol(vartotrna[2]) > 0 &
                                                          !!as.symbol(varsigtotrna[2]) == "yes" &
                                                          !!as.symbol(varribo[2]) < 0 &
                                                          !!as.symbol(varsigribo[2]) == "yes")
  
  return(list(Prot_mRNA = A,
              Prot_mRNA_RPF = B))
}

# computing for every contrast we are interested in
contrasts = c("TP3_vs_TP2",
              "TP4_vs_TP3",
              "TP2_vs_TP2",
              "TP3_vs_TP3",
              "TP4_vs_TP4")

regRules = list()
for(i in contrasts){
  regRules[[i]] = findNumbersReg(alldfs, i)
}

# function to extract the lists containing locus_tags
# and return euler plots
# for each of the time point comparison
# eulerPlot = function(bigList){
#   contrasts = c("TP3_vs_TP2",
#                 "TP4_vs_TP3",
#                 "TP2_vs_TP2",
#                 "TP3_vs_TP3",
#                 "TP4_vs_TP4")
#   
#   ltList = list()
#   for(i in contrasts){
#     ltList[[i]] = lapply(bigList[[i]][["Prot_mRNA_RPF"]], function(x) dplyr::select(x, locus_tag)) %>%
#       unlist() %>%
#       unname()
#   }
#   
#   p = plot(euler(compact(ltList)), quantities = T)
#   
#   return(p)
# }
# 
# # calling function to plot euler
# eulerPlot(regRules)
