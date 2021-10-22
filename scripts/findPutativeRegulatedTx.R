# alorenzetti 20201230

# description ####
# this script will find potentially post-transcriptionally
# and/or translationally regulated transcripts

# main function ####
# provide abund dataset and TPi (options TP2, TP3, TP4)
findPutativeRegulated = function(df, tp){
  tpi = tp %>%
    str_replace(., "^TP(.*)$", "\\1") %>% 
    as.numeric()
    
  rnatpi = paste0("mean_abundance_rna_total_TP", tpi)
  rnatpiminus1 = paste0("mean_abundance_rna_total_TP", tpi-1)
  proteintpi = paste0("mean_abundance_protein_lysate_TP", tpi)
  proteintpiminus1 = paste0("mean_abundance_protein_lysate_TP", tpi-1)
  ribotpi = paste0("mean_abundance_rna_ribofraction_TP", tpi)
  ribotpiminus1 = paste0("mean_abundance_rna_ribofraction_TP", tpi-1)
  
  # interesting cases
  who = list()
  results = list()
  
  # mRNA UP and protein FLAT
  who[["mrna_up_protein_flat"]] = (df[,rnatpi]/df[,rnatpiminus1] >= 2) &
    (df[,proteintpi]/df[,proteintpiminus1] >= 1 & df[,proteintpi]/df[,proteintpiminus1] <= 1.5)
  results[["mrna_up_protein_flat"]] = df$locus_tag[which(who[["mrna_up_protein_flat"]])]
  
  # mRNA DOWN and protein FLAT and ribo UP
  who[["mrna_down_protein_flat_ribo_up"]] = (df[,rnatpi]/df[,rnatpiminus1] <= 0.5) &
    (df[,proteintpi]/df[,proteintpiminus1] >= 1 & df[,proteintpi]/df[,proteintpiminus1] <= 1.5) &
    (df[,ribotpi]/df[,ribotpiminus1] >= 2)
  results[["mrna_down_protein_flat_ribo_up"]] = df$locus_tag[which(who[["mrna_down_protein_flat_ribo_up"]])]
  
  return(results)
    
}

# running for normalized and non-normalized datasets ####
# all possible timepoints
datasets = list()
datasets[["nonnorm"]] = abund
datasets[["norm"]] = abundNorm

output = list()

for(i in names(datasets)){
  for(j in paste0("TP", 2:4)){
    output[[i]][[j]] = findPutativeRegulated(datasets[[i]], j)
  }
}