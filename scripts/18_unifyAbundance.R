# alorenzetti 202008

# description #####
# this script is intended to
# unify and parse tables generated
# on the previous scripts

# unifying protein counts and rna counts ####
# this includes transposase protein entries not having valid
# measures giving the value 1 instead of NA (low abundance)
abund = left_join(tpm, spectroWide, by = "locus_tag") %>% 
  filter(str_detect(string = locus_tag,
                    pattern = "VNG_t|VNG_s|VNG_r",
                    negate = T))

# these objects are going to be used by the heat map script
# we've got to provide a pseudocount for proteins not
# detected by the proteome survey
mobilome = dictFunCat %>% 
  filter(str_detect(cog_category, "Mobilome")) %>% 
  pull(pfeiLocusTag)

abundAlt = abund %>% 
  mutate(across(.cols = starts_with("se_abundance_protein"),
                .fns = ~ case_when((locus_tag %in% locus_tag) & # replace 2nd locus_tag for mobilome if required only for mobilome proteins
                                     is.na(.x) ~ 0,
                                   TRUE ~ as.numeric(.x)))) %>% 
  mutate(across(.cols = starts_with("mean_abundance_protein"),
                .fns = ~ case_when((locus_tag %in% locus_tag) & # replace 2nd locus_tag for mobilome if required only for mobilome proteins
                                     is.na(.x) ~ 1,
                                   TRUE ~ as.numeric(.x))))

# creating derived variables from original
# abundance variables
abundDer = abund %>% 
  mutate(mean_abundance_rna_occupancy_TP1 = mean_abundance_rna_ribofraction_TP1 / mean_abundance_rna_total_TP1,
         mean_abundance_rna_occupancy_TP2 = mean_abundance_rna_ribofraction_TP2 / mean_abundance_rna_total_TP2,
         mean_abundance_rna_occupancy_TP3 = mean_abundance_rna_ribofraction_TP3 / mean_abundance_rna_total_TP3,
         mean_abundance_rna_occupancy_TP4 = mean_abundance_rna_ribofraction_TP4 / mean_abundance_rna_total_TP4) %>%
  
  mutate(mean_abundance_rna_psiTE_TP1 = mean_abundance_protein_lysate_TP1 / mean_abundance_rna_total_TP1,
         mean_abundance_rna_psiTE_TP2 = mean_abundance_protein_lysate_TP2 / mean_abundance_rna_total_TP2,
         mean_abundance_rna_psiTE_TP3 = mean_abundance_protein_lysate_TP3 / mean_abundance_rna_total_TP3,
         mean_abundance_rna_psiTE_TP4 = mean_abundance_protein_lysate_TP4 / mean_abundance_rna_total_TP4) %>%
  
  mutate(mean_abundance_rna_tlr_TP1 = mean_abundance_rna_psiTE_TP1 / mean_abundance_rna_occupancy_TP1,
         mean_abundance_rna_tlr_TP2 = mean_abundance_rna_psiTE_TP2 / mean_abundance_rna_occupancy_TP2,
         mean_abundance_rna_tlr_TP3 = mean_abundance_rna_psiTE_TP3 / mean_abundance_rna_occupancy_TP3,
         mean_abundance_rna_tlr_TP4 = mean_abundance_rna_psiTE_TP4 / mean_abundance_rna_occupancy_TP4)

# creating a quantile normalized version of
# abund dataset
abundmean = abund %>%
  dplyr::select(starts_with("mean")) %>%
  as.matrix() %>% 
  normalize.quantiles() %>% 
  as_tibble(rownames = NULL)
colnames(abundmean) = abund %>%
  dplyr::select(starts_with("mean")) %>% 
  colnames()

abundse = abund %>%
  dplyr::select(starts_with("se")) %>%
  as.matrix() %>% 
  normalize.quantiles() %>% 
  as_tibble(rownames = NULL)
colnames(abundse) = abund %>%
  dplyr::select(starts_with("se")) %>% 
  colnames()

abundNorm = bind_cols(abund[,"locus_tag"],
                      abundmean,
                      abundse)

abundNormDer = bind_cols(abund[,"locus_tag"],
                         abundmean,
                         abundse) %>% 
  mutate(mean_abundance_rna_occupancy_TP1 = mean_abundance_rna_ribofraction_TP1 / mean_abundance_rna_total_TP1,
         mean_abundance_rna_occupancy_TP2 = mean_abundance_rna_ribofraction_TP2 / mean_abundance_rna_total_TP2,
         mean_abundance_rna_occupancy_TP3 = mean_abundance_rna_ribofraction_TP3 / mean_abundance_rna_total_TP3,
         mean_abundance_rna_occupancy_TP4 = mean_abundance_rna_ribofraction_TP4 / mean_abundance_rna_total_TP4) %>%
  
  mutate(mean_abundance_rna_psiTE_TP1 = mean_abundance_protein_lysate_TP1 / mean_abundance_rna_total_TP1,
         mean_abundance_rna_psiTE_TP2 = mean_abundance_protein_lysate_TP2 / mean_abundance_rna_total_TP2,
         mean_abundance_rna_psiTE_TP3 = mean_abundance_protein_lysate_TP3 / mean_abundance_rna_total_TP3,
         mean_abundance_rna_psiTE_TP4 = mean_abundance_protein_lysate_TP4 / mean_abundance_rna_total_TP4) %>%
  
  mutate(mean_abundance_rna_tlr_TP1 = mean_abundance_rna_psiTE_TP1 / mean_abundance_rna_occupancy_TP1,
         mean_abundance_rna_tlr_TP2 = mean_abundance_rna_psiTE_TP2 / mean_abundance_rna_occupancy_TP2,
         mean_abundance_rna_tlr_TP3 = mean_abundance_rna_psiTE_TP3 / mean_abundance_rna_occupancy_TP3,
         mean_abundance_rna_tlr_TP4 = mean_abundance_rna_psiTE_TP4 / mean_abundance_rna_occupancy_TP4)

# pivoting a version of abund non normalized dataset
abundLong = abund %>% 
  pivot_longer(cols = contains("abundance"),
               names_pattern = "^(.*)_.*_(.*_.*)_(.*)$",
               names_to = c("measure", "libtype", "timepoint"),
               values_to = "abundance") %>% 
  pivot_wider(names_from = c("measure"),
              values_from = "abundance")

# pivoting a version of abund dataset with derived vars
abundDerLong = abundDer %>% 
  mutate(mean_abundance_rna_occupancy_TP0 = 1,
         mean_abundance_rna_psiTE_TP0 = 1,
         mean_abundance_rna_tlr_TP0 = 1) %>%
  pivot_longer(cols = contains("abundance"),
               names_pattern = "^(.*)_.*_(.*_.*)_(.*)$",
               names_to = c("measure", "libtype", "timepoint"),
               values_to = "abundance") %>% 
  pivot_wider(names_from = c("measure"),
              values_from = "abundance")

# pivoting a version of abund normalized dataset
abundNormLong = abundNorm %>% 
  mutate(mean_abundance_protein_lysate_TP0 = 0,
         mean_abundance_rna_total_TP0 = 0,
         mean_abundance_rna_ribofraction_TP0 = 0) %>%
  pivot_longer(cols = contains("abundance"),
               names_pattern = "^(.*)_.*_(.*_.*)_(.*)$",
               names_to = c("measure", "libtype", "timepoint"),
               values_to = "abundance") %>% 
  pivot_wider(names_from = c("measure"),
              values_from = "abundance")

# pivoting a version of abund normalized dataset with derived vars
abundNormDerLong = abundNormDer %>% 
  mutate(mean_abundance_rna_occupancy_TP0 = 1,
         mean_abundance_rna_psiTE_TP0 = 1,
         mean_abundance_rna_tlr_TP0 = 1) %>%
  pivot_longer(cols = contains("abundance"),
               names_pattern = "^(.*)_.*_(.*_.*)_(.*)$",
               names_to = c("measure", "libtype", "timepoint"),
               values_to = "abundance") %>% 
  pivot_wider(names_from = c("measure"),
              values_from = "abundance")

# getting genes pertaining to the
# top mRNA stratum and bottom portein stratum
# for each time point
ptgsAbund = list()

# setting up timepoints
tp = c("TP1", "TP2", "TP3", "TP4")

for(i in tp){
  upperthr = 0.8
  bottomthr = 0.2
  
  protcolname = paste0("mean_abundance_protein_lysate_", i)
  mrnacolname = paste0("mean_abundance_rna_total_", i)
  
  prottopthr = abund[, protcolname] %>% quantile(probs = upperthr, na.rm = T)
  protbotthr = abund[, protcolname] %>% quantile(probs = bottomthr, na.rm = T)
    
  mrnatopthr = abund[, mrnacolname] %>% quantile(probs = upperthr, na.rm = T)
  mrnabotthr = abund[, mrnacolname] %>% quantile(probs = bottomthr, na.rm = T)
  
  ptgsAbund[[i]][["prot_top"]] = abund %>%
    filter(abund[, protcolname] >= prottopthr) %>% 
    pull(locus_tag)
  
  ptgsAbund[[i]][["prot_bot"]] = abund %>%
    filter(abund[, protcolname] < protbotthr) %>% 
    pull(locus_tag)
  
  ptgsAbund[[i]][["mrna_top"]] = abund %>%
    filter(abund[, mrnacolname] >= mrnatopthr) %>% 
    pull(locus_tag)
  
  ptgsAbund[[i]][["mrna_bot"]] = abund %>%
    filter(abund[, mrnacolname] < mrnabotthr) %>% 
    pull(locus_tag)
  
  # using an alternative version of abund
  # to get info about non observed proteins
  # that have mRNAs in the upper stratum
  ptgsAbund[[i]][["prot_non"]] = abundAlt %>%
    filter(abundAlt[, protcolname] == 1) %>% 
    pull(locus_tag)
  
  ptgsAbund[[i]][["prot_bot_mrna_top"]] = base::intersect(ptgsAbund[[i]][["prot_bot"]], ptgsAbund[[i]][["mrna_top"]])
  ptgsAbund[[i]][["prot_top_mrna_bot"]] = base::intersect(ptgsAbund[[i]][["prot_top"]], ptgsAbund[[i]][["mrna_bot"]])
  ptgsAbund[[i]][["prot_non_mrna_top"]] = base::intersect(ptgsAbund[[i]][["prot_non"]], ptgsAbund[[i]][["mrna_top"]])
  ptgsAbund[[i]][["prot_bot_prot_non_mrna_top"]] = c(ptgsAbund[[i]][["prot_bot_mrna_top"]],
                                                     ptgsAbund[[i]][["prot_non_mrna_top"]]) %>% 
    unique()
}

# getting the union of prot_bot_mrna_top
ptgsAbund$union$prot_bot_mrna_top = c(ptgsAbund$TP1$prot_bot_mrna_top,
                                                 ptgsAbund$TP2$prot_bot_mrna_top,
                                                 ptgsAbund$TP3$prot_bot_mrna_top,
                                                 ptgsAbund$TP4$prot_bot_mrna_top) %>%
  unique()

# getting the union of prot_non_mrna_top
ptgsAbund$union$prot_non_mrna_top = c(ptgsAbund$TP1$prot_non_mrna_top,
                                                 ptgsAbund$TP2$prot_non_mrna_top,
                                                 ptgsAbund$TP3$prot_non_mrna_top,
                                                 ptgsAbund$TP4$prot_non_mrna_top) %>%
  unique()

# getting the union of prot_bot_prot_non_mrna_top
ptgsAbund$union$prot_bot_prot_non_mrna_top = c(ptgsAbund$union$prot_bot_mrna_top,
                                               ptgsAbund$union$prot_non_mrna_top) %>% 
  unique()

# getting the union of abundance and fold change approaches
# for putative post-transcriptionally repressed transcripts
ptgsAbundFcQ4 = c(ptgsAbund$union$prot_bot_prot_non_mrna_top,
                  ptgs$union$Q4$short_and_long_transition_locus_tag) %>% 
  unique()
