# alorenzetti 202109

# description ####
# this script will gather and parse
# relevant info in order to output
# supp material spreadsheets

# getting started ####
# creating an output directory
if(!dir.exists("results/supp_tables")){dir.create("results/supp_tables", recursive = T)}

# all the atlas data as a supplemental table ####
# also with nrtx table
atlasdata = list("atlas" = hmaFuncat,
                 "non_redundant_tx_dictionary" = nrtx,
                 "non_normalized_atlas_tidy_data" = abundDerLongFuncat)

write.xlsx(atlasdata,
           file = "results/supp_tables/Table_S_atlasData.xlsx",
           overwrite = T)

# generating the table containing clusters extracted from fold change analysis ####
# creating a list of data frames to store locus tags
outfcclusters = list()

# getting data
for(i in names(ptgs)){
  spdsheet = tibble()
  
  for(j in names(ptgs[[i]])){
    line = tibble(cluster = j,
                  locus_tag = paste0(ptgs[[i]][[j]]$locus_tag, collapse = ","))
    spdsheet = bind_rows(spdsheet, line)
  }
  
  outfcclusters[[i]] = spdsheet %>% 
    mutate(color = case_when(cluster == "Q1" ~ "Red",
                             cluster == "Q2" ~ "Orange",
                             cluster == "Q3" ~ "Blue",
                             cluster == "Q4" ~ "Green",
                             cluster == "BL12" ~ "Pink",
                             cluster == "BL23" ~ "Light teal",
                             cluster == "BL34" ~ "Brown",
                             cluster == "BL41" ~ "Purple")) %>% 
    select(cluster,
           color,
           locus_tag)
}

# ordering according to panel of plots
outfcclusters = outfcclusters[c("TP2_vs_TP1",
                                "TP3_vs_TP2",
                                "TP4_vs_TP3",
                                "TP3_vs_TP1",
                                "TP4_vs_TP1")]

# writing table with multi sheets
write.xlsx(outfcclusters,
           file = "results/supp_tables/Table_S_clusters.xlsx",
           overwrite = T)

# creating a table of potentially post-transcriptionally repressed genes ####
# with enrichment p-values
# abundance-based
timepoints = names(ptgsAbund)
modus = c("prot_bot_mrna_top", "prot_non_mrna_top", "prot_bot_prot_non_mrna_top")
testcols = c("lsmSense", "asRNA", "tps")

ptgsAbundRepTib = tibble()
for(i in timepoints){
  for(j in modus){
    ptgsAbundRepTib = bind_rows(ptgsAbundRepTib,
                                tibble(criteria = j,
                                       timepoint = i,
                                       count = ptgsAbund[[i]][[j]] %>% length(),
                                       smap1_enrich_p = enrichAnalysis(hmaFuncat, ptgsAbund[[i]][[j]], "lsmSense", "yes"),
                                       asrna_enrich_p = enrichAnalysis(hmaFuncat, ptgsAbund[[i]][[j]], "asRNA", "yes"),
                                       tps_enrich_p = enrichAnalysis(hmaFuncat, ptgsAbund[[i]][[j]], "tps", "yes"),
                                       locus_tag = ptgsAbund[[i]][[j]] %>% paste0(collapse= ",")))
  }
}

ptgsAbundRepTib = ptgsAbundRepTib %>% 
  arrange(criteria)

# fc-based
timepoints = names(ptgs)
modus = names(ptgs$union$Q4)

ptgsFcRepTib = tibble()
for(i in timepoints){
    ptgsFcRepTib = bind_rows(ptgsFcRepTib,
                             tibble(criteria = "prot_downregulated_mrna_upregulated",
                                    timepoint = i,
                                    count = ptgs[[i]][["Q4"]][["locus_tag"]] %>% length(),
                                    smap1_enrich_p = enrichAnalysis(hmaFuncat, ptgs[[i]][["Q4"]][["locus_tag"]], "lsmSense", "yes"),
                                    asrna_enrich_p = enrichAnalysis(hmaFuncat, ptgs[[i]][["Q4"]][["locus_tag"]], "asRNA", "yes"),
                                    tps_enrich_p = enrichAnalysis(hmaFuncat, ptgs[[i]][["Q4"]][["locus_tag"]], "tps", "yes"),
                                    locus_tag = ptgs[[i]][["Q4"]][["locus_tag"]] %>% paste0(collapse= ",")))
}

# finding potential mechanisms ruling the post-transcriptional regulation ####

# abundance approach
# categorical variables
abundStats = tibble()
criteria = c("prot_bot_mrna_top", "prot_non_mrna_top", "prot_bot_prot_non_mrna_top", "prot_top_mrna_bot")
testVars = c("lsmSense", "asRNA", "tps")
testLevels = c("yes")
minCount = 3
pvalthr = 0.01

for(tp in names(ptgsAbund)){
  for(crit in criteria){
    worthClasses = hmaFuncat %>%
      filter(hmaFuncat$locus_tag %in% ptgsAbund[[tp]][[crit]]) %>%
      select(cog_category) %>%
      group_by(cog_category) %>%
      summarise(count = n()) %>%
      filter(count > minCount) %>% 
      pull(cog_category)
    
    for(class in c("all", worthClasses)){
      if(class == "all") {
        loci = ptgsAbund[[tp]][[crit]]
      } else {
        loci = hmaFuncat %>%
          filter(locus_tag %in% ptgsAbund[[tp]][[crit]] & cog_category == class) %>% 
          pull(locus_tag)
      }
      
      if(length(loci) == 0){next}
    
      for(var in testVars){
        for(lvl in testLevels){
          p = enrichAnalysis(df = hmaFuncat,
                             subvec = loci,
                             testcol = var,
                             lev = lvl)
          
          tibLine = tibble(criteria = crit,
                           timepoint = tp,
                           count = length(loci),
                           functional_class = class,
                           var_type = "categorical",
                           tested_var = var,
                           tested_level = lvl,
                           pvalue = p,
                           locus_tag = paste0(loci, collapse = ","))
          
          abundStats = bind_rows(abundStats,
                                 tibLine)
        }
      }
    }
  }
}
# abundStats$qvalue = p.adjust(p = abundStats$pvalue, method = "BH")

# quantitative variables
testVars = c("HL", "cai", "lfc2647", "lfc2099")

for(tp in names(ptgsAbund)){
  for(crit in criteria){
    worthClasses = hmaFuncat %>%
      filter(hmaFuncat$locus_tag %in% ptgsAbund[[tp]][[crit]]) %>%
      select(cog_category) %>%
      group_by(cog_category) %>%
      summarise(count = n()) %>%
      filter(count > minCount) %>% 
      pull(cog_category)
    
    for(class in c("all", worthClasses)){
      if(class == "all") {
        loci = ptgsAbund[[tp]][[crit]]
      } else {
        loci = hmaFuncat %>%
          filter(locus_tag %in% ptgsAbund[[tp]][[crit]] & cog_category == class) %>% 
          pull(locus_tag)
      }
      
      if(length(loci) == 0){next}
      
      for(var in testVars){
        xtest = hmaFuncat %>%
          filter(locus_tag %in% loci) %>% 
          pull(!!var) %>% 
          unname()
        ytest = hmaFuncat %>%
          filter(! (locus_tag %in% loci)) %>%
          pull(!!var) %>% 
          unname()
        
        if(is.na(mean(xtest, na.rm = T)) | is.na(mean(ytest, na.rm = T))){next}
        if(mean(xtest, na.rm = T) > mean(ytest, na.rm = T)){lvl = "higher"}else{lvl = "lower"}
        
        wstat = wilcox.test(x = xtest,
                            y = ytest,
                            paired = F)
        
        tibLine = tibble(criteria = crit,
                         timepoint = tp,
                         count = length(loci),
                         functional_class = class,
                         var_type = "quantitative",
                         tested_var = var,
                         tested_level = lvl,
                         pvalue = wstat$p.value,
                         locus_tag = paste0(loci, collapse = ","))
        
        abundStats = bind_rows(abundStats,
                               tibLine)
      }
    }
  }
}

# learning potential post-transcriptional regulation
# mechanisms for each group
ptgsAbundRepTib = abundStats %>%
  filter(pvalue < pvalthr) %>% 
  mutate(tested_var = case_when(tested_var == "lsmSense" ~ "SmAP1",
                                tested_var == "tps" ~ "TPS",
                                tested_var == "cai" ~ "CAI",
                                TRUE ~ tested_var)) %>% 
  mutate(varlvl = paste0(tested_var, ":", tested_level)) %>% 
  group_by(criteria,timepoint,count,functional_class,locus_tag) %>%
  summarise(mechanism = paste0(varlvl, collapse = ","),
            pvalues = paste0(pvalue, collapse = ",")) %>% 
  select(criteria,
         timepoint,
         count,
         functional_class,
         mechanism,
         pvalues,
         locus_tag)

# performing the same operations for fold change
# based analysis

# abundance approach
# categorical variables
fcStats = tibble()
criteria = c("Q1", "Q2", "Q3", "Q4", "BL41")
testVars = c("lsmSense", "asRNA", "tps")
testLevels = c("yes")
minCount = 3

for(tp in names(ptgs)){
  for(crit in criteria){
    worthClasses = hmaFuncat %>%
      filter(hmaFuncat$locus_tag %in% ptgs[[tp]][[crit]]$locus_tag) %>%
      select(cog_category) %>%
      group_by(cog_category) %>%
      summarise(count = n()) %>%
      filter(count > minCount) %>% 
      pull(cog_category)
    
    for(class in c("all", worthClasses)){
      if(class == "all") {
        loci = ptgs[[tp]][[crit]]$locus_tag
      } else {
        loci = hmaFuncat %>%
          filter(locus_tag %in% ptgs[[tp]][[crit]]$locus_tag & cog_category == class) %>% 
          pull(locus_tag)
      }
      
      if(length(loci) == 0){next}
      
      for(var in testVars){
        for(lvl in testLevels){
          p = enrichAnalysis(df = hmaFuncat,
                             subvec = loci,
                             testcol = var,
                             lev = lvl)
          
          tibLine = tibble(criteria = crit,
                           timepoint = tp,
                           count = length(loci),
                           functional_class = class,
                           var_type = "categorical",
                           tested_var = var,
                           tested_level = lvl,
                           pvalue = p,
                           locus_tag = paste0(loci, collapse = ","))
          
          fcStats = bind_rows(fcStats,
                              tibLine)
        }
      }
    }
  }
}

# quantitative variables
testVars = c("HL", "cai", "lfc2647", "lfc2099")

for(tp in names(ptgs)){
  for(crit in criteria){
    worthClasses = hmaFuncat %>%
      filter(hmaFuncat$locus_tag %in% ptgs[[tp]][[crit]]$locus_tag) %>%
      select(cog_category) %>%
      group_by(cog_category) %>%
      summarise(count = n()) %>%
      filter(count > minCount) %>% 
      pull(cog_category)
    
    for(class in c("all", worthClasses)){
      if(class == "all") {
        loci = ptgs[[tp]][[crit]]$locus_tag
      } else {
        loci = hmaFuncat %>%
          filter(locus_tag %in% ptgsAbund[[tp]][[crit]]$locus_tag & cog_category == class) %>% 
          pull(locus_tag)
      }
      
      if(length(loci) == 0){next}
      
      for(var in testVars){
        xtest = hmaFuncat %>%
          filter(locus_tag %in% loci) %>% 
          pull(!!var) %>% 
          unname()
        ytest = hmaFuncat %>%
          filter(! (locus_tag %in% loci)) %>%
          pull(!!var) %>% 
          unname()
        
        if(is.na(mean(xtest, na.rm = T)) | is.na(mean(ytest, na.rm = T))){next}
        if(mean(xtest, na.rm = T) > mean(ytest, na.rm = T)){lvl = "higher"}else{lvl = "lower"}
        
        wstat = wilcox.test(x = xtest,
                            y = ytest,
                            paired = F)
        
        tibLine = tibble(criteria = crit,
                         timepoint = tp,
                         count = length(loci),
                         functional_class = class,
                         var_type = "quantitative",
                         tested_var = var,
                         tested_level = lvl,
                         pvalue = wstat$p.value,
                         locus_tag = paste0(loci, collapse = ","))
        
        fcStats = bind_rows(fcStats,
                            tibLine)
      }
    }
  }
}

# learning potential post-transcriptional regulation
# mechanisms for each group
ptgsFcRepTib = fcStats %>%
  filter(pvalue < pvalthr) %>% 
  mutate(tested_var = case_when(tested_var == "lsmSense" ~ "SmAP1",
                                tested_var == "tps" ~ "TPS",
                                tested_var == "cai" ~ "CAI",
                                TRUE ~ tested_var)) %>% 
  mutate(varlvl = paste0(tested_var, ":", tested_level)) %>% 
  group_by(criteria,timepoint,count,functional_class,locus_tag) %>%
  summarise(mechanism = paste0(varlvl, collapse = ","),
            pvalues = paste0(pvalue, collapse = ",")) %>% 
  select(criteria,
         timepoint,
         count,
         functional_class,
         mechanism,
         pvalues,
         locus_tag)

# creating a list to store tables
ptgsEnrich = list("Abundance_based_analysis" = ptgsAbundRepTib,
                  "FC_based_analysis" = ptgsFcRepTib)

# writing table with multi sheets
write.xlsx(ptgsEnrich,
           file = "results/supp_tables/Table_S_ptgsEnrich.xlsx",
           overwrite = T)
