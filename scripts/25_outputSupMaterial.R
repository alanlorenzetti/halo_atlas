# alorenzetti 202109

# description ####
# this script will gather and parse
# relevant info in order to output
# supp material spreadsheets

# getting started ####
# creating an output directory
if(!dir.exists("results/supp_tables")){dir.create("results/supp_tables", recursive = T)}

# summary of post-transcriptional related elements for coding genes ####
# table containing general counts of asRNAs, TPS, SmAP1 binding, RNase differential expression
summaryTPelements = list()

# 2099 mutant
summaryTPelements$`2099`$up = res2099 %>% 
  filter(logFC >= log2fcthreshold & adj.P.Val < padjthreshold) %>% 
  filter(str_detect(representative, "VNG_[0-9]")) %>% 
  pull(representative)

summaryTPelements$`2099`$down = res2099 %>% 
  filter(logFC <= -log2fcthreshold & adj.P.Val < padjthreshold) %>% 
  filter(str_detect(representative, "VNG_[0-9]")) %>% 
  pull(representative)

summaryTPelements$SmAP1 = dictFunCat %>% 
  filter(lsmSense == "yes") %>% 
  filter(str_detect(pfeiLocusTag, "VNG_[0-9]")) %>% 
  pull(pfeiLocusTag)

summaryTPelements$asRNA = dictFunCat %>% 
  filter(asRNA == "yes") %>% 
  filter(str_detect(pfeiLocusTag, "VNG_[0-9]")) %>% 
  pull(pfeiLocusTag)

summaryTPelements$tps$tps1 = tpscount %>%
  filter(tps == 1) %>% 
  filter(str_detect(representative, "VNG_[0-9]")) %>% 
  pull(representative)

summaryTPelements$tps$tps2to5 = tpscount %>%
  filter(tps > 1 & tps <= 5) %>% 
  filter(str_detect(representative, "VNG_[0-9]")) %>% 
  pull(representative)

summaryTPelements$tps$tps5 = tpscount %>%
  filter(tps > 5) %>% 
  filter(str_detect(representative, "VNG_[0-9]")) %>% 
  pull(representative)

tablePTelements = tibble(feature = c("2099 Upregulated",
                                     "2099 Downregulated",
                                     "SmAP1",
                                     "asRNA",
                                     "TPS (1)",
                                     "TPS (2-5)",
                                     "TPS (>5)"),
                         count = c(summaryTPelements$`2099`$up %>% length(),
                                   summaryTPelements$`2099`$down %>% length(),
                                   summaryTPelements$SmAP1 %>% length(),
                                   summaryTPelements$asRNA %>% length(),
                                   summaryTPelements$tps$tps1 %>% length(),
                                   summaryTPelements$tps$tps2to5 %>% length(),
                                   summaryTPelements$tps$tps5 %>% length()))

# all the atlas data as a supplemental table ####
# description of whole halo atlas dataset columns
atlasDesc = tibble(column_name = colnames(atlasdata$atlas),
                   description = c("locus_tag" = "Locus tag of a given instance according to Pfeiffer et al. (2019). This locus tag may be a representative of many others if they were collapsed in our non-redundant transcriptome.",
                                   structure(names = paste0("mean_abundance_protein_lysate_", paste0("TP", 1:4)), .Data = paste0("Proteome quantitative measure for ", paste0("TP", 1:4), ".")),
                                   structure(names = paste0("mean_abundance_rna_total_", paste0("TP", 1:4)), .Data = paste0("Transcriptome quantitative measure for ", paste0("TP", 1:4), ".")),
                                   structure(names = paste0("mean_abundance_rna_ribofraction_", paste0("TP", 1:4)), .Data = paste0("Ribo-Seq quantitative measure for ", paste0("TP", 1:4), ".")),
                                   structure(names = paste0("TE_", paste0("TP", 1:4)), .Data = paste0("Translational efficiency for ", paste0("TP", 1:4), ". Given by mean_abundance_protein_lysate divided by mean_abundance_rna_total for each time point sample.")),
                                   structure(names = paste0("RO_", paste0("TP", 1:4)), .Data = paste0("Ribosome occupancy for ", paste0("TP", 1:4), ". Given by mean_abundance_rna_ribofraction divided by mean_abundance_rna_total for each time point sample.")),
                                   "product" = "Protein product given by Pfeiffer et al. (2019).",
                                   "gene_symbol" = "Gene symbol according to COG 2020.",
                                   "cog_id" = "ID according to COG 2020.",
                                   "cog_name" = "Protein product according to COG 2020.",
                                   "cog_category" = "Category according to COG 2020.",
                                   "functional_pathway" = "Functional pathway according to COG 2020.",
                                   "lsmSense" = "Whether there is at least one SmAP1 (formerly known as LSm) binding site on the same strand of a given gene.",
                                   "lsmAntiSense" = "Whether there is at least one SmAP1 (formerly known as LSm) binding site on the opposite strand of a given gene.",
                                   "asRNA" = "Whether there is at least one annotated antisense RNA (asRNA) according to de Almeida et al. (2019).",
                                   "GC" = "GC content of a given gene.",
                                   "GCdev" = "Difference between GC content of a given gene and the mean GC content considering all the genes.",
                                   "HL" = "Experimentally determined half-life of a transcript according to Hundt et al. (2007).",
                                   "cai" = "Codon adaptation index computed using as reference set the 5% most abundant proteins in this study.",
                                   structure(names = paste0("ChIPSeq_", paste0("Tfb", c("B", "D", "G"))), .Data = paste0("Whether there is at least one ", paste0("Tfb", c("B", "D", "G")), " binding site 150 nt upstream or downstream to the first codon of a gene. Binding sites were inferred using ChIP-Seq.")),
                                   structure(names = paste0("ChIPChip_", paste0("Tbp", c("B", "C", "E", "F"))), .Data = paste0("Whether there is at least one ", paste0("Tbp", c("B", "C", "E", "F")), " binding site 150 nt upstream or downstream to the first codon of a gene. Binding sites were inferred using ChIP-Chip.")),
                                   structure(names = paste0("ChIPChip_", paste0("Tfb", c("A", "B", "C", "D", "E", "F", "G"))), .Data = paste0("Whether there is at least one ", paste0("Tfb", c("A", "B", "C", "D", "E", "F", "G")), " binding site 150 nt upstream or downstream to the first codon of a gene. Binding sites were inferred using ChIP-Chip.")),
                                   "Chromosome" = "The represented instance is located within the chromosome (NC_002607.1).",
                                   "pNRC100" = "The representative instance is located within the plasmid pNRC100 (NC_001869.1).",
                                   "pNRC200" = "The representative instance is located within the plasmid pNRC200 (NC_002608.1).",
                                   "tps" = "Whether there is at least one transcript processing site (TPS) on the same strand of a given gene.",
                                   "IRs" = "Number of inverted repeats found within a given gene.",
                                   "lfc2099" = "Log2-transformed fold change of RNase VNG2099 knockout vs. the control strain.",
                                   "ISFamily" = "If the representative instance is located within an insertion sequence, this field will display the insertion sequence family.",
                                   "protein_coding" = "Whether the representative instance encodes a protein."))

# also with nrtx table
atlasdata = list("atlas" = hmaFuncat,
                 "atlasColumnDescription" = atlasDesc,
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
pvalthr = 0.05

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
testVars = c("HL", "cai", "lfc2099")

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
        meanx = mean(xtest, na.rm = T)
        meany = mean(ytest, na.rm = T)
        if(meanx > meany){
          lvl = paste0("higher", ":", round(meanx, digits = 4), ":", round(meany, digits = 4))
        }else{
          lvl = paste0("lower", ":", round(meanx, digits = 4), ":", round(meany, digits = 4))
        }
        
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
testVars = c("HL", "cai", "lfc2099")

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
        meanx = mean(xtest, na.rm = T)
        meany = mean(ytest, na.rm = T)
        if(meanx > meany){
          lvl = paste0("higher", ":", round(meanx, digits = 4), ":", round(meany, digits = 4))
        }else{
          lvl = paste0("lower", ":", round(meanx, digits = 4), ":", round(meany, digits = 4))
        }
        
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

