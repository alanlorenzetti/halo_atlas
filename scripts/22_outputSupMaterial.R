# alorenzetti 202109

# description ####
# this script will gather and parse
# relevant info in order to output
# supp material spreadsheets

# getting started ####
# creating an output directory
if(!dir.exists("results/supp_tables")){dir.create("results/supp_tables", recursive = T)}

# creating a directory for figures
if(!dir.exists("results/figures")){dir.create("results/figures", recursive = T)}

# summary of post-transcriptional related elements for coding genes ####
# table containing general counts of asRNAs, TPS, SmAP1 binding,
# RNase differential expression
# main-text table
summaryTPelements = list()

# 2099 mutant
summaryTPelements$`2099`$up = res2099 %>% 
  filter(logFC >= log2fcthreshold & P.Value < padjthreshold) %>% 
  filter(str_detect(representative, "VNG_[0-9]")) %>% 
  pull(representative)

summaryTPelements$`2099`$down = res2099 %>% 
  filter(logFC <= -log2fcthreshold & P.Value < padjthreshold) %>% 
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
atlasColDescription = c("locus_tag" = "Locus tag for a given instance according to Pfeiffer et al. (2019). This locus tag may be a representative of many others if they were collapsed in our non-redundant transcriptome.",
                structure(names = paste0("mean_abundance_protein_lysate_", paste0("TP", 1:4)), .Data = paste0("Proteome quantitative measure (a pseudocount was imputed for missing values) for ", paste0("TP", 1:4), ".")),
                structure(names = paste0("mean_abundance_rna_total_", paste0("TP", 1:4)), .Data = paste0("Transcriptome quantitative measure (TPM+1) for ", paste0("TP", 1:4), ".")),
                structure(names = paste0("mean_abundance_rna_ribofraction_", paste0("TP", 1:4)), .Data = paste0("Ribo-Seq quantitative measure (TPM+1) for ", paste0("TP", 1:4), ".")),
                structure(names = paste0("TE_", paste0("TP", 1:4)), .Data = paste0("Translational efficiency for ", paste0("TP", 1:4), ". Given by mean_abundance_protein_lysate divided by mean_abundance_rna_total for each time point sample.")),
                structure(names = paste0("RO_", paste0("TP", 1:4)), .Data = paste0("Ribosome occupancy for ", paste0("TP", 1:4), ". Given by mean_abundance_rna_ribofraction divided by mean_abundance_rna_total for each time point sample.")),
                "product" = "Protein product given by Pfeiffer et al. (2019).",
                "gene_symbol" = "Gene symbol according to COG 2020.",
                "cog_id" = "ID according to COG 2020.",
                "cog_name" = "Protein product according to COG 2020.",
                "cog_category" = "Category according to COG 2020.",
                "functional_pathway" = "Functional pathway according to COG 2020.",
                "smap1Sense" = "Whether there is at least one SmAP1 binding site on the same strand for a given gene.",
                # "smap1AntiSense" = "Whether there is at least one SmAP1 binding site on the opposite strand for a given gene.",
                "asRNA" = "Whether there is at least one annotated antisense RNA (asRNA) according to de Almeida et al. (2019).",
                "tps" = "Whether there is at least one transcript processing site (TPS) on the same strand for a given gene.",
                "GC" = "GC content for a given gene.",
                "GCdev" = "Difference between GC content for a given gene and the mean GC content considering all the genes.",
                "HL" = "Experimentally determined half-life of a transcript according to Hundt et al. (2007).",
                "cai" = "Codon adaptation index computed using as reference set the 5% most abundant proteins in this study.",
                # structure(names = paste0("ChIPSeq_", paste0("Tfb", c("B", "D", "G"))), .Data = paste0("Whether there is at least one ", paste0("Tfb", c("B", "D", "G")), " binding site 150 nt upstream or downstream to the first codon of a gene. Binding sites were inferred using ChIP-Seq.")),
                # structure(names = paste0("ChIPChip_", paste0("Tbp", c("B", "C", "E", "F"))), .Data = paste0("Whether there is at least one ", paste0("Tbp", c("B", "C", "E", "F")), " binding site 150 nt upstream or downstream to the first codon of a gene. Binding sites were inferred using ChIP-Chip.")),
                # structure(names = paste0("ChIPChip_", paste0("Tfb", c("A", "B", "C", "D", "E", "F", "G"))), .Data = paste0("Whether there is at least one ", paste0("Tfb", c("A", "B", "C", "D", "E", "F", "G")), " binding site 150 nt upstream or downstream to the first codon of a gene. Binding sites were inferred using ChIP-Chip.")),
                "Chromosome" = "The represented instance is located within the chromosome (NC_002607.1).",
                "pNRC100" = "The representative instance is located within the plasmid pNRC100 (NC_001869.1).",
                "pNRC200" = "The representative instance is located within the plasmid pNRC200 (NC_002608.1).",
                # "IRs" = "Number of inverted repeats found within a given gene.",
                "lfc2099" = "Log2-transformed fold change of RNase VNG2099 knockout vs. the control strain.",
                "ISFamily" = "If the representative instance is located within an insertion sequence, this field will display the insertion sequence family."
                # "protein_coding" = "Whether the representative instance encodes a protein."
                )
                
atlasDesc = tibble(column_name = names(atlasColDescription),
                   description = atlasColDescription)


atlasTidyColDescription = c("locus_tag" = "Locus tag for a given instance according to Pfeiffer et al. (2019). This locus tag may be a representative of many others if they were collapsed in our non-redundant transcriptome.",
                            "libtype" = "Measured variable. rna_total: mRNA level measured by RNA-Seq (TPM+1); rna_ribo: ribosome protected mRNA fragments (RPF; TPM+1) measured by Ribo-Seq; protein_lysate: protein abundance measured by SWATH-MS.",
                            "timepoint" = "Time point for which libtype variables were measured (TP1, TP2, TP3, and TP4).",
                            "mean" = "Mean values for libtype variables. RNA-Seq: n = 3; Ribo-Seq: n = 3; SWATH-MS: n ≥ 6.",
                            "se" = "Standard error of the mean (SEM) for mean values",
                            "product" = "Protein product given by Pfeiffer et al. (2019).",
                            "gene_symbol" = "Gene symbol according to COG 2020.",
                            "cog_id" = "ID according to COG 2020.",
                            "cog_name" = "Protein product according to COG 2020.",
                            "cog_category" = "Category according to COG 2020.",
                            "functional_pathway" = "Functional pathway according to COG 2020.",
                            "smap1Sense" = "Whether there is at least one SmAP1 binding site on the same strand for a given gene.",
                            # "smap1AntiSense" = "Whether there is at least one SmAP1 binding site on the opposite strand for a given gene.",
                            "asRNA" = "Whether there is at least one annotated antisense RNA (asRNA) according to de Almeida et al. (2019).",
                            "tps" = "Whether there is at least one transcript processing site (TPS) on the same strand for a given gene.",
                            "GC" = "GC content for a given gene.",
                            "GCdev" = "Difference between GC content for a given gene and the mean GC content considering all the genes.",
                            "HL" = "Experimentally determined half-life of a transcript according to Hundt et al. (2007).",
                            "cai" = "Codon adaptation index computed using as reference set the 5% most abundant proteins in this study.",
                            # structure(names = paste0("ChIPSeq_", paste0("Tfb", c("B", "D", "G"))), .Data = paste0("Whether there is at least one ", paste0("Tfb", c("B", "D", "G")), " binding site 150 nt upstream or downstream to the first codon of a gene. Binding sites were inferred using ChIP-Seq.")),
                            # structure(names = paste0("ChIPChip_", paste0("Tbp", c("B", "C", "E", "F"))), .Data = paste0("Whether there is at least one ", paste0("Tbp", c("B", "C", "E", "F")), " binding site 150 nt upstream or downstream to the first codon of a gene. Binding sites were inferred using ChIP-Chip.")),
                            # structure(names = paste0("ChIPChip_", paste0("Tfb", c("A", "B", "C", "D", "E", "F", "G"))), .Data = paste0("Whether there is at least one ", paste0("Tfb", c("A", "B", "C", "D", "E", "F", "G")), " binding site 150 nt upstream or downstream to the first codon of a gene. Binding sites were inferred using ChIP-Chip.")),
                            "Chromosome" = "The represented instance is located within the chromosome (NC_002607.1).",
                            "pNRC100" = "The representative instance is located within the plasmid pNRC100 (NC_001869.1).",
                            "pNRC200" = "The representative instance is located within the plasmid pNRC200 (NC_002608.1).",
                            # "IRs" = "Number of inverted repeats found within a given gene.",
                            "lfc2099" = "Log2-transformed fold change of RNase VNG2099 knockout vs. the control strain.",
                            "ISFamily" = "If the representative instance is located within an insertion sequence, this field will display the insertion sequence family."
                            # "protein_coding" = "Whether the representative instance encodes a protein."
                            )
  
atlasTidyDesc = tibble(column_name = names(atlasTidyColDescription),
                       description = atlasTidyColDescription)

nrtxDesc = tibble(column_name = c("representative",
                                  "product",
                                  "locus_tag"),
                  description = c("Locus tag for a given instance according to Pfeiffer et al. (2019). This locus tag may be a representative of many others if they were collapsed in our non-redundant transcriptome.",
                                  "Gene product given by Pfeiffer et al. (2019).",
                                  "Locus tags represented by the representative field. Might be a synonym or a locus tag for an almost identical gene that was collapsed by our non-redundant transcriptome approach."))

# creating a tab description object
atlasTabDesc = tibble(tab_name = c("non_redundant_tx_dictionary",
                                   "nrtx_column_description",
                                   "atlas_normalized_wide_data",
                                   "atlas_column_description_wide",
                                   "atlas_non_normalized_tidy_data",
                                   "atlas_column_description_tidy"),
                      description = c("The non-redundant transcriptome dictionary generated for this study. Source code and further info at: https://github.com/alanlorenzetti/halo_nr_tx .",
                                      "Column description for the non_redundant_tx_dictionary tab.",
                                      "Wide data table used to generate the interactive heat maps for this study.",
                                      "Column description for the atlas_normalized_wide_data tab.",
                                      "A tidy data table used to generate atlas_normalized_wide_data. This table contains the most essential (non-normalized) data for this study.",
                                      "Column description for the atlas_non_normalized_tidy_data tab."))

# also with nrtx table
atlasdata = list("tab_guide_readme" = atlasTabDesc,
                 "non_redundant_tx_dictionary" = nrtx,
                 "nrtx_column_description" = nrtxDesc,
                 "atlas_normalized_wide_data" = hmaFuncat %>% 
                   dplyr::rename(smap1Sense = "lsmSense",
                                 smap1AntiSense = "lsmAntiSense") %>% 
                   mutate(across(.cols = where(is.numeric),
                                 .fns = ~ str_replace(format(round(.x, 3), nsmall = 3), "^ *", "") %>% as.numeric())) %>% 
                   select(-starts_with("ChIP")) %>% 
                   select(-c(smap1AntiSense, protein_coding)) %>% 
                   mutate(gene_symbol = case_when(gene_symbol == "Undefined" ~ NA_character_,
                                                  TRUE ~ gene_symbol),
                          functional_pathway = case_when(functional_pathway == "Undefined" ~ NA_character_,
                                                  TRUE ~ functional_pathway)) %>% 
                   relocate(tps, .after = asRNA),
                 "atlas_column_description_wide" = atlasDesc,
                 "atlas_non_normalized_tidy_data" = abundLongFuncat %>% 
                   dplyr::rename(locus_tag_synonyms = "locus_tag.y",
                                 smap1Sense = "lsmSense",
                                 smap1AntiSense = "lsmAntiSense") %>% 
                   mutate(across(.cols = where(is.numeric),
                                 .fns = ~ str_replace(format(round(.x, 3), nsmall = 3), "^ *", "") %>% as.numeric())) %>% 
                   select(-starts_with("ChIP")) %>% 
                   select(-c(smap1AntiSense, protein_coding)) %>% 
                   mutate(gene_symbol = case_when(gene_symbol == "Undefined" ~ NA_character_,
                                                  TRUE ~ gene_symbol),
                          functional_pathway = case_when(functional_pathway == "Undefined" ~ NA_character_,
                                                         TRUE ~ functional_pathway)) %>% 
                   relocate(tps, .after = asRNA),
                 "atlas_column_description_tidy" = atlasTidyDesc
                 )

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
                  count = ptgs[[i]][[j]]$locus_tag %>% length(),
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
                             cluster == "BL41" ~ "Purple"),
           change_status = case_when(cluster == "Q1" ~ "Protein Up & mRNA Up",
                                     cluster == "Q2" ~ "Protein Up & mRNA Down",
                                     cluster == "Q3" ~ "Protein Down & mRNA Down",
                                     cluster == "Q4" ~ "Protein Down & mRNA Up",
                                     cluster == "BL12" ~ "Protein Up & mRNA Flat",
                                     cluster == "BL23" ~ "Protein Flat & mRNA Down",
                                     cluster == "BL34" ~ "Protein Down & mRNA Flat",
                                     cluster == "BL41" ~ "Protein Flat & mRNA Up")
    ) %>% 
    select(cluster,
           color,
           change_status,
           count,
           locus_tag)
}

# setting up a tab description tab
outfcclusters$tab_guide_readme = tibble(tab_name = c("TP2_vs_TP1",
                                                     "TP3_vs_TP2",
                                                     "TP4_vs_TP3",
                                                     "TP3_vs_TP1",
                                                     "TP4_vs_TP1",
                                                     "column_description"),
                                        description = c("Clusters identified in the transition from early exponential (TP1) to mid-exponential (TP2) growth phase.",
                                                        "Clusters identified in the transition from mid-exponential (TP2) to late exponential (TP3) growth phase.",
                                                        "Clusters identified in the transition from late exponential growth (TP3) to stationary (TP4) phase.",
                                                        "Clusters identified in the transition from early exponential (TP1) to late exponential (TP3) growth phase.",
                                                        "Clusters identified in the transition from early exponential (TP1) to stationary (TP4) phase.",
                                                        "Column description for the mentioned tabs."))

# setting up a column description tab
outfcclusters$column_description = tibble(column_name = c("cluster",
                                                          "color",
                                                          "change_status",
                                                          "count",
                                                          "locus_tag"),
                                          description = c("Cluster identifier according to the quadrant nomenclature convention using the Cartesian coordinate system. E.g., Q1: quadrant I; BL14: Border line between quadrant I and quadrant IV.",
                                                          "Color used in our figures to highlight the cluster.",
                                                          "Change status considering the transition from a physiological state to another. Take as an example the transition from TP1 to TP2: if the protein is upregulated and the mRNA is upregulate, the change status is Protein Up & mRNA Up. We use the term *Flat* to represent unchanging protein or mRNA levels.",
                                                          "Number of instances (locus tags) within a given cluster.",
                                                          "Locus tag for a given instance according to Pfeiffer et al. (2019). This locus tag may be a representative of many others if they were collapsed in our non-redundant transcriptome."))

# ordering according to panel of plots
outfcclusters = outfcclusters[c("tab_guide_readme",
                                "TP2_vs_TP1",
                                "TP3_vs_TP2",
                                "TP4_vs_TP3",
                                "TP3_vs_TP1",
                                "TP4_vs_TP1",
                                "column_description")]

# writing table containing multiple sheets
write.xlsx(outfcclusters,
           file = "results/supp_tables/Table_S_clusters.xlsx",
           overwrite = T)

# generating the table containing clusters extracted from abundance analysis ####
# creating a list of data frames to store locus tags
outabundclusters = list()

# getting data
for(i in names(ptgsAbund)){
  spdsheet = tibble()
  
  for(j in c("prot_bot_mrna_top", "prot_non_mrna_top", "prot_bot_prot_non_mrna_top")){
    line = tibble(cluster = j,
                  count = ptgsAbund[[i]][[j]] %>% length(),
                  locus_tag = paste0(ptgsAbund[[i]][[j]], collapse = ","))
    spdsheet = bind_rows(spdsheet, line)
  }
  
  outabundclusters[[i]] = spdsheet %>% 
    mutate(status = case_when(cluster == "prot_bot_mrna_top" ~ "Protein low level & mRNA high level",
                                     cluster == "prot_non_mrna_top" ~ "Protein undetected & mRNA high level",
                                     cluster == "prot_bot_prot_non_mrna_top" ~ "Protein low level or undetected & mRNA high level")
    ) %>% 
    select(cluster,
           status,
           count,
           locus_tag)
}

# setting up a tab description tab
outabundclusters$tab_guide_readme = tibble(tab_name = c("TP1",
                                                        "TP2",
                                                        "TP3",
                                                        "TP4",
                                                        "union",
                                                        "column_description"),
                                           description = c("Clusters identified using protein and mRNA data obtained from early exponential (TP1) growth phase.",
                                                           "Clusters identified using protein and mRNA data obtained from mid-exponential (TP2) growth phase.",
                                                           "Clusters identified using protein and mRNA data obtained from late exponential (TP3) growth phase.",
                                                           "Clusters identified using protein and mRNA data obtained from stationary phase (TP4).",
                                                           "Union of locus tags identified across the growth curve (TP1, TP2, TP3, and TP4).",
                                                           "Column description for the mentioned tabs."))

# setting up a column description tab
outabundclusters$column_description = tibble(column_name = c("cluster",
                                                             "status",
                                                             "count",
                                                             "locus_tag"),
                                             description = c("Cluster identifier according to the expression status of protein and mRNA. prot_bot: 20% lowest abundance proteins (bottom); prot_non: Undetected proteins; mrna_top: 20% greatest abundance mRNAs (top).",
                                                             "A more intuitive description of the cluster field.",
                                                             "Number of instances (locus tags) within a given cluster.",
                                                             "Locus tag for a given instance according to Pfeiffer et al. (2019). This locus tag may be a representative of many others if they were collapsed in our non-redundant transcriptome."))

# ordering according to panel of plots
outabundclusters = outabundclusters[c("tab_guide_readme",
                                      "TP1",
                                      "TP2",
                                      "TP3",
                                      "TP4",
                                      "union",
                                      "column_description")]

# writing table containing multiple sheets
write.xlsx(outabundclusters,
           file = "results/supp_tables/Table_S_clustersAbund.xlsx",
           overwrite = T)

# # creating a table of potentially post-transcriptionally repressed genes ####
# # with enrichment p-values
# # abundance-based
# timepoints = names(ptgsAbund)
# modus = c("prot_bot_mrna_top", "prot_non_mrna_top", "prot_bot_prot_non_mrna_top")
# testcols = c("lsmSense", "asRNA", "tps")
# 
# ptgsAbundRepTib = tibble()
# for(i in timepoints){
#   for(j in modus){
#     ptgsAbundRepTib = bind_rows(ptgsAbundRepTib,
#                                 tibble(criteria = j,
#                                        timepoint = i,
#                                        count = ptgsAbund[[i]][[j]] %>% length(),
#                                        smap1_enrich_p = enrichAnalysis(hmaFuncat, ptgsAbund[[i]][[j]], "lsmSense", "yes"),
#                                        asrna_enrich_p = enrichAnalysis(hmaFuncat, ptgsAbund[[i]][[j]], "asRNA", "yes"),
#                                        tps_enrich_p = enrichAnalysis(hmaFuncat, ptgsAbund[[i]][[j]], "tps", "yes"),
#                                        locus_tag = ptgsAbund[[i]][[j]] %>% paste0(collapse= ",")))
#   }
# }
# 
# ptgsAbundRepTib = ptgsAbundRepTib %>% 
#   arrange(criteria)
# 
# # fc-based
# timepoints = names(ptgs)
# modus = names(ptgs$union$Q4)
# 
# ptgsFcRepTib = tibble()
# for(i in timepoints){
#     ptgsFcRepTib = bind_rows(ptgsFcRepTib,
#                              tibble(criteria = "prot_downregulated_mrna_upregulated",
#                                     timepoint = i,
#                                     count = ptgs[[i]][["Q4"]][["locus_tag"]] %>% length(),
#                                     smap1_enrich_p = enrichAnalysis(hmaFuncat, ptgs[[i]][["Q4"]][["locus_tag"]], "lsmSense", "yes"),
#                                     asrna_enrich_p = enrichAnalysis(hmaFuncat, ptgs[[i]][["Q4"]][["locus_tag"]], "asRNA", "yes"),
#                                     tps_enrich_p = enrichAnalysis(hmaFuncat, ptgs[[i]][["Q4"]][["locus_tag"]], "tps", "yes"),
#                                     locus_tag = ptgs[[i]][["Q4"]][["locus_tag"]] %>% paste0(collapse= ",")))
# }

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
  mutate(functional_class = case_when(functional_class == "all" ~ "All",
                                      TRUE ~ functional_class)) %>% 
  mutate(varlvl = paste0(tested_var, ":", tested_level)) %>% 
  mutate(pvalue = paste0(tested_var, ":", format(pvalue, scientific = T, digits = 3))) %>% 
  group_by(criteria,timepoint,count,functional_class,locus_tag) %>%
  summarise(mechanism = paste0(varlvl, collapse = ","),
            pvalues = paste0(pvalue, collapse = ",")) %>% 
  select(cluster = criteria,
         timepoint,
         cog_category = functional_class,
         significant_features = mechanism,
         pvalues,
         count,
         locus_tag)

# performing the same operations for fold change
# based analysis

# abundance approach
# categorical variables
fcStats = tibble()
criteria = c("Q1", "Q2", "Q3", "Q4", "BL12", "BL23", "BL34", "BL41")
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
  mutate(functional_class = case_when(functional_class == "all" ~ "All",
                                      TRUE ~ functional_class)) %>% 
  mutate(varlvl = paste0(tested_var, ":", tested_level)) %>% 
  mutate(pvalue = paste0(tested_var, ":", format(pvalue, scientific = T, digits = 3))) %>% 
  group_by(criteria,timepoint,count,functional_class,locus_tag) %>%
  summarise(mechanism = paste0(varlvl, collapse = ","),
            pvalues = paste0(pvalue, collapse = ",")) %>% 
  select(cluster = criteria,
         timepoint,
         cog_category = functional_class,
         significant_features = mechanism,
         pvalues,
         count,
         locus_tag)

# creating a list to store tables
ptgsEnrich = list("abs_abundance_based_analysis" = ptgsAbundRepTib,
                  "rel_abundance_based_analysis" = ptgsFcRepTib)

# creating a tab guide
ptgsEnrich$tab_guide_readme = tibble(tab_name = c("abs_abundance_based_analysis",
                                                  "rel_abundance_based_analysis",
                                                  "column_description"),
                                     description = c("Feature enrichment analysis for clusters of potentially post-transcriptionally regulated genes defined by the absolute abundance-based approach.",
                                                     "Feature enrichment analysis for clusters of potentially post-transcriptionally regulated genes defined by the relative abundance-based approach.",
                                                     "Column description for the mentioned tabs."))

# setting up a column description tab
ptgsEnrich$column_description = tibble(column_name = c("cluster",
                                                       "timepoint",
                                                       "cog_category",
                                                       "significant_features",
                                                       "pvalues",
                                                       "count",
                                                       "locus_tag"),
                                       description = c("Cluster identifier according to the absolute abundance- or relative abundance-based analysis. Please, refer to the supplemental files for each one of the approaches for more information.",
                                                       "Time point for the absolute abundance-based approach and time point contrast for the relative abundance-based approach.",
                                                       "Category according to COG 2020. 'All' is an arbitrary category including all the categories.",
                                                       "Significantly enriched/different features and levels (feature:level format) that could help us understand why protein and mRNA levels are incoherent for a given cluster. E.g., 'SmAP1:yes' indicates that the cluster is enriched for SmAP1 binding. E.g., 'CAI:higher:0.911:0.773' indicates that the CAI of the cluster (0.911) is significantly higher than the CAI of the backgroun genes (0.773). Different features are separated by commas.",
                                                       "P-values for the significant features and levels. E.g., 'SmAP1:6.41e-03' means that SmAP1 binding is enriched in a given cluster with a p-value of 6.41e-03. Different p-values are separated by commas.",
                                                       "Number of instances (locus tags) within a given cluster.",
                                                       "Locus tag for a given instance according to Pfeiffer et al. (2019). This locus tag may be a representative of many others if they were collapsed in our non-redundant transcriptome."))

# ordering according to panel of plots
ptgsEnrich = ptgsEnrich[c("tab_guide_readme",
                          "abs_abundance_based_analysis",
                          "rel_abundance_based_analysis",
                          "column_description")]
  
# writing table with multi sheets
write.xlsx(ptgsEnrich,
           file = "results/supp_tables/Table_S_ptgsEnrich.xlsx",
           overwrite = T)

# mobilization events supp ####
mobSupp = list()

mobSupp$mobilization_events = dfinsdel %>% 
  mutate(strain = case_when(strain == "&Delta;*ura3* A" ~ "dura3_A",
                            strain == "&Delta;*ura3* B" ~ "dura3_B",
                            strain == "&Delta;*ura3* C" ~ "dura3_C",
                            strain == "&Delta;*ura3* &Delta;*smap1* A" ~ "dura3dsmap1_A",
                            strain == "&Delta;*ura3* &Delta;*smap1* B" ~ "dura3dsmap1_B",
                            strain == "&Delta;*ura3* &Delta;*smap1* C" ~ "dura3dsmap1_C",
                            TRUE ~ strain),
         ISFamily = case_when(ISFamily == "Outras famílias" ~ "Other families",
                              TRUE ~ ISFamily),
         cluster = case_when(svType == "Insertion" ~ str_replace(cluster, ":", ":I:"),
                             svType == "Excision" ~ str_replace(cluster, ":", ":E:"),
                             TRUE ~ cluster),
         ISName = str_replace_all(ISName, "\\*", ""),
         ISFamily = str_replace_all(ISFamily, "\\*", ""))

mobSupp$column_description = tibble(column_name = colnames(mobSupp$mobilization_events),
                                   description = c("strain" = "Halobacterium salinarum NRC-1 strain and replicate (A, B, and C). dura3 is the parent Ura3 knockout strain (control). dura3dsmap1 is the SmAP1 knockout strain.",
                                                   "cluster" = "Structural variant cluster identifier. Unique within sample. I: insertion cluster. E: excision cluster.",
                                                   "replicon" = "Replicon where the structural variant cluster was found. NC_002607.1: Chromosome; NC_001869.1: Plasmid pNRC100; NC_002608.1: Plasmid pNRC200.",
                                                   "ISName" = "Name of the insertion sequence matching the the structural variant cluster.",
                                                   "ISFamily" = "Insertion sequence family for a given insertion sequence.",
                                                   "meanStart" = "Position of the detected cluster. Mean is applied if a cluster is supported by many mobilization events. The position is relative to the coordinates of a modified reference genome that excludes long duplications. See methods for details.",
                                                   "sdStart" = "Uncertainty of position of the detected cluster. Standad deviation is applied if a cluster is supported by many mobilization events. The position is relative to the coordinates of a modified reference genome that excludes long duplications. See methods for details.",
                                                   "meanLength" = "Length of the detected cluster. Mean is applied if a cluster is supported by many mobilization events.",
                                                   "sdLength" = "Length of the detected cluster. Standard deviation is applied if a cluster is supported by many mobilization events.",
                                                   "count" = "Number of reads supporting the mobilization cluster, that is, number of mobilization events observed within a cluster.",
                                                   "rnames" = "Name of the reads supporting the mobilization events for a given cluster.",
                                                   "status" = "Occurency status of a mobilization cluster inferred by *count* divided by the local sequencing depth (events by coverage; e/c). e/c ≤ 0.1: Rare; 0.1 < e/c ≤ 0.5: Common; e/c > 0.5: Predominant.",
                                                   "svType" = "Type of structural variant cluster: Insertion or Excision."))

# writing table with multi sheets
write.xlsx(mobSupp,
           file = "results/supp_tables/Table_S_mob_events.xlsx",
           overwrite = T)
  
# membrane proteins that are likely undetected ####
# checking the proteins
undetectedRemovDF = nrtx %>%
  filter(representative %in%
           ptgsAbund$union$prot_non_mrna_top[ptgsAbund$union$prot_non_mrna_top %in%
                                               c(membProtsFinal, topsconMemb$two)])

# updating our final dataframe
# removing proteins 
# that are likely transmembrane
undetectedRemovDF = undetectedRemovDF %>% 
  filter(!representative %in% membFP)

# due to their transmembrane nature
transmembSupp = list()

transmembSupp$transmembrane_proteins = undetectedRemovDF %>% 
  left_join(x = .,
            y = topsconRes,
            by = c("representative" = "locus_tag")) %>% 
  select(-c(Length, SignalPeptide)) %>% 
  mutate(experimental_evidence = case_when(representative %in% membProtsFinal ~ "yes",
                                           TRUE ~ "no")) %>% 
  relocate(locus_tag,
           .after = experimental_evidence)
  
transmembSupp$column_description =  tibble(column_name = c("representative",
                                                           "product",
                                                           "numTM",
                                                           "experimental_evidence",
                                                           "locus_tag"),
                                           description = c("Locus tag for a given instance according to Pfeiffer et al. (2019). This locus tag may be a representative of many others if they were collapsed in our non-redundant transcriptome.",
                                                           "Gene product given by Pfeiffer et al. (2019).",
                                                           "Number of transmembrane domains predicted by TOPCONS webserver.",
                                                           "Whether there is proteome experimental evidence for membrane samples in Goo et al. (2003) or Klein et al. (2005).",
                                                           "Locus tags represented by the representative field. Might be a synonym or a locus tag for an almost identical gene that was collapsed by our non-redundant transcriptome approach."))

write.xlsx(transmembSupp,
           file = "results/supp_tables/Table_S_transmemb.xlsx",
           overwrite = T)

# saving 2099 differential expression analysis ####
DE2099 = list()

DE2099$DE_analysis_2099 = res2099 %>% 
  mutate(differentially_expressed = case_when(abs(logFC) >= log2fcthreshold & P.Value < padjthreshold ~ "yes",
                                              TRUE ~ "no"),
         direction_of_change = case_when(logFC >= log2fcthreshold & differentially_expressed == "yes" ~ "upregulated",
                                         logFC <= log2fcthreshold & differentially_expressed == "yes" ~ "downregulated",
                                         TRUE ~ "non-differentially_expressed")) %>% 
  arrange(desc(differentially_expressed), logFC) %>% 
  mutate(across(.cols = where(is.numeric),
                .fns = ~ str_replace(format(round(.x, 3), nsmall = 3), "^ *", "") %>% as.numeric()))

DE2099$column_description =  tibble(column_name = c("representative",
                                                    "adj.P.Val",
                                                    "P.Value",
                                                    "t",
                                                    "B",
                                                    "logFC",
                                                    "differentially_expressed",
                                                    "direction_of_change"),
                                    description = c("Locus tag for a given instance according to Pfeiffer et al. (2019). This locus tag may be a representative of many others if they were collapsed in our non-redundant transcriptome.",
                                                    "limma::topTable output field: adjusted p-value or q-value (Benjamini-Hochberg procedure). ",
                                                    "limma::topTable output field: raw p-value.",
                                                    "limma::topTable output field: moderated t-statistic.",
                                                    "limma::topTable output field: log-odds that the gene is differentially expressed.",
                                                    "limma::topTable output field: estimate of the log2-fold-change corresponding to the contrast (mutant/control).",
                                                    "Whether the gene is differentially expressed according to our criteria.",
                                                    "Direction of change for differentially expressed genes (upregulated or downregulated)."))

write.xlsx(DE2099,
           file = "results/supp_tables/Table_S_2099.xlsx",
           overwrite = T)

# differential expression analysis for growth curve ####

# this function will convert protein ids
# of the original de analysis to the
# ones we anchoring in this study
prepareLocusTag = function(dictobj = dictobj,
                           inputDF = inputDF) {
  
  originalLocusTags = inputDF$name %>% unique()
  
  myNewDF = left_join(x = dictobj,
                      y = inputDF,
                      by = c("locus_tag" = "name")) %>% 
    dplyr::filter(locus_tag %in% originalLocusTags) %>% 
    dplyr::select(-c(product,locus_tag))
  
  return(myNewDF)
}

# this function will create new variables
# for DE dataset
createNewDEvars = function(inputDF = inputDF) {
  
  if("adj_pval" %in% colnames(inputDF)){
    rna_or_protein = "protein"
  }else if("padj" %in% colnames(inputDF)){
    rna_or_protein = "mRNA"
  }
  
  if(rna_or_protein == "mRNA"){
    myNewDF = inputDF %>% 
      mutate(differentially_expressed = case_when(abs(log2FoldChange) >= log2fcthreshold & padj < padjthreshold ~ "yes",
                                                  TRUE ~ "no"),
             direction_of_change = case_when(log2FoldChange >= log2fcthreshold & differentially_expressed == "yes" ~ "upregulated",
                                             log2FoldChange <= log2fcthreshold & differentially_expressed == "yes" ~ "downregulated",
                                             TRUE ~ "non-differentially_expressed")) %>% 
      mutate(across(.cols = where(is.numeric),
                    .fns = ~ str_replace(format(round(.x, 3), nsmall = 3), "^ *", "") %>% as.numeric())) %>% 
      dplyr::rename("representative" = target_id)
    
  }else if(rna_or_protein == "protein"){
    myNewDF = inputDF %>% 
      mutate(differentially_expressed = case_when(abs(diff) >= log2fcthreshold & adj_pval < padjthreshold ~ "yes",
                                                  TRUE ~ "no"),
             direction_of_change = case_when(diff >= log2fcthreshold & differentially_expressed == "yes" ~ "upregulated",
                                             diff <= log2fcthreshold & differentially_expressed == "yes" ~ "downregulated",
                                             TRUE ~ "non-differentially_expressed")) %>% 
      mutate(across(.cols = where(is.numeric),
                    .fns = ~ str_replace(format(round(.x, 3), nsmall = 3), "^ *", "") %>% as.numeric()))
  }
  
  return(myNewDF)
}

DEgrowth = list()

DEgrowth$tab_guide_readme = tibble(tab_name = c("mrna_TP2_vs_TP1",
                                                "mrna_TP3_vs_TP2",
                                                "mrna_TP4_vs_TP3",
                                                "mrna_TP3_vs_TP1",
                                                "mrna_TP4_vs_TP1",
                                                "mrna_column_description",
                                                "protein_TP2_vs_TP1",
                                                "protein_TP3_vs_TP2",
                                                "protein_TP4_vs_TP3",
                                                "protein_TP3_vs_TP1",
                                                "protein_TP4_vs_TP1",
                                                "protein_column_description"),
                                   description = c("mRNA differential expression analysis. Contrast TP2 vs. TP1. DESeq2::results output.",
                                                   "mRNA differential expression analysis. Contrast TP3 vs. TP2. DESeq2::results output.",
                                                   "mRNA differential expression analysis. Contrast TP4 vs. TP3. DESeq2::results output.",
                                                   "mRNA differential expression analysis. Contrast TP3 vs. TP1. DESeq2::results output.",
                                                   "mRNA differential expression analysis. Contrast TP4 vs. TP1. DESeq2::results output.",
                                                   "Column description for mRNA differential expression analysis tabs.",
                                                   "Protein differential expression analysis. Contrast TP2 vs. TP1. proDA::test_diff output.",
                                                   "Protein differential expression analysis. Contrast TP3 vs. TP2. proDA::test_diff output.",
                                                   "Protein differential expression analysis. Contrast TP4 vs. TP3. proDA::test_diff output.",
                                                   "Protein differential expression analysis. Contrast TP3 vs. TP1. proDA::test_diff output.",
                                                   "Protein differential expression analysis. Contrast TP4 vs. TP1. proDA::test_diff output.",
                                                   "Column description for Protein differential expression analysis tabs."))

DEgrowth$mrna_TP2_vs_TP1 = results$totrna_TP2_vs_TP1 %>% createNewDEvars()
DEgrowth$mrna_TP3_vs_TP2 = results$totrna_TP3_vs_TP2 %>% createNewDEvars()
DEgrowth$mrna_TP4_vs_TP3 = results$totrna_TP4_vs_TP3 %>% createNewDEvars()
DEgrowth$mrna_TP3_vs_TP1 = results$totrna_TP3_vs_TP1 %>% createNewDEvars()
DEgrowth$mrna_TP4_vs_TP1 = results$totrna_TP4_vs_TP1 %>% createNewDEvars()

DEgrowth$mrna_column_description = tibble(column_name = c("representative",
                                                          "baseMean",
                                                          "log2FoldChange",
                                                          "lfcSE",
                                                          "stat",
                                                          "pvalue",
                                                          "padj",
                                                          "differentially_expressed",
                                                          "direction_of_change"),
                                          description = c("Locus tag for a given instance according to Pfeiffer et al. (2019). This locus tag may be a representative of many others if they were collapsed in our non-redundant transcriptome.",
                                                          "DESeq2::results output field: mean of normalized counts for all samples.",
                                                          "DESeq2::results output field: log2 fold change (MLE).",
                                                          "DESeq2::results output field: standard error.",
                                                          "DESeq2::results output field: Wald statistic.",
                                                          "DESeq2::results output field: Wald test p-value.",
                                                          "DESeq2::results output field: BH adjusted p-values.",
                                                          "Whether the gene is differentially expressed according to our criteria.",
                                                          "Direction of change for differentially expressed genes (upregulated or downregulated)."))

DEgrowth$protein_TP2_vs_TP1 = resultsDEProt$protein_TP2_vs_TP1$results %>% prepareLocusTag(dictobj = nrtxsep, inputDF = .) %>% createNewDEvars()
DEgrowth$protein_TP3_vs_TP2 = resultsDEProt$protein_TP3_vs_TP2$results %>% prepareLocusTag(dictobj = nrtxsep, inputDF = .) %>% createNewDEvars()
DEgrowth$protein_TP4_vs_TP3 = resultsDEProt$protein_TP4_vs_TP3$results %>% prepareLocusTag(dictobj = nrtxsep, inputDF = .) %>% createNewDEvars()
DEgrowth$protein_TP3_vs_TP1 = resultsDEProt$protein_TP3_vs_TP1$results %>% prepareLocusTag(dictobj = nrtxsep, inputDF = .) %>% createNewDEvars()
DEgrowth$protein_TP4_vs_TP1 = resultsDEProt$protein_TP4_vs_TP1$results %>% prepareLocusTag(dictobj = nrtxsep, inputDF = .) %>% createNewDEvars()

DEgrowth$protein_column_description = tibble(column_name = c("representative",
                                                             "pval",
                                                             "adj_pval",
                                                             "diff",
                                                             "t_statistic",
                                                             "se",
                                                             "df",
                                                             "avg_abundance",
                                                             "n_approx",
                                                             "n_obs",
                                                             "differentially_expressed",
                                                             "direction_of_change"),
                                             description = c("Locus tag for a given instance according to Pfeiffer et al. (2019). This locus tag may be a representative of many others if they were collapsed in our non-redundant transcriptome.",
                                                             "proDA::test_diff output field: the p-value of the statistical test.",
                                                             "proDA::test_diff output field: the multiple testing adjusted p-value.",
                                                             "proDA::test_diff output field: the difference that particular coefficient makes. In differential expression analysis this value is also called log fold change, which is equivalent to the difference on the log scale.",
                                                             "proDA::test_diff output field: the diff divided by the standard error se",
                                                             "proDA::test_diff output field: the standard error associated with the diff",
                                                             "proDA::test_diff output field: the degrees of freedom, which describe the amount of available information for estimating the se. They are the sum of the number of samples the protein was observed in, the amount of information contained in the missing values, and the degrees of freedom of the variance prior.",
                                                             "proDA::test_diff output field: the estimate of the average abundance of the protein across all samples.",
                                                             "proDA::test_diff output field: the approximated information available for estimating the protein features, expressed as multiple of the information contained in one observed value.",
                                                             "proDA::test_diff output field: the number of samples a protein was observed in.",
                                                             "Whether the gene is differentially expressed according to our criteria.",
                                                             "Direction of change for differentially expressed genes (upregulated or downregulated)."))

write.xlsx(DEgrowth,
           file = "results/supp_tables/Table_S_DEgrowth.xlsx",
           overwrite = T)

# copying figures ####
# made for this paper to the 
# appropriate directory
figFiles = c(
"plots/17_ptgsFeatures_venn.png" = "M_Genes potentially subject to post-transcriptional regulation",
"plots/abundanceHeatmap_en.png" = "M_An atlas of transcriptome, ribosome profile, and proteome for Halobacterium salinarum NRC-1",
"plots/18_abund_fc_prot_mrna_panel.png" = "M_Genes following patterns compatible with post-transcriptional regulation",
"plots/gvp_traj.png" = "dependencyFig_Post-transcriptional_regulation_of_gvp1_operons",
"plots/gvp_traj.svg" = "dependencyFig_Post-transcriptional_regulation_of_gvp1_operons",
"plots/mobileElPanelFeatures_en.png" = "M_Protein and mRNA levels of mobile elements",
"plots/mobilizationComparisonPerFamily_en.png" = "M_Detected mobilizations for decomposed insertion sequence families",
"plots/13_prot_vs_gc_tpwise.png" = "S_Protein levels are associated with transcript GC content",
"plots/15_prot_botOrNon_mrna_top_venn.png" = "S_Venn diagrams of potentially post-transcriptionally regulated genes shared among different physiological states",
"plots/abundanceHeatmap_ptgs/TP2_vs_TP1_Q4.png" = "S_Atlas section of potentially post-transcriptionally regulated genes in the transition from TP1 to TP2",
"plots/upset_plot_fcbased.png" = "S_UpSet plot of potentially post-transcriptionally regulated genes shared among different physiological states transitions",
"plots/gvp_heat.png" = "S_Protein-mRNA dynamics and various features of genes encoding gas vesicle biogenesis proteins",
"plots/smap1_panel.png" = "S_SmAP1 descriptive features",
"plots/mob_observed_events_panel.png" = "S_Detected insertion and excision events",
"plots/smap1_ura3_growth_curve.png" = "S_Growth_curve_of_ura3_and_ura3smap1_strains_v2",
"plots/abundanceHeatmap_expanded_en.pdf" = "S_Atlas heat map expanded")

extensions = names(figFiles) %>% str_replace(".*(\\..*)$", "\\1")
orifile = names(figFiles)
destfile = paste0("results/figures/", figFiles %>% str_replace_all(" ", "_"), extensions)

file.copy(from = orifile, to = destfile, overwrite = T)
file.copy(from = "plots/abundanceHeatmap_expanded_en.pdf", to = "results/supp_tables/Figure_S_heatmap_expanded.pdf", overwrite = T)

# sessionInfo object
seshInfo = list(sys_time = Sys.time(),
                session_info = sessionInfo())

save(seshInfo,
     file = "scripts/session_info_latest_run.RData")
