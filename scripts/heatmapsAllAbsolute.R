# alorenzetti 20201104

# description ####
# this script will plot a heatmap
# based on absolute counts
# for proteins, mRNAs and RPFs
# functional categorization
# here I am trying to integrate
# the time series into a single heatmap

# preparing main abundance dataset ####
# dataset will be normalized within each timepoint
hma = abund %>% 
  dplyr::select(locus_tag,
                starts_with("mean_abundance_protein_lysate"),
                starts_with("mean_abundance_rna_total"),
                starts_with("mean_abundance_rna_ribofraction")) %>% 
  dplyr::select(locus_tag,
                ends_with("TP1"),
                ends_with("TP2"),
                ends_with("TP3"),
                ends_with("TP4")) %>% 
  drop_na()

# normalizing TP1
hmaM1 = hma %>% 
  dplyr::select(ends_with("TP1")) %>% 
  as.matrix() %>% 
  normalize.quantiles() %>% 
  as_tibble()

# normalizing TP2
hmaM2 = hma %>% 
  dplyr::select(ends_with("TP2")) %>% 
  as.matrix() %>% 
  normalize.quantiles() %>% 
  as_tibble()

# normalizing TP3
hmaM3 = hma %>% 
  dplyr::select(ends_with("TP3")) %>% 
  as.matrix() %>% 
  normalize.quantiles() %>% 
  as_tibble()

# normalizing TP4
hmaM4 = hma %>% 
  dplyr::select(ends_with("TP4")) %>% 
  as.matrix() %>% 
  normalize.quantiles() %>% 
  as_tibble()

# unifying matrices
hmaM = base::cbind(hmaM1, hmaM2, hmaM3, hmaM4)
colnames(hmaM) = colnames(hma)[-1]
rownames(hmaM) = hma$locus_tag

# reordering cols
cols = paste0(c(rep("mean_abundance_protein_lysate", 4),
                rep("mean_abundance_rna_total", 4),
                rep("mean_abundance_rna_ribofraction", 4)),
              rep(paste0("_TP", 1:4), 3))

hmaM = hmaM[,cols] %>%
  as.matrix()

# creating object with functional categories
hmaFuncat = left_join(hma, dictFunCat,
                      by = c("locus_tag" = "pfeiLocusTag")) %>% 
  dplyr::select(-locus_tag.y)

# preparing datasets for additional matrices ####
# translational efficiency and
# ribosome occupancy
hma = hmaM %>% 
  as_tibble()

hma$locus_tag = rownames(hmaM)
hma = hma %>% 
  relocate(locus_tag)
hma = hma %>% 
  mutate(TE_TP1 = mean_abundance_protein_lysate_TP1 / mean_abundance_rna_total_TP1,
         TE_TP2 = mean_abundance_protein_lysate_TP2 / mean_abundance_rna_total_TP2,
         TE_TP3 = mean_abundance_protein_lysate_TP3 / mean_abundance_rna_total_TP3,
         TE_TP4 = mean_abundance_protein_lysate_TP4 / mean_abundance_rna_total_TP4) %>% 
  mutate(RO_TP1 = mean_abundance_rna_ribofraction_TP1 / mean_abundance_rna_total_TP1,
         RO_TP2 = mean_abundance_rna_ribofraction_TP2 / mean_abundance_rna_total_TP2,
         RO_TP3 = mean_abundance_rna_ribofraction_TP3 / mean_abundance_rna_total_TP3,
         RO_TP4 = mean_abundance_rna_ribofraction_TP4 / mean_abundance_rna_total_TP4)

# translational efficiency
hmaTE = hma %>% 
  select(starts_with("TE")) %>% 
  as.matrix()
rownames(hmaTE) = rownames(hmaM)

# ribosome occupancy
hmaRO = hma %>% 
  select(starts_with("RO")) %>% 
  as.matrix()
rownames(hmaRO) = rownames(hmaM)

# defining colors and color functions ####
# 21 manual colors;
# mostly extracted from
# ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10` and
# ggthemes_data$tableau$`color-palettes`$regular$`Tableau 20`
arCOGcols = c(
  "Amino acid transport and metabolism" = "#4E79A7", # blue
  "Carbohydrate transport and metabolism" = "#A0CBE8", # light blue
  "Cell cycle control, cell division, chromosome partitioning" = "#F28E2B", # orange
  "Cell motility" = "#FFBE7D", # light orange
  "Cell wall/membrane/envelope biogenesis" = "#59A14F", # green
  "Chromatin structure and dynamics" = "yellow", # 
  "Coenzyme transport and metabolism" = "#8CD17D", # light green
  "Defense mechanisms" = "#B6992D", # yellow green
  "Energy production and conversion" = "#F1CE63", # yellow
  "Extracellular structures" = "blue", # 
  "General function prediction only" = "grey70", 
  "Inorganic ion transport and metabolism" = "#86BCB6", # light teal
  "Function unknown" = "#79706E", # dark grey
  "Intracellular trafficking, secretion, and vesicular transport" = "#E15759", # red
  "Lipid transport and metabolism" = "#FF9D9A", # pink
  "Mobilome: prophages, transposons" = "#D37295", # pink
  "Nucleotide transport and metabolism" = "orchid1", # orchid1
  "Posttranslational modification, protein turnover, chaperones" = "darkturquoise", # darkturquoise
  "Replication, recombination and repair" = "skyblue2", # skyblue2,
  "RNA processing and modification" = "#B07AA1", # purple
  "Secondary metabolites biosynthesis, transport and catabolism" = "#9D7660", # brown
  "Signal transduction mechanisms" = "#D7B5A6", # light orange
  "Transcription" = "#499894", # teal
  "Translation, ribosomal structure and biogenesis" = "maroon" # maroon
)

# getting classes included in arCOG
# this will be used when setting the legend for arCOG
arCOGClasses = hmaFuncat$cog_category %>% sort() %>% unique()
arCOGcols = arCOGcols[names(arCOGcols) %in% arCOGClasses]

# defining colors for the heatmap and annots
heatCols = list(
  lsmCol = c("no" = "white",
             "yes" = "#E15759"),
  arCOGCol = arCOGcols,
  utrSizeCol = c("Unknown" = "grey70",
                 "leaderless" = "#4E79A7",
                 "short" = "#EDC948",
                 "mid" = "#F28E2B",
                 "long" = "#E15759"),
  asRNACol = c("no" = "white",
               "yes" = "#B07AA1"),
  mfeCol = colorRamp2(breaks = c(0, hmaFuncat$utrmfe %>% min(na.rm = T) %>% floor()),
                      colors = c("white", "#4E79A7")),
  HLCol = colorRamp2(breaks = c(0, hmaFuncat$HL %>% max(na.rm = T) %>% ceiling()),
                     colors = c("white", "#4E79A7")),
  caiCol = colorRamp2(breaks = c(0.5, hmaFuncat$cai %>% max(na.rm = T)),
                      colors = c("white", "#4E79A7")),
  GCdevcol = colorRamp2(breaks = c(-0.1, 0, 0.1),
                        colors = c("#4E79A7", "white", "#E15759")),
  TFCol = c("no" = "white",
            "yes" = "black")
)

# defining colors and values for legends ####
heatLegs = list(
  arCOG = Legend(title = "arCOG",
                 at = names(heatCols$arCOGCol),
                 legend_gp  = gpar(fill = heatCols$arCOGCol %>% unname()),
                 border="black"),
  lsmSense = Legend(title = "LSm at Mid. Exponential",
                    at = c("no", "yes"),
                    labels = c("No", "Yes"),
                    legend_gp  = gpar(fill = heatCols$lsmCol),
                    border="black"),
  asRNA = Legend(title = "asRNA",
                 at = c("no", "yes"),
                 labels = c("No", "Yes"),
                 legend_gp  = gpar(fill = heatCols$asRNACol),
                 border="black"),
  utrSize = Legend(title = "5' UTR Length",
                   at = c("Unknown", "leaderless", "short", "mid", "long"),
                   labels = c("Unknown", "Absent", "Short (1-10 nt)", "Mid (11-100 nt)", "Long (101-250 nt)"),
                   legend_gp  = gpar(fill = heatCols$utrSizeCol),
                   border="black"),
  mfe = Legend(title = "5' UTR MFE",
               col_fun = heatCols$mfeCol,
               border="black"),
  HL = Legend(title = "Half life (min)",
              col_fun = heatCols$HLCol,
              border="black"),
  cai = Legend(title = "CAI",
               col_fun = heatCols$caiCol,
               border="black"),
  GCdev = Legend(title = "GC deviation",
                 col_fun = heatCols$GCdevcol,
                 border = "black"),
  TFs = Legend(title = "Transcription Factors",
               at = c("no", "yes"),
               labels = c("No", "Yes"),
               legend_gp = gpar(fill = heatCols$TFCol),
               border = "black")
)

# defining annotation columns ####
row_ha = HeatmapAnnotation(which = "row",
                           arCOG = anno_simple(hmaFuncat$arCOG,
                                               border = T,
                                               col = heatCols$arCOGCol),
                           lsmSense = anno_simple(hmaFuncat$lsmSense,
                                                  col = heatCols$lsmCol,
                                                  border = T),
                           lsmAntiSense = anno_simple(hmaFuncat$lsmAntiSense,
                                                      col = heatCols$lsmCol,
                                                      border = T),
                           asRNA = anno_simple(hmaFuncat$asRNA,
                                               col = heatCols$asRNACol,
                                               border = T),
                           utrSize = anno_simple(hmaFuncat$utrSize,
                                                 col = heatCols$utrSizeCol,
                                                 border = T),
                           mfe = anno_simple(hmaFuncat$utrmfe,
                                             col = heatCols$mfeCol,
                                             border = T),
                           HL = anno_simple(hmaFuncat$HL,
                                            col = heatCols$HLCol,
                                            border = T),
                           cai = anno_simple(hmaFuncat$cai,
                                             col = heatCols$caiCol,
                                             border = T),
                           GCdev = anno_simple(hmaFuncat$GCdev,
                                               col = heatCols$GCdevcol,
                                               border = T),
                           cstfbB = anno_simple(hmaFuncat$ChIPSeq_tfbB,
                                                col = heatCols$TFCol,
                                                border = T),
                           cstfbD = anno_simple(hmaFuncat$ChIPSeq_tfbD,
                                                col = heatCols$TFCol,
                                                border = T),
                           cstfbG = anno_simple(hmaFuncat$ChIPSeq_tfbG,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctbpB = anno_simple(hmaFuncat$ChIPChip_tbpB,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctbpC = anno_simple(hmaFuncat$ChIPChip_tbpC,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctbpE = anno_simple(hmaFuncat$ChIPChip_tbpE,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctbpF = anno_simple(hmaFuncat$ChIPChip_tbpF,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctfbA = anno_simple(hmaFuncat$ChIPChip_tfbA,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctfbB = anno_simple(hmaFuncat$ChIPChip_tfbB,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctfbC = anno_simple(hmaFuncat$ChIPChip_tfbC,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctfbD = anno_simple(hmaFuncat$ChIPChip_tfbD,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctfbE = anno_simple(hmaFuncat$ChIPChip_tfbE,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctfbF = anno_simple(hmaFuncat$ChIPChip_tfbF,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctfbG = anno_simple(hmaFuncat$ChIPChip_tfbG,
                                                col = heatCols$TFCol,
                                                border = T),
                           annotation_label = c("arCOG",
                                                "LSm Sense",
                                                "LSm AntiSense",
                                                "asRNA",
                                                "5' UTR Length",
                                                "5' UTR MFE",
                                                "Half life",
                                                "CAI",
                                                "GC deviation",
                                                paste0("ChIP-Seq tfb", c("B","D", "G")),
                                                paste0("ChIP-Chip tbp", c("B","C", "E", "F")),
                                                paste0("ChIP-Chip tfb", LETTERS[1:7])
                           )
)


# plotting heatmap for protein abundance #####
colors = colorRamp2(c(0, 3, 6), viridis(3))
htProt = Heatmap(log10(hmaM[,1:4]),
                 name = "Protein",
                 col = colors,
                 show_heatmap_legend = T,
                 row_names_side = "right",
                 show_row_names = F,
                 row_names_gp = gpar(fontsize = 6),
#                 row_order = hmaFuncat %>% dplyr::arrange(GCdev) %>% select(locus_tag) %>% unlist(use.names = F),
                 column_order = colnames(hmaM[,1:4]),
                 column_labels = c("TP1", "TP2", "TP3", "TP4"),
                 column_title = "Protein",
                 left_annotation = row_ha,
                 row_split = factor(hmaFuncat$arCOG),
                 #row_split = factor(hmaFuncat$asRNA),
                 #row_split = factor(hmaFuncat$lsmSense),
                 #row_split = factor(hmaFuncat$utrSize),
                 #row_split = factor(hmaFuncat$ChIPSeq_tfbD),
                 cluster_row_slices = F,
                 row_title = NULL,
                 border = T,
                 heatmap_legend_param = list(
                   title = "log10(Abund.)",
                   at = c(0, 3, 6),
                   border = "black")
)

# plotting heatmap for mrna abundance #####
htmRNA = Heatmap(log10(hmaM[,5:8]),
                 name = "mRNA",
                 col = colors,
                 show_heatmap_legend = F,
                 row_names_side = "right",
                 show_row_names = F,
                 row_names_gp = gpar(fontsize = 6),
                 column_order = colnames(hmaM[,5:8]),
                 column_labels = c("TP1", "TP2", "TP3", "TP4"),
                 column_title = "mRNA",
                 cluster_row_slices = F,
                 row_title = NULL,
                 border = T,
                 heatmap_legend_param = list(
                   title = "log10(Abund.)",
                   at = c(0, 3, 6),
                   border = "black")
)

# plotting heatmap for RPF abundance #####
htRPF = Heatmap(log10(hmaM[,9:12]),
                 name = "RPF",
                 col = colors,
                 show_heatmap_legend = F,
                 row_names_side = "right",
                 show_row_names = F,
                 row_names_gp = gpar(fontsize = 6),
                 column_order = colnames(hmaM[,9:12]),
                 column_labels = c("TP1", "TP2", "TP3", "TP4"),
                 column_title = "RPF",
                 cluster_row_slices = F,
                 row_title = NULL,
                 border = T,
                 heatmap_legend_param = list(
                   title = "log10(Abund.)",
                   at = c(0, 3, 6),
                   border = "black")
)

# plotting heatmap for TE and RO #####
colors2 = colorRamp2(c(-8, 0, 8),
                     c("#4E79A7", "white", "#E15759"))

# TE
htTE = Heatmap(log2(hmaTE),
               name = "TE",
               col = colors2,
               show_heatmap_legend = T,
               row_names_side = "right",
               show_row_names = F,
               row_names_gp = gpar(fontsize = 6),
               column_order = colnames(hmaTE),
               column_labels = c("TP1", "TP2", "TP3", "TP4"),
               column_title = "TE",
               cluster_row_slices = F,
               row_title = NULL,
               border = T,
               heatmap_legend_param = list(
                 title = "log2(Ratio)",
                 at = c(-8, 0, 8),
                 border = "black"))

# RO
htRO = Heatmap(log2(hmaRO),
               name = "RO",
               col = colors2,
               show_heatmap_legend = F,
               row_names_side = "right",
               show_row_names = T,
               row_names_gp = gpar(fontsize = 6),
               column_order = colnames(hmaRO),
               column_labels = c("TP1", "TP2", "TP3", "TP4"),
               column_title = "RO",
               cluster_row_slices = F,
               row_title = NULL,
               border = T,
               heatmap_legend_param = list(
                 title = "log2(Ratio)",
                 at = c(-8, 0, 8),
                 border = "black"))

# unifying heatmaps ####
htComplete = htProt + htmRNA + htRPF + htTE + htRO
draw(htComplete,
     annotation_legend_list = heatLegs,
     main_heatmap = "Protein")

# # # saving objx
# hma %>%
#   as.data.frame() %>%
#   write.xlsx(., file = "~/gdrive/documentos/doutorado/20201217-reuniao-vencio/dataMatrix_v2.xlsx")
# # 
# hmaFuncat %>%
#   select(-starts_with("mean")) %>%
#   as.data.frame() %>%
#   write.xlsx(., file = "~/gdrive/documentos/doutorado/20201217-reuniao-vencio/funcatMatrix_v2.xlsx")

# saving non normalized dataset for rvencio
# abund %>% 
#   dplyr::select(locus_tag,
#                 starts_with("mean_abundance_protein_lysate"),
#                 starts_with("mean_abundance_rna_total"),
#                 starts_with("mean_abundance_rna_as"),
#                 starts_with("mean_abundance_rna_ribofraction")) %>% 
#   dplyr::select(locus_tag,
#                 ends_with("TP1"),
#                 ends_with("TP2"),
#                 ends_with("TP3"),
#                 ends_with("TP4")) %>% 
#   write.xlsx(., file = "~/gdrive/documentos/doutorado/20201217-reuniao-vencio/dataMatrix_v2_nonNorm.xlsx")
