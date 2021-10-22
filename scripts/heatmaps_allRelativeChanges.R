# alorenzetti 20201104

# description ####
# this script will plot a heatmap
# based on fold change data
# for proteins, mRNAs and RPFs
# functional categorization
# here I am trying to integrate
# the time series into a single heatmap

# loading libs #####
source("scripts/loadingLibs.R")

# unifying all fold changes in a same heatmap figure #####
# preparing data
hmfc = left_join(x = alldfs$TP2 %>% dplyr::select(locus_tag, protein, mRNA, RPF),
                 y = alldfs$TP3 %>% dplyr::select(locus_tag, protein, mRNA, RPF),
                 by = "locus_tag",
                 suffix = c("_TP2", "_TP3"))

hmfc = left_join(x = hmfc,
                 y = alldfs$TP4 %>%
                   dplyr::select(locus_tag,
                                 protein_TP4 = protein,
                                 mRNA_TP4 = mRNA,
                                 RPF_TP4 = RPF),
                 by = "locus_tag") %>% 
  arrange(locus_tag)

hmfcM = hmfc %>%
  dplyr::select(-locus_tag) %>% 
  as.matrix()

rownames(hmfcM) = hmfc$locus_tag

# creating object with functional categories
hmfcFuncat = left_join(hmfc, dictFunCat,
                       by = c("locus_tag" = "pfeiLocusTag")) %>% 
  dplyr::select(-locus_tag.y)

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
  "Coenzyme transport and metabolism" = "#8CD17D", # light green
  "Defense mechanisms" = "#B6992D", # yellow green
  "Energy production and conversion" = "#F1CE63", # yellow
  "General function prediction only" = "grey70", 
  "Inorganic ion transport and metabolism" = "#86BCB6", # light teal
  "Function unknown" = "#79706E", # dark grey
  "Intracellular trafficking, secretion, and vesicular transport" = "#E15759", # red
  "Lipid transport and metabolism" = "#FF9D9A", # pink
  "Mobilome: prophages, transposons" = "#D37295", # pink
  "Nucleotide transport and metabolism" = "orchid1", # orchid1
  "Posttranslational modification, protein turnover, chaperones" = "darkturquoise", # darkturquoise
  "Replication, recombination and repair" = "skyblue2", # skyblue2
  "Secondary metabolites biosynthesis, transport and catabolism" = "#9D7660", # brown
  "Signal transduction mechanisms" = "#D7B5A6", # light orange
  "Transcription" = "#499894", # teal
  "Translation, ribosomal structure and biogenesis" = "maroon" # maroon
)

# getting classes included in arCOG
# this will be used when setting the legend for arCOG
arCOGClasses = hmfcFuncat$arCOG %>% sort() %>% unique()
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
  mfeCol = colorRamp2(breaks = c(0,hmfcFuncat$utrmfe %>% min(na.rm = T) %>% floor()),
                      colors = c("white", "#4E79A7")),
  HLCol = colorRamp2(breaks = c(0,hmfcFuncat$HL %>% max(na.rm = T) %>% ceiling()),
                     colors = c("white", "#4E79A7")),
  caiCol = colorRamp2(breaks = c(0.5,hmfcFuncat$cai %>% max(na.rm = T)),
                      colors = c("white", "#4E79A7")),
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
  TFs = Legend(title = "Transcription Factors",
               at = c("no", "yes"),
               labels = c("No", "Yes"),
               legend_gp = gpar(fill = heatCols$TFCol),
               border = "black")
)

# defining annotation columns ####
row_ha = HeatmapAnnotation(which = "row",
                           arCOG = anno_simple(hmfcFuncat$arCOG,
                                               border = T,
                                               col = heatCols$arCOGCol),
                           lsmSense = anno_simple(hmfcFuncat$lsmSense,
                                                  col = heatCols$lsmCol,
                                                  border = T),
                           lsmAntiSense = anno_simple(hmfcFuncat$lsmAntiSense,
                                                      col = heatCols$lsmCol,
                                                      border = T),
                           asRNA = anno_simple(hmfcFuncat$asRNA,
                                               col = heatCols$asRNACol,
                                               border = T),
                           utrSize = anno_simple(hmfcFuncat$utrSize,
                                                 col = heatCols$utrSizeCol,
                                                 border = T),
                           mfe = anno_simple(hmfcFuncat$utrmfe,
                                             col = heatCols$mfeCol,
                                             border = T),
                           HL = anno_simple(hmfcFuncat$HL,
                                            col = heatCols$HLCol,
                                            border = T),
                           cai = anno_simple(hmfcFuncat$cai,
                                             col = heatCols$caiCol,
                                             border = T),
                           cstfbB = anno_simple(hmfcFuncat$ChIPSeq_tfbB,
                                                col = heatCols$TFCol,
                                                border = T),
                           cstfbD = anno_simple(hmfcFuncat$ChIPSeq_tfbD,
                                                col = heatCols$TFCol,
                                                border = T),
                           cstfbG = anno_simple(hmfcFuncat$ChIPSeq_tfbG,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctbpB = anno_simple(hmfcFuncat$ChIPChip_tbpB,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctbpC = anno_simple(hmfcFuncat$ChIPChip_tbpC,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctbpE = anno_simple(hmfcFuncat$ChIPChip_tbpE,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctbpF = anno_simple(hmfcFuncat$ChIPChip_tbpF,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctfbA = anno_simple(hmfcFuncat$ChIPChip_tfbA,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctfbB = anno_simple(hmfcFuncat$ChIPChip_tfbB,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctfbC = anno_simple(hmfcFuncat$ChIPChip_tfbC,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctfbD = anno_simple(hmfcFuncat$ChIPChip_tfbD,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctfbE = anno_simple(hmfcFuncat$ChIPChip_tfbE,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctfbF = anno_simple(hmfcFuncat$ChIPChip_tfbF,
                                                col = heatCols$TFCol,
                                                border = T),
                           cctfbG = anno_simple(hmfcFuncat$ChIPChip_tfbG,
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
                                                paste0("ChIP-Seq tfb", c("B","D", "G")),
                                                paste0("ChIP-Chip tbp", c("B","C", "E", "F")),
                                                paste0("ChIP-Chip tfb", LETTERS[1:7])
                                                )
)

# plotting heatmap #####
colors = colorRamp2(c(-6, 0, 6), c("#4E79A7", "white", "#E15759"))
ht = Heatmap(hmfcM,
             name = "LFC",
             col = colors,
             show_heatmap_legend = T,
             row_names_side = "right",
             show_row_names = F,
             row_names_gp = gpar(fontsize = 6),
             column_order = colnames(hmfcM),
             column_split = factor(c(rep("TP2 vs. TP1", 3),
                                     rep("TP3 vs. TP1", 3),
                                     rep("TP4 vs. TP1", 3)),
                                   levels = c("TP2 vs. TP1",
                                              "TP3 vs. TP1",
                                              "TP4 vs. TP1")),
             column_labels = rep(c("Protein", "mRNA", "RPF"), 3),
             left_annotation = row_ha,
             row_split = factor(hmfcFuncat$arCOG),
             cluster_row_slices = F,
             row_title = NULL,
             border = T,
             heatmap_legend_param = list(
               title = "LFC",
               at = c(-6, 0, 6),
               border = T))

draw(ht, annotation_legend_list = heatLegs)

#ht_shiny(ht)
