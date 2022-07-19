# alorenzetti 20210429

# description ####
# this script will plot heat maps
# based on normalized absolute counts
# for proteins, mRNAs and RPFs

# preparing main abundance dataset ####
# dataset will be normalized within each timepoint
# establishing a function that will convert pseudocounts
# to NA before normalization
one2NA = function(x) {
  x[x == 1] = NA_integer_
  return(x)
}

# afterwards, we will also convert NA
# back to 1
NA2one = function(x) {
  x[is.na(x)] = 1
  return(x)
}

hma = abundAlt %>% 
  dplyr::select(locus_tag,
                starts_with("mean_abundance_protein_lysate"),
                starts_with("mean_abundance_rna_total"),
                starts_with("mean_abundance_rna_ribofraction")) %>% 
  dplyr::select(locus_tag,
                ends_with("TP1"),
                ends_with("TP2"),
                ends_with("TP3"),
                ends_with("TP4")) %>% 
  dplyr::mutate(across(.cols = where(is.numeric),
                       .fns = ~ one2NA(.x)))

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

# assigning a pseudocount to NA values
# in order to avoid misrepresentation in
# the heat map
hmaM[is.na(hmaM)] = 1

# reordering cols
cols = paste0(c(rep("mean_abundance_protein_lysate", 4),
                rep("mean_abundance_rna_total", 4),
                rep("mean_abundance_rna_ribofraction", 4)),
              rep(paste0("_TP", 1:4), 3))

hmaM = hmaM[,cols] %>%
  as.matrix()

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

# creating object with functional categories
hmaFuncat = left_join(hma, dictFunCat,
                      by = c("locus_tag" = "pfeiLocusTag")) %>% 
  dplyr::select(-locus_tag.y)

# defining colors and color functions ####
# 24 manual colors;
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
  "Function unknown" = "#79706E", # dark grey
  "General function prediction only" = "grey70", 
  "Inorganic ion transport and metabolism" = "#86BCB6", # light teal
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

# getting classes included in COG
# this will be used when setting the legend for COG
arCOGClasses = hmaFuncat$cog_category %>% sort() %>% unique()
arCOGcols = arCOGcols[names(arCOGcols) %in% arCOGClasses]

# defining colors for IS families
isCols = c("IS4" = "#4E79A7", # blue
           "IS5" = "#9C755F", # brown
           "IS6" = "#59A14F", # green
           "IS66" = "#FF9DA7", # pink
           "IS200/IS605" = "#E15759", #red
           "IS630" = "#B07AA1", # purple
           "IS1595" = "#EDC948", # yellow
           "ISH3" = "#F28E2B", # orange
           "ISNCY" = "#BAB0AC") # grey

isClasses = hmaFuncat$ISFamily %>% sort() %>% unique()
isCols = isCols[names(isCols) %in% isClasses]

# defining colors for the heatmap and annots
heatCols = list(
  lsmCol = c("no" = "white",
             "yes" = "#E15759"),
  arCOGCol = arCOGcols,
  isCol = isCols,
  asRNACol = c("no" = "white",
               "yes" = "#B07AA1"),
  tpsCol = c("no" = "white",
             "yes" = "#59A14F"),
  chrCol = c("no" = "white",
             "yes" = "#79706E"),
  HLCol = colorRamp2(breaks = c(0, hmaFuncat$HL %>% max(na.rm = T) %>% ceiling()),
                     colors = c("white", "#4E79A7")),
  caiCol = colorRamp2(breaks = c(0.5, hmaFuncat$cai %>% max(na.rm = T)),
                      colors = c("white", "#4E79A7")),
  # irCol = colorRamp2(breaks = c(1, hmaFuncat$IRs %>% max(na.rm = T)),
  #                    colors = c("white", "#4E79A7")),
  GCdevcol = colorRamp2(breaks = c(-0.1, 0, 0.1),
                        colors = c("#4E79A7", "white", "#E15759")),
  lfcCol = colorRamp2(breaks = c(-4, 0, 4),
                      colors = c("#4E79A7", "white", "#E15759"))
)

# defining colors and values for legends ####
heatLegs = list(
  arCOG = Legend(title = "COG",
                 at = names(heatCols$arCOGCol),
                 legend_gp  = gpar(fill = heatCols$arCOGCol %>% unname()),
                 border="black"),
  chr = Legend(title = "Location",
               at = c("no", "yes"),
               labels = c("No", "Yes"),
               legend_gp  = gpar(fill = heatCols$chrCol),
               border="black"),
  lsmSense = Legend(title = "SmAP1 Interaction",
                    at = c("no", "yes"),
                    labels = c("No", "Yes"),
                    legend_gp  = gpar(fill = heatCols$lsmCol),
                    border="black"),
  asRNA = Legend(title = "asRNA",
                 at = c("no", "yes"),
                 labels = c("No", "Yes"),
                 legend_gp  = gpar(fill = heatCols$asRNACol),
                 border="black"),
  tps = Legend(title = "TPS",
               at = c("no", "yes"),
               labels = c("No", "Yes"),
               legend_gp  = gpar(fill = heatCols$tpsCol),
               border="black"),
  lfc = Legend(title = "Mutant LFC",
               col_fun = heatCols$lfcCol,
               border = "black"),
  HL = Legend(title = "Half-life (min)",
              col_fun = heatCols$HLCol,
              border="black"),
  cai = Legend(title = "CAI",
               col_fun = heatCols$caiCol,
               border="black"),
  # ir = Legend(title = "# IRs",
  #             col_fun = heatCols$irCol,
  #             border="black"),
  GCdev = Legend(title = "GC Content",
                 col_fun = heatCols$GCdevcol,
                 border = "black")
)

# defining annotation columns ####
row_ha = HeatmapAnnotation(which = "row",
                           cog_category = anno_simple(hmaFuncat$cog_category,
                                               border = T,
                                               col = heatCols$arCOGCol),
                           Chromosome = anno_simple(hmaFuncat$Chromosome,
                                             col = heatCols$chrCol,
                                             border = T),
                           pNRC100 = anno_simple(hmaFuncat$pNRC100,
                                             col = heatCols$chrCol,
                                             border = T),
                           pNRC200 = anno_simple(hmaFuncat$pNRC200,
                                             col = heatCols$chrCol,
                                             border = T),
                           smap1Sense = anno_simple(hmaFuncat$lsmSense,
                                                  col = heatCols$lsmCol,
                                                  border = T),
                           asRNA = anno_simple(hmaFuncat$asRNA,
                                               col = heatCols$asRNACol,
                                               border = T),
                           tps = anno_simple(hmaFuncat$tps,
                                             col = heatCols$tpsCol,
                                             border = T),
                           lfc2099 = anno_simple(hmaFuncat$lfc2099,
                                               col = heatCols$lfcCol,
                                               border = T),
                           HL = anno_simple(hmaFuncat$HL,
                                            col = heatCols$HLCol,
                                            border = T),
                           cai = anno_simple(hmaFuncat$cai,
                                             col = heatCols$caiCol,
                                             border = T),
                           # ir = anno_simple(hmaFuncat$IRs,
                           #                  col = heatCols$irCol,
                           #                  border = T),
                           GCdev = anno_simple(hmaFuncat$GCdev,
                                               col = heatCols$GCdevcol,
                                               border = T),
                           annotation_label = c("COG",
                                                "Chromosome",
                                                "pNRC100",
                                                "pNRC200",
                                                "SmAP1",
                                                "asRNA",
                                                "TPS",
                                                "2099",
                                                "Half-life",
                                                "CAI",
                                                # "IRs",
                                                "GC Content"
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
                 row_split = factor(hmaFuncat$cog_category),
                 #row_split = factor(hmaFuncat$asRNA),
                 #row_split = factor(hmaFuncat$lsmSense),
                 cluster_row_slices = F,
                 row_title = NULL,
                 border = T,
                 heatmap_legend_param = list(
                   title = expression(Log[10](Abund.)),
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
                   title = expression(Log[10](Abund.)),
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
                  title = expression(Log[10](Abund.)),
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
                 title = expression(Log[2](Ratio)),
                 at = c(-8, 0, 8),
                 border = "black"))

# RO
htRO = Heatmap(log2(hmaRO),
               name = "RO",
               col = colors2,
               show_heatmap_legend = F,
               row_names_side = "right",
               show_row_names = F,
               row_names_gp = gpar(fontsize = 6),
               column_order = colnames(hmaRO),
               column_labels = c("TP1", "TP2", "TP3", "TP4"),
               column_title = "RO",
               cluster_row_slices = F,
               row_title = NULL,
               border = T,
               heatmap_legend_param = list(
                 title = expression(Log[2](Ratio)),
                 at = c(-8, 0, 8),
                 border = "black"))

# unifying heatmaps and saving the compact version of figure ####
htComplete = htProt + htmRNA + htRPF + htTE + htRO

# saving
pdf(file = "plots/abundanceHeatmap_en.pdf",
       width = 12.5,
       height = 7)
draw(htComplete,
     annotation_legend_list = heatLegs,
     main_heatmap = "Protein")
dev.off()

tiff(file = "plots/abundanceHeatmap_en.tiff",
    width = 12.5,
    height = 7,
    units = "in",
    res = 600)
draw(htComplete,
     annotation_legend_list = heatLegs,
     main_heatmap = "Protein")
dev.off()

tiff(file = "plots/abundanceHeatmap_en.png",
     width = 12.5,
     height = 7,
     units = "in",
     res = 600)
draw(htComplete,
     annotation_legend_list = heatLegs,
     main_heatmap = "Protein")
dev.off()

# saving the expanded version of figure with gene names and products ####
htRO = Heatmap(log2(hmaRO),
               name = "RO",
               col = colors2,
               show_heatmap_legend = F,
               row_names_side = "right",
               row_labels = paste0(hmaFuncat$locus_tag, "; ", hmaFuncat$product),
               show_row_names = T,
               row_names_gp = gpar(fontsize = 6),
               row_names_max_width = unit(10, "cm"),
               column_order = colnames(hmaRO),
               column_labels = c("TP1", "TP2", "TP3", "TP4"),
               column_title = "RO",
               cluster_row_slices = F,
               row_title = NULL,
               border = T,
               heatmap_legend_param = list(
                 title = expression(Log[2](Ratio)),
                 at = c(-8, 0, 8),
                 border = "black"))

htComplete = htProt + htmRNA + htRPF + htTE + htRO

pdf(file = "plots/abundanceHeatmap_expanded_en.pdf",
    width = 15,
    height = 180)
draw(htComplete,
     annotation_legend_list = heatLegs,
     main_heatmap = "Protein")
dev.off()

# saving the expanded version for 269 putative
# post-transcriptionally regulated genes
# pdf(file = "plots/abundanceHeatmap_expanded_269_genes_en.pdf",
#     width = 15,
#     height = 40)
# draw(htComplete[hmaFuncat$locus_tag %in%
#                   unique(c(ptgs$union$Q4$locus_tag,
#                            ptgsAbund$union$prot_bot_prot_non_mrna_top))],
#      annotation_legend_list = heatLegs,
#      main_heatmap = "Protein")
# dev.off()

# preparing a version for interactive heatmap
htRO = Heatmap(log2(hmaRO),
               name = "RO",
               col = colors2,
               show_heatmap_legend = F,
               row_names_side = "right",
               show_row_names = F,
               row_labels = hmaFuncat$locus_tag,
               row_names_gp = gpar(fontsize = 10),
               row_names_max_width = unit(10, "cm"),
               column_order = colnames(hmaRO),
               column_labels = c("TP1", "TP2", "TP3", "TP4"),
               column_title = "RO",
               cluster_row_slices = F,
               row_title = NULL,
               border = T,
               heatmap_legend_param = list(
                 title = expression(Log[2](Ratio)),
                 at = c(-8, 0, 8),
                 border = "black"))

htComplete = htProt + htmRNA + htRPF + htTE + htRO

# saving an object to use in shiny
save(htComplete,
     file = "results/ht_with_names_for_shiny.RData")

# saving legends to use in shiny
save(heatLegs,
     file = "results/heat_legs_for_shiny.RData")

png(file = paste0("plots/heatmap_legends.png"),
    width = 25,
    height = 12.5,
    units = "cm",
    res = 600)
draw(packLegend(list = heatLegs,
                max_height = unit(10, "cm"),
                column_gap = unit(1, "cm")))
dev.off()

# storing and saving subset versions
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
                 row_split = factor(hmaFuncat$cog_category),
                 #row_split = factor(hmaFuncat$asRNA),
                 #row_split = factor(hmaFuncat$lsmSense),
                 cluster_row_slices = F,
                 row_title = NULL,
                 border = T,
                 heatmap_legend_param = list(
                   title = expression(Log[10](Abund.)),
                   at = c(0, 3, 6),
                   border = "black",
                   title_position = "lefttop",
                   direction = "horizontal")
)

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
                 title = expression(Log[2](Ratio)),
                 at = c(-8, 0, 8),
                 border = "black",
                 title_position = "lefttop",
                 direction = "horizontal"))

htRO = Heatmap(log2(hmaRO),
               name = "RO",
               col = colors2,
               show_heatmap_legend = F,
               row_names_side = "right",
               row_labels = paste0(hmaFuncat$locus_tag, "; ", hmaFuncat$product),
               show_row_names = T,
               row_names_gp = gpar(fontsize = 10),
               row_names_max_width = unit(10, "cm"),
               column_order = colnames(hmaRO),
               column_labels = c("TP1", "TP2", "TP3", "TP4"),
               column_title = "RO",
               cluster_row_slices = F,
               row_title = NULL,
               border = T,
               heatmap_legend_param = list(
                 title = expression(Log[2](Ratio)),
                 at = c(-8, 0, 8),
                 border = "black"))

htComplete = htProt + htmRNA + htRPF + htTE + htRO

# including only putatively ptgs genes
htptgs = list()
# if(!dir.exists("plots/abundanceHeatmap_ptgs/")){dir.create("plots/abundanceHeatmap_ptgs/")}
# for(i in names(ptgs)){
#   for(j in names(ptgs[[i]])){
#     htptgs[[i]][[j]] = htComplete[hmaFuncat$locus_tag %in% ptgs[[i]][[j]]$locus_tag]
#     
#     if(length(ptgs[[i]][[j]]$locus_tag) > 0){
#       fight = 1 + length(ptgs[[i]][[j]]$locus_tag)/3
#       
#       png(file = paste0("plots/abundanceHeatmap_ptgs/", i, "_", j, ".png"),
#            width = 12,
#            height = fight,
#            units = "in",
#            res = 600)
#       draw(htptgs[[i]][[j]],
#            show_heatmap_legend = F,
#            column_title = paste(i, j))
#       dev.off()
#     }
#   }
# }

# saving the heatmap for TP1 vs TP2
htptgs$TP2_vs_TP1$Q4 = htComplete[hmaFuncat$locus_tag %in% ptgs$TP2_vs_TP1$Q4$locus_tag]
categories = hmaFuncat %>%
  dplyr::filter(locus_tag %in% ptgs$TP2_vs_TP1$Q4$locus_tag) %>% 
  pull(cog_category) %>% 
  unique()
indexes = which(names(heatCols$arCOGCol) %in% categories)

heatLegsSub = heatLegs
heatLegsSub$arCOG = Legend(title = "COG",
                           at = names(heatCols$arCOGCol)[indexes],
                           legend_gp  = gpar(fill = heatCols$arCOGCol %>% unname() %>% .[indexes]),
                           border="black")

png(file = paste0("plots/tp1_vs_tp2_q4_heatmap.png"),
    width = 14,
    height = 3.5,
    units = "in",
    res = 600)
draw(htptgs$TP2_vs_TP1$Q4,
     heatmap_legend_side = "top",
     annotation_legend_list = heatLegsSub,
     annotation_legend_side = "bottom")
dev.off()

# checking out gvp1 cluster (VNG_7015-7028)
htgvp = htComplete[hmaFuncat$locus_tag %in% paste0("VNG_", 7015:7028)]

categories = hmaFuncat %>%
  dplyr::filter(locus_tag %in% paste0("VNG_", 7015:7028)) %>% 
  pull(cog_category) %>% 
  unique()
indexes = which(names(heatCols$arCOGCol) %in% categories)

heatLegsSub = heatLegs
heatLegsSub$arCOG = Legend(title = "COG",
                           at = names(heatCols$arCOGCol)[indexes],
                           legend_gp  = gpar(fill = heatCols$arCOGCol %>% unname() %>% .[indexes]),
                           border="black")

png(file = paste0("plots/gvp_heat.png"),
    width = 14,
    height = 5,
    units = "in",
    res = 600)
draw(htgvp,
     heatmap_legend_side = "top",
     annotation_legend_list = heatLegsSub,
     annotation_legend_side = "bottom")
dev.off()

# # arranging and saving a version ordered by GC ####
# htProt = Heatmap(log10(hmaM[,1:4]),
#                  name = "Protein",
#                  col = colors,
#                  show_heatmap_legend = T,
#                  row_names_side = "right",
#                  show_row_names = F,
#                  row_names_gp = gpar(fontsize = 6),
#                  row_order = hmaFuncat %>% dplyr::arrange(GCdev) %>% select(locus_tag) %>% unlist(use.names = F),
#                  column_order = colnames(hmaM[,1:4]),
#                  column_labels = c("TP1", "TP2", "TP3", "TP4"),
#                  column_title = "Protein",
#                  left_annotation = row_ha,
#                  cluster_row_slices = F,
#                  row_title = NULL,
#                  border = T,
#                  heatmap_legend_param = list(
#                    title = expression(Log[10](Abund.)),
#                    at = c(0, 3, 6),
#                    border = "black")
# )
# 
# htRO = Heatmap(log2(hmaRO),
#                name = "RO",
#                col = colors2,
#                show_heatmap_legend = F,
#                row_names_side = "right",
#                show_row_names = F,
#                row_names_gp = gpar(fontsize = 6),
#                column_order = colnames(hmaRO),
#                column_labels = c("TP1", "TP2", "TP3", "TP4"),
#                column_title = "RO",
#                cluster_row_slices = F,
#                row_title = NULL,
#                border = T,
#                heatmap_legend_param = list(
#                  title = expression(Log[2](Ratio)),
#                  at = c(-8, 0, 8),
#                  border = "black"))
# 
# htComplete = htProt + htmRNA + htRPF + htTE + htRO
# 
# # saving
# pdf(file = "plots/abundanceHeatmap_orderedGCdev_en.pdf",
#     width = 11.5,
#     height = 7)
# draw(htComplete,
#      annotation_legend_list = heatLegs,
#      main_heatmap = "Protein")
# dev.off()

# creating and saving a version emphasizing mobilome cluster ####
# adding IS family info for this block only
row_ha2 = HeatmapAnnotation(which = "row",
                            isFam = anno_simple(hmaFuncat$ISFamily,
                                                border = T,
                                                col = heatCols$isCol),
                            annotation_label = "IS Family")

row_haMob = c(row_ha, row_ha2)

# making a simplified version of legends because
# not all need legends again
simp_legs = list(lfc = heatLegs$lfc,
                 isFam = Legend(title = "IS Family",
                                at = names(heatCols$isCol),
                                legend_gp  = gpar(fill = heatCols$isCol %>% unname()),
                                border="black"))

htProt = Heatmap(log10(hmaM[,1:4]),
                 name = "Protein",
                 col = colors,
                 show_heatmap_legend = T,
                 row_names_side = "right",
                 show_row_names = F,
                 row_names_gp = gpar(fontsize = 6),
                 #row_order = hmaFuncat %>% dplyr::arrange(GCdev) %>% select(locus_tag) %>% unlist(use.names = F),
                 column_order = colnames(hmaM[,1:4]),
                 column_labels = c("TP1", "TP2", "TP3", "TP4"),
                 column_title = "Protein",
                 left_annotation = row_haMob,
                 row_split = factor(hmaFuncat$cog_category),
                 #row_split = factor(hmaFuncat$asRNA),
                 #row_split = factor(hmaFuncat$lsmSense),
                 cluster_row_slices = F,
                 row_title = NULL,
                 border = T,
                 heatmap_legend_param = list(
                   title = expression(Log[10](Abund.)),
                   at = c(0, 3, 6),
                   border = "black",
                   title_position = "lefttop",
                   direction = "horizontal")
)

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
                 title = expression(Log[2](Ratio)),
                 at = c(-8, 0, 8),
                 border = "black",
                 title_position = "lefttop",
                 direction = "horizontal"))

htRO = Heatmap(log2(hmaRO),
               name = "RO",
               col = colors2,
               show_heatmap_legend = F,
               row_names_side = "right",
               show_row_names = T,
               row_names_gp = gpar(fontsize = 10),
               column_order = colnames(hmaRO),
               column_labels = c("TP1", "TP2", "TP3", "TP4"),
               column_title = "RO",
               cluster_row_slices = F,
               row_title = NULL,
               border = T,
               heatmap_legend_param = list(
                 title = expression(Log[2](Ratio)),
                 at = c(-8, 0, 8),
                 border = "black"))

htComplete = htProt + htmRNA + htRPF + htTE + htRO
htComplete = htComplete[hmaFuncat$cog_category == "Mobilome: prophages, transposons",]
htComplete = grid.grabExpr(draw(htComplete,
                                heatmap_legend_side = "top",
                                annotation_legend_side = "right",
                                annotation_legend_list = simp_legs),
                           wrap = T)

# ggarrange(plotlist = list(htComplete))

# comparing mobile elements to everything else 
pcomp = list()
al = 0.1
sz = 0.3

# protein levels
pcomp[["prot"]] = hmaFuncat %>% 
  mutate(cog_category = case_when(cog_category != "Mobilome: prophages, transposons" ~ "Other classes",
                                  TRUE ~ "Mobilome")) %>% 
  rowwise() %>% 
  mutate(prot = mean(mean_abundance_protein_lysate_TP1,
                     mean_abundance_protein_lysate_TP2,
                     mean_abundance_protein_lysate_TP3,
                     mean_abundance_protein_lysate_TP4)) %>% 
  ggplot(aes(y = log10(prot), x = cog_category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size = sz, alpha = al) +
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
  theme_pubr() +
  ylab("Log<sub>10</sub>(Protein abundance)") +
  xlab(NULL) +
  ylim(c(0.001,6)) +
  #  ggtitle("B") +
  theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_markdown())

# protein levels
pcomp[["mrna"]] = hmaFuncat %>% 
  mutate(cog_category = case_when(cog_category != "Mobilome: prophages, transposons" ~ "Other classes",
                                  TRUE ~ "Mobilome")) %>% 
  rowwise() %>% 
  mutate(mrna = mean(mean_abundance_rna_total_TP1,
                     mean_abundance_rna_total_TP2,
                     mean_abundance_rna_total_TP3,
                     mean_abundance_rna_total_TP4)) %>% 
  ggplot(aes(y = log10(mrna), x = cog_category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size = sz, alpha = al) +
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
  theme_pubr() +
  ylab("Log<sub>10</sub>(mRNA abundance)") +
  xlab(NULL) +
  #  ggtitle("C") +
  theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_markdown())

# TE
pcomp[["TE"]] = hmaFuncat %>% 
  mutate(cog_category = case_when(cog_category != "Mobilome: prophages, transposons" ~ "Other classes",
                                  TRUE ~ "Mobilome")) %>% 
  rowwise() %>% 
  mutate(TE = mean(mean_abundance_protein_lysate_TP1 / mean_abundance_rna_total_TP1,
                   mean_abundance_protein_lysate_TP2 / mean_abundance_rna_total_TP2,
                   mean_abundance_protein_lysate_TP3 / mean_abundance_rna_total_TP3,
                   mean_abundance_protein_lysate_TP4 / mean_abundance_rna_total_TP4)) %>% 
  ggplot(aes(y = log2(TE), x = cog_category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size = sz, alpha = al) +
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
  theme_pubr() +
  ylab("Log<sub>2</sub>(TE)") +
  xlab(NULL) +
  #  ggtitle("D") +
  theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_markdown())

# codon adaptation index 
pcomp[["CAI"]] = hmaFuncat %>% 
  mutate(cog_category = case_when(cog_category != "Mobilome: prophages, transposons" ~ "Other classes",
                                  TRUE ~ "Mobilome")) %>% 
  ggplot(aes(y = cai, x = cog_category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size = sz, alpha = al) +
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
  theme_pubr() +
  ylab("CAI") +
  xlab(NULL) +
  #  ggtitle("E") +
  theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_markdown())

# GC
pcomp[["GC"]] = hmaFuncat %>% 
  mutate(cog_category = case_when(cog_category != "Mobilome: prophages, transposons" ~ "Other classes",
                                  TRUE ~ "Mobilome")) %>% 
  ggplot(aes(y = GC, x = cog_category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size = sz, alpha = al) +
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
  theme_pubr() +
  ylab("GC content") +
  xlab(NULL) +
  #  ggtitle("F") +
  theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_markdown())

# arranging plots
panelHeatmap = ggarrange(plotlist = list(htComplete),
                         labels = "AUTO")

panelBoxplots = ggarrange(plotlist = list(pcomp$GC,
                                          pcomp$prot,
                                          pcomp$CAI,
                                          pcomp$mrna),
                          nrow = 1,
                          ncol = 4,
                          labels = LETTERS[2:5])

finalPanel = ggarrange(plotlist = list(panelHeatmap,
                                       panelBoxplots),
                       nrow = 2,
                       heights = c(1.5,1))

ggsave(filename = "plots/mobileElPanelFeatures_en.png",
       plot = finalPanel,
       units = "in",
       width = 8.5,
       height = 10,
       dpi = 300)
