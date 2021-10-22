# alorenzetti 20210217

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
  asRNACol = c("no" = "white",
               "yes" = "#B07AA1"),
  HLCol = colorRamp2(breaks = c(0, hmaFuncat$HL %>% max(na.rm = T) %>% ceiling()),
                     colors = c("white", "#4E79A7")),
  caiCol = colorRamp2(breaks = c(0.5, hmaFuncat$cai %>% max(na.rm = T)),
                      colors = c("white", "#4E79A7")),
  GCdevcol = colorRamp2(breaks = c(-0.1, 0, 0.1),
                        colors = c("#4E79A7", "white", "#E15759"))
)

# defining colors and values for legends ####
heatLegs = list(
  arCOG = Legend(title = "COG",
                 at = names(heatCols$arCOGCol),
                 legend_gp  = gpar(fill = heatCols$arCOGCol %>% unname()),
                 border="black"),
  lsmSense = Legend(title = "Interação SmAP1",
                    at = c("no", "yes"),
                    labels = c("Não", "Sim"),
                    legend_gp  = gpar(fill = heatCols$lsmCol),
                    border="black"),
  asRNA = Legend(title = "asRNA",
                 at = c("no", "yes"),
                 labels = c("Não", "Sim"),
                 legend_gp  = gpar(fill = heatCols$asRNACol),
                 border="black"),
  HL = Legend(title = "Meia-vida (min)",
              col_fun = heatCols$HLCol,
              border="black"),
  cai = Legend(title = "CAI",
               col_fun = heatCols$caiCol,
               border="black"),
  GCdev = Legend(title = "Resíduo GC",
                 col_fun = heatCols$GCdevcol,
                 border = "black")
)

# defining annotation columns ####
row_ha = HeatmapAnnotation(which = "row",
                           arCOG = anno_simple(hmaFuncat$cog_category,
                                               border = T,
                                               col = heatCols$arCOGCol),
                           lsmSense = anno_simple(hmaFuncat$lsmSense,
                                                  col = heatCols$lsmCol,
                                                  border = T),
                           asRNA = anno_simple(hmaFuncat$asRNA,
                                               col = heatCols$asRNACol,
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
                           annotation_label = c("COG",
                                                "SmAP1",
                                                "asRNA",
                                                "Meia-vida",
                                                "CAI",
                                                "Resíduo GC"
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
                 column_title = "Proteína",
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
                 title = expression(Log[2](Razão)),
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
                 title = expression(Log[2](Razão)),
                 at = c(-8, 0, 8),
                 border = "black"))

# unifying heatmaps and saving the compact version of figure ####
htComplete = htProt + htmRNA + htRPF + htTE + htRO

# saving
if(!dir.exists("plots/tese")){dir.create("plots/tese")}
ggsave(filename = "plots/tese/abundanceHeatmap.png",
       plot = grid.grabExpr(draw(htComplete,
                                 annotation_legend_list = heatLegs,
                                 main_heatmap = "Protein")),
       width = 11.5,
       height = 7,
       dpi = 300)

# saving the expanded version of figure with gene names ####
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
                 title = expression(Log[2](Razão)),
                 at = c(-8, 0, 8),
                 border = "black"))

htComplete = htProt + htmRNA + htRPF + htTE + htRO

pdf(file = "plots/tese/abundanceHeatmap_expanded.pdf",
    width = 11.5,
    height = 100)
draw(htComplete,
     annotation_legend_list = heatLegs,
     main_heatmap = "Protein")
dev.off()

# arranging and saving a version ordered by GC ####
htProt = Heatmap(log10(hmaM[,1:4]),
                 name = "Protein",
                 col = colors,
                 show_heatmap_legend = T,
                 row_names_side = "right",
                 show_row_names = F,
                 row_names_gp = gpar(fontsize = 6),
                 row_order = hmaFuncat %>% dplyr::arrange(GCdev) %>% select(locus_tag) %>% unlist(use.names = F),
                 column_order = colnames(hmaM[,1:4]),
                 column_labels = c("TP1", "TP2", "TP3", "TP4"),
                 column_title = "Proteína",
                 left_annotation = row_ha,
                 cluster_row_slices = F,
                 row_title = NULL,
                 border = T,
                 heatmap_legend_param = list(
                   title = expression(Log[10](Abund.)),
                   at = c(0, 3, 6),
                   border = "black")
)

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
                 title = expression(Log[2](Razão)),
                 at = c(-8, 0, 8),
                 border = "black"))

htComplete = htProt + htmRNA + htRPF + htTE + htRO

# saving
ggsave(filename = "plots/tese/abundanceHeatmap_orderedGCdev.png",
       plot = grid.grabExpr(draw(htComplete,
                                 annotation_legend_list = heatLegs,
                                 main_heatmap = "Protein")),
       width = 11.5,
       height = 7,
       dpi = 300)

# creating and saving a version emphasizing mobilome cluster ####
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
                 column_title = "Proteína",
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
                 title = expression(Log[2](Razão)),
                 at = c(-8, 0, 8),
                 border = "black",
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
                 title = expression(Log[2](Razão)),
                 at = c(-8, 0, 8),
                 border = "black"))

htComplete = htProt + htmRNA + htRPF + htTE + htRO
htComplete = htComplete[hmaFuncat$cog_category == "Mobilome: prophages, transposons",]
htComplete = grid.grabExpr(draw(htComplete,
                                heatmap_legend_side = "top"))

# comparing mobile elements to everything else 
pcomp = list()
al = 0.25
sz = 0.3
# protein levels
pcomp[["prot"]] = hmaFuncat %>% 
  mutate(cog_category = case_when(cog_category != "Mobilome: prophages, transposons" ~ "Outras classes",
                           TRUE ~ "Mobiloma")) %>% 
  rowwise() %>% 
  mutate(prot = mean(mean_abundance_protein_lysate_TP1,
                     mean_abundance_protein_lysate_TP2,
                     mean_abundance_protein_lysate_TP3,
                     mean_abundance_protein_lysate_TP4)) %>% 
  ggplot(aes(y = log10(prot), x = cog_category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = sz, alpha = al) +
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
  theme_pubr() +
  ylab("Log<sub>10</sub>(Nível de proteína)") +
  xlab(NULL) +
#  ggtitle("B") +
  theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_markdown())

# protein levels
pcomp[["mrna"]] = hmaFuncat %>% 
  mutate(cog_category = case_when(cog_category != "Mobilome: prophages, transposons" ~ "Outras classes",
                           TRUE ~ "Mobiloma")) %>% 
  rowwise() %>% 
  mutate(mrna = mean(mean_abundance_rna_total_TP1,
                     mean_abundance_rna_total_TP2,
                     mean_abundance_rna_total_TP3,
                     mean_abundance_rna_total_TP4)) %>% 
  ggplot(aes(y = log10(mrna), x = cog_category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = sz, alpha = al) +
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
  theme_pubr() +
  ylab("Log<sub>10</sub>(Nível de mRNA)") +
  xlab(NULL) +
#  ggtitle("C") +
  theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_markdown())

# TE
pcomp[["TE"]] = hmaFuncat %>% 
  mutate(cog_category = case_when(cog_category != "Mobilome: prophages, transposons" ~ "Outras classes",
                           TRUE ~ "Mobiloma")) %>% 
  rowwise() %>% 
  mutate(TE = mean(mean_abundance_protein_lysate_TP1 / mean_abundance_rna_total_TP1,
                   mean_abundance_protein_lysate_TP2 / mean_abundance_rna_total_TP2,
                   mean_abundance_protein_lysate_TP3 / mean_abundance_rna_total_TP3,
                   mean_abundance_protein_lysate_TP4 / mean_abundance_rna_total_TP4)) %>% 
  ggplot(aes(y = log2(TE), x = cog_category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = sz, alpha = al) +
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
  theme_pubr() +
  ylab("Log<sub>2</sub>(TE)") +
  xlab(NULL) +
#  ggtitle("D") +
  theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_markdown())

# codon adaptation index 
pcomp[["CAI"]] = hmaFuncat %>% 
  mutate(cog_category = case_when(cog_category != "Mobilome: prophages, transposons" ~ "Outras classes",
                           TRUE ~ "Mobiloma")) %>% 
  ggplot(aes(y = cai, x = cog_category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = sz, alpha = al) +
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
  theme_pubr() +
  ylab("CAI") +
  xlab(NULL) +
#  ggtitle("E") +
  theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_markdown())

# GC
pcomp[["GC"]] = hmaFuncat %>% 
  mutate(cog_category = case_when(cog_category != "Mobilome: prophages, transposons" ~ "Outras classes",
                           TRUE ~ "Mobiloma")) %>% 
  ggplot(aes(y = GC, x = cog_category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = sz, alpha = al) +
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
  theme_pubr() +
  ylab("Conteúdo GC") +
  xlab(NULL) +
#  ggtitle("F") +
  theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_markdown())

# arranging plots
panelHeatmap = ggarrange(plotlist = list(htComplete),
                         labels = "AUTO")

panelBoxplots = ggarrange(plotlist = pcomp,
                          nrow = 1,
                          ncol = 5,
                          labels = LETTERS[2:6])

finalPanel = ggarrange(plotlist = list(panelHeatmap,
                                       panelBoxplots),
                       nrow = 2,
                       heights = c(1,1.35))

ggsave(filename = "plots/tese/mobileElPanelFeatures.png",
       plot = finalPanel,
       units = "in",
       width = 7.25,
       height = 6.6,
       dpi = 300)

# saving hmaFuncat object
# for posterior use
write_tsv(x = hmaFuncat, 
          file = "results/hmaFuncat.tsv")
