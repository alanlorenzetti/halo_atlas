# alorenzetti 20201118

# description ####
# this script will generate
# figures for the final bepe report

# loading libs ####
source("scripts/loadingLibs.R")

# plotting ####
# general abundance heatmap ####
# saving the abundance heatmap
svglite(file = "~/gdrive/documentos/fapesp/bepe/report/img/abundanceHeatmap.svg",
        width = 15.5,
        height = 8)
draw(htComplete,
     annotation_legend_list = heatLegs,
     main_heatmap = "Protein")
dev.off()

# GCdev ordinated heatmap ####
# saving the abundance heatmap
svglite(file = "~/gdrive/documentos/fapesp/bepe/report/img/abundanceHeatmapGCdevOrd.svg",
        width = 15.5,
        height = 8)
draw(htComplete,
     annotation_legend_list = heatLegs,
     main_heatmap = "Protein")
dev.off()

# zoomed heatmap ####
# creating zoomed fig
svglite(file = "~/gdrive/documentos/fapesp/bepe/report/img/abundanceHeatmapZoomed.svg",
        width = 10,
        height = 90)
draw(htComplete,
     main_heatmap = "Protein")
dev.off()

# general protein-mRNA trends plot ####
breaks = 10^(-10:10)
minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))

# protein abundance in function of mRNA ####
# no colorspace  
protvsmrna = abundNormLongFuncat %>% 
  filter(timepoint != "TP0") %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = protein_lysate,
             x = rna_total)) +
  #    geom_density2d() +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks()

# saving
ggsave(filename = "~/gdrive/documentos/fapesp/bepe/report/img/mRNAproteinCorrelation.svg",
       plot = protvsmrna,
       units = "in",
       width = 6,
       height = 5)

# comparing mobile elements to everything else ####
pcomp = list()
al = 0.25
sz = 0.3
# protein levels
pcomp[["prot"]] = hmaFuncat %>% 
  mutate(arCOG = case_when(arCOG != "Mobilome: prophages, transposons" ~ "Other classes",
                           TRUE ~ "Mobilome")) %>% 
  rowwise() %>% 
  mutate(prot = mean(mean_abundance_protein_lysate_TP1,
                     mean_abundance_protein_lysate_TP2,
                     mean_abundance_protein_lysate_TP3,
                     mean_abundance_protein_lysate_TP4)) %>% 
  ggplot(aes(y = log10(prot), x = arCOG)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = sz, alpha = al) +
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
  theme_pubr() +
  ylab("log10(Protein level)") +
  xlab(NULL) +
  ggtitle("B") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# protein levels
pcomp[["mrna"]] = hmaFuncat %>% 
  mutate(arCOG = case_when(arCOG != "Mobilome: prophages, transposons" ~ "Other classes",
                           TRUE ~ "Mobilome")) %>% 
  rowwise() %>% 
  mutate(mrna = mean(mean_abundance_rna_total_TP1,
                     mean_abundance_rna_total_TP2,
                     mean_abundance_rna_total_TP3,
                     mean_abundance_rna_total_TP4)) %>% 
  ggplot(aes(y = log10(mrna), x = arCOG)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = sz, alpha = al) +
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
  theme_pubr() +
  ylab("log10(mRNA level)") +
  xlab(NULL) +
  ggtitle("C") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# TE
pcomp[["TE"]] = hmaFuncat %>% 
  mutate(arCOG = case_when(arCOG != "Mobilome: prophages, transposons" ~ "Other classes",
                           TRUE ~ "Mobilome")) %>% 
  rowwise() %>% 
  mutate(TE = mean(mean_abundance_protein_lysate_TP1 / mean_abundance_rna_total_TP1,
                   mean_abundance_protein_lysate_TP2 / mean_abundance_rna_total_TP2,
                   mean_abundance_protein_lysate_TP3 / mean_abundance_rna_total_TP3,
                   mean_abundance_protein_lysate_TP4 / mean_abundance_rna_total_TP4)) %>% 
  ggplot(aes(y = log2(TE), x = arCOG)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = sz, alpha = al) +
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
  theme_pubr() +
  ylab("log2(Translational Efficiency)") +
  xlab(NULL) +
  ggtitle("D") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# codon adaptation index 
pcomp[["CAI"]] = hmaFuncat %>% 
  mutate(arCOG = case_when(arCOG != "Mobilome: prophages, transposons" ~ "Other classes",
                           TRUE ~ "Mobilome")) %>% 
  ggplot(aes(y = cai, x = arCOG)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = sz, alpha = al) +
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
  theme_pubr() +
  ylab("Codon Adaptation Index") +
  xlab(NULL) +
  ggtitle("E") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# GC
pcomp[["GC"]] = hmaFuncat %>% 
  mutate(arCOG = case_when(arCOG != "Mobilome: prophages, transposons" ~ "Other classes",
                           TRUE ~ "Mobilome")) %>% 
  ggplot(aes(y = GC, x = arCOG)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = sz, alpha = al) +
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
  theme_pubr() +
  ylab("GC content") +
  xlab(NULL) +
  ggtitle("F") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# arranging plots
panel = ggarrange(plotlist = pcomp,
                  nrow = 1,
                  ncol = 5)

ggsave("~/gdrive/documentos/fapesp/bepe/report/img/mobileElFeatComparison.svg")


