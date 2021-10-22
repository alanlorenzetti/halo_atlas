# 20210813

# description ####
# a couple of charts to show the correlation between 
# transcripts, RPF, and protein levels. Ideally, we should be able to
# identify cases of transcriptional regulation, translational efficiency
# and translational regulation

# getting started ####
# adding functional categorization to
# long version of tables
abundLongFuncat = left_join(abundLong,
                            dictFunCat,
                            by = c("locus_tag" = "pfeiLocusTag"))

abundDerLongFuncat = left_join(abundDerLong,
                               dictFunCat,
                               by = c("locus_tag" = "pfeiLocusTag"))


abundNormLongFuncat = left_join(abundNormLong,
                                dictFunCat,
                                by = c("locus_tag" = "pfeiLocusTag"))

abundNormDerLongFuncat = left_join(abundNormDerLong,
                                   dictFunCat,
                                   by = c("locus_tag" = "pfeiLocusTag"))


# plotting
breaks = 10^(-10:10)
minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))

# axis names
yname = "Protein Abundance"
xname = "mRNA Abundance (TPM)"

# protein abundance in function of mRNA ####
# plot function
# one should provide tp for protein and for mRNA
plotpvsm = function(df, tpprot, tpmrna){
  # getting ptgsAbund locus tags
  prot_top_mrna_bot = base::intersect(ptgsAbund[[tpprot]][["prot_top"]], ptgsAbund[[tpmrna]][["mrna_bot"]])
  prot_bot_mrna_top = base::intersect(ptgsAbund[[tpprot]][["prot_bot"]], ptgsAbund[[tpmrna]][["mrna_top"]])
  
  # setting ggtitle
  if(tpprot == tpmrna){
    tptitle = tpmrna
  } else {
    tptitle = paste("Protein", tpprot, "vs.",
                     "mRNA", tpmrna)
  }
  
  # setting colnames
  protcol = paste0("mean_abundance_protein_lysate_", tpprot)
  mrnacol = paste0("mean_abundance_rna_total_", tpmrna)
  
  plot = df %>% 
    mutate(ptgsclass = case_when(locus_tag %in% prot_top_mrna_bot ~ "prot_top_mrna_bot",
                                 locus_tag %in% prot_bot_mrna_top ~ "prot_bot_mrna_top",
                                 TRUE ~ "regular")) %>% 
    ggplot(aes(y = get(protcol),
               x = get(mrnacol),
               color = ptgsclass)) +
    geom_point(alpha = 0.75, stroke = 0,
               show.legend = F) +
    geom_smooth(aes(y = get(protcol),
                    x = get(mrnacol)),
                method = "lm",
                color = tab10$value[1]) +
    stat_cor(aes(y = get(protcol),
                 x = get(mrnacol)),
             inherit.aes = F,
             method = "pearson") +
    ylab(yname) +
    xlab(xname) +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown()) +
    scale_x_log10(breaks = breaks,
                  minor_breaks = minor_breaks,
                  labels = function(x) format(x, scientific = F)) +
    scale_y_log10(breaks = breaks,
                  minor_breaks = minor_breaks,
                  labels = function(x) format(x, scientific = F)) +
    annotation_logticks() +
    scale_color_manual(values = c("prot_top_mrna_bot" = "#F28E2B",
                                  "prot_bot_mrna_top" = "#59A14F",
                                  "regular" = "black")) +
    ggtitle(tptitle)
  
  return(plot)
}

# getting plots
protvsmrna = list()
for(i in abundLongFuncat$timepoint %>% unique){
  protvsmrna[[i]] = plotpvsm(abund, i, i)
}

protvsmrnapanel = ggarrange(plotlist = protvsmrna,
                            ncol = 2, nrow = 2)

# saving
ggsave(filename = "plots/11_prot_vs_mrna_tpwise.png",
       plot = protvsmrnapanel,
       units = "in",
       width = 6,
       height = 5)

# protein levels seem to be related
# to GC content
protvsgc = abundDerLongFuncat %>% 
  filter(timepoint != "TP0") %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = protein_lysate,
             x = GC,
             color = cai)) +
  ylab("Protein Abundance") +
  xlab("GC content") +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown()) +
  geom_point(alpha = 0.25, stroke = 0) +
  geom_smooth(method = "loess",
              color = "black",
              alpha = 0.75) +
  facet_wrap(~ timepoint) +
  geom_vline(xintercept = mean(dictFunCat$GC),
             linetype = "dashed",
             size = 1,
             alpha = 0.75) +
  scale_color_viridis(name = "CAI") +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  annotation_logticks(sides = "l") +
  xlim(c(0.39, 0.8))

# saving
ggsave(filename = "plots/13_prot_vs_gc_tpwise.png",
       plot = protvsgc,
       units = "in",
       width = 6,
       height = 5)

# slide prot from TPi+1 mRNA from TPi
# TP2 vs TP1
slide1 = plotpvsm(abund, "TP2", "TP1")

# TP3 vs TP2
slide2 = plotpvsm(abund, "TP3", "TP2")

# TP4 vs TP3
slide3 = plotpvsm(abund, "TP4", "TP3")

# arranging plots
protvsmrnaslide = ggarrange(slide1, slide2, slide3,
                            ncol = 2, nrow = 2)

# saving
ggsave(filename = "plots/12_prot_vs_mrna_slides.png",
       plot = protvsmrnaslide,
       units = "in",
       width = 7,
       height = 6)

# arranging combined version
# of tp-wise plots and
# time-lag plots
protvsmrnatpwisetimelag = ggarrange(plotlist = list(protvsmrna$TP1,
                                                    protvsmrna$TP2,
                                                    protvsmrna$TP3,
                                                    protvsmrna$TP4,
                                                    slide1,
                                                    slide2,
                                                    slide3),
                                    ncol = 2, nrow = 4,
                                    labels = "AUTO")

ggsave(filename = "plots/14_prot_vs_mrna_panel.png",
       plot = protvsmrnatpwisetimelag,
       units = "in",
       dpi = 300,
       width = 7,
       height = 10)

# plotting venn diagrams ####
# containing intersections of prot_bot_mrna_top
# for each time point
vennplot = list()
vennplot$venn_prot_bot_mrna_top$preplot = venn(combinations = list("TP1" = ptgsAbund$TP1$prot_bot_mrna_top,
                                                                   "TP2" = ptgsAbund$TP2$prot_bot_mrna_top,
                                                                   "TP3" = ptgsAbund$TP3$prot_bot_mrna_top,
                                                                   "TP4" = ptgsAbund$TP4$prot_bot_mrna_top))

vennplot$venn_prot_bot_mrna_top$plot = plot(vennplot$venn_prot_bot_mrna_top$preplot,
                                            fills = "white",
                                            edges = list(lex = 2))

# containing intersections of prot_bot_mrna_top
# for each time point
vennplot$venn_prot_non_mrna_top$preplot = venn(combinations = list("TP1" = ptgsAbund$TP1$prot_non_mrna_top,
                                                                   "TP2" = ptgsAbund$TP2$prot_non_mrna_top,
                                                                   "TP3" = ptgsAbund$TP3$prot_non_mrna_top,
                                                                   "TP4" = ptgsAbund$TP4$prot_non_mrna_top))

vennplot$venn_prot_non_mrna_top$plot = plot(vennplot$venn_prot_non_mrna_top$preplot,
                                            fills = "white",
                                            edges = list(lex = 2))

# containing intersections of prot_bot_prot_non_mrna_top
# for each time point
vennplot$venn_prot_bot_prot_non_mrna_top$preplot = venn(combinations = list("TP1" = ptgsAbund$TP1$prot_bot_prot_non_mrna_top,
                                                                            "TP2" = ptgsAbund$TP2$prot_bot_prot_non_mrna_top,
                                                                            "TP3" = ptgsAbund$TP3$prot_bot_prot_non_mrna_top,
                                                                            "TP4" = ptgsAbund$TP4$prot_bot_prot_non_mrna_top))

vennplot$venn_prot_bot_prot_non_mrna_top$plot = plot(vennplot$venn_prot_bot_prot_non_mrna_top$preplot,
                                                     fills = "white",
                                                     edges = list(lex = 2))

# creating a panel
vennplot$panel_abund = ggarrange(plotlist = list(vennplot$venn_prot_bot_mrna_top$plot,
                                           vennplot$venn_prot_non_mrna_top$plot,
                                           vennplot$venn_prot_bot_prot_non_mrna_top$plot),
                           labels = "AUTO",
                           nrow = 1,
                           ncol = 3)

ggsave(filename = "plots/15_prot_botOrNon_mrna_top_venn.png",
       plot = vennplot$panel_abund,
       units = "in",
       dpi = 300,
       width = 11,
       height = 4)

# venn diagrams for putative post-transcriptionally
# regulated genes using fc approach
vennplot$venn_fc_q4$preplot = venn(combinations = list("TP2 vs. TP1" = ptgs$TP2_vs_TP1$Q4$locus_tag,
                                                       "TP3 vs. TP2" = ptgs$TP3_vs_TP2$Q4$locus_tag,
                                                       "TP4 vs. TP3" = ptgs$TP4_vs_TP3$Q4$locus_tag,
                                                       "TP3 vs. TP1" = ptgs$TP3_vs_TP1$Q4$locus_tag,
                                                       "TP4 vs. TP1" = ptgs$TP4_vs_TP1$Q4$locus_tag))

vennplot$venn_fc_q4$plot = plot(vennplot$venn_fc_q4$preplot,
                                fills = "white",
                                edges = list(lex = 2))

ggsave(filename = "plots/16_prod_down_mrna_up_venn.png",
       plot = vennplot$venn_fc_q4$plot,
       units = "in",
       dpi = 300,
       width = 6.5,
       height = 6.5)

# inspecting gvp1 trajectories
# protein abundance in function of mRNA
# no colorspace
breaksbase2 = 2^(-32:32)

# declaring gvp1a
gvp1a = paste0("VNG_", 7015:7028)
gvp1anames = paste0("(Gvp", c(LETTERS[13:4], "A", "C", "N", "O"), ")")

othercols = c(ggthemes_data$tableau$`color-palettes`$`ordered-sequential`$Gray$value[c(1,5,15,20) %>% rev()],
              ggthemes_data$tableau$`color-palettes`$`ordered-sequential`$Blue$value[c(1,5,15,20) %>% rev()],
              ggthemes_data$tableau$`color-palettes`$`ordered-sequential`$Green$value[c(1,20)] %>% rev())
reds = ggthemes_data$tableau$`color-palettes`$`ordered-sequential`$Red$value[c(1,5,15,20) %>% rev()]

cols = c(othercols, reds)
names(cols) = gvp1a

labs = paste0(gvp1a, " ", gvp1anames)
names(labs) = gvp1a

# getting the locus tags
# for those proteins with at
# least one missing value (for label drawing)
mp = abund %>%
  mutate(missingprot = case_when(is.na(mean_abundance_protein_lysate_TP1) |
                                   is.na(mean_abundance_protein_lysate_TP2) |
                                   is.na(mean_abundance_protein_lysate_TP3) |
                                   is.na(mean_abundance_protein_lysate_TP4) ~ TRUE,
                                 TRUE ~ FALSE)) %>%
  pull(missingprot)

mp = abund$locus_tag[mp] %>% unique()

p_gvp_traj = abundLongFuncat %>% 
  dplyr::select(-se) %>% 
  filter(locus_tag %in% gvp1a) %>%
  filter(timepoint != "TP0") %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  mutate(repellabs = case_when(locus_tag == "VNG_7025" ~ timepoint,
                               #                               locus_tag %in% mp ~ timepoint,
                               TRUE ~ NA_character_)) %>% 
  ggplot(aes(y = protein_lysate,
             x = rna_total,
             group = locus_tag,
             color = locus_tag)) +
  geom_path(arrow = arrow(ends = "last",
                          type = "closed",
                          length = unit(0.1, "inches"),
                          angle = 20),
            size = 1.25,
            alpha = 1,
            show.legend = T) +
  geom_label_repel(aes(label = repellabs),
                   show.legend = F) +
  ylab(paste0("Log<sub>2</sub> ", yname)) +
  xlab(paste0("Log<sub>2</sub> ", xname))  +
  scale_x_continuous(breaks = breaksbase2,
                     minor_breaks = NULL,
                     trans = "log2",
                     labels = function(x) format(log2(x), scientific = F),
                     limits = (c(2^5,2^18))) +
  scale_y_continuous(breaks = breaksbase2,
                     minor_breaks = NULL,
                     trans = "log2",
                     labels = function(x) format(log2(x), scientific = F),
                     limits = (c(2^5,2^18))) +
  annotation_logticks(base = 2) +
  scale_color_manual("Locus tag",
                     labels = labs,
                     values = cols) +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown())

# saving
ggsave(filename = "plots/gvp_traj.png",
       plot = p_gvp_traj,
       width = 7,
       height = 5.5,
       units = "in",
       dpi = 300)

# plotting all genes for manual inspection: abundance #####
allLocusTags = abundLongFuncat$locus_tag %>%
  sort() %>%
  unique()

# trajectories of genes showing
# interesting patterns of translational regulation
# melements = abundLongFuncat %>% 
#   filter(cog_category == "Mobilome: prophages, transposons") %>% 
#   dplyr::select(locus_tag) %>% 
#   unlist(use.names = F) %>% 
#   unique()

# mcluster = c("VNG_0013C",
#              "VNG_0042G",
#              "VNG_0044H",
#              "VNG_0210H",
#              "VNG_2652H")
# 
# ftsz = c("VNG_0192G",
#          "VNG_0265G",
#          "VNG_0376G",
#          "VNG_1933G",
#          "VNG_6260G")

# # setting vars
# i=1 # init
# f=30 # final
# 
# # starting loop
# while(f <= length(allLocusTags)){
#   
#   if(!dir.exists("plots/abundance_v3")){
#     dir.create("plots/abundance_v3")
#   }
#   
#   protSet = allLocusTags[i:f]
#   
#   # raw variables
#   cols = c("protein_lysate" = "#E15759",
#            "rna_ribofraction" = "#59A14F",
#            "rna_total" = "#4E79A7")
#   breaks = 10^(-10:10)
#   minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))
#   
#   abundTrajectories = abundLongFuncat %>% 
#     filter(locus_tag %in% protSet) %>% 
#     filter(libtype != "protein_ribo" & libtype != "rna_occupancy" & libtype != "rna_psiTE") %>% 
#     ggplot(aes(x = timepoint,
#                y = mean,
#                color = libtype,
#                group = libtype)) +
#     geom_line() +
#     geom_linerange(aes(ymin = mean-se, ymax = mean+se)) +
#     facet_wrap(~ locus_tag) +
#     scale_y_log10(breaks = breaks,
#                   minor_breaks = minor_breaks,
#                   labels = function(x) format(log10(x), scientific = F)) +
#     ylab("Log<sub>10</sub> (Abundance)") +
#     xlab("Time point") +
#     theme(axis.title.x = element_markdown(),
#           axis.title.y = element_markdown()) +
#     annotation_logticks(sides = "l") +
#     scale_color_manual(name = "Lib. Type",
#                        values = cols,
#                        breaks = c("protein_lysate", "rna_ribofraction", "rna_total"),
#                        labels = c("Protein", "RPFs", "mRNA"))
#   
#   # saving
#   filename = paste0("plots/abundance_v3/71_", allLocusTags[i], "-", allLocusTags[f], "_abundTrajectories_tpwise.png")
#   ggsave(filename = filename,
#          plot = abundTrajectories,
#          units = "in",
#          width = 10,
#          height = 8)
#   
#   if(!dir.exists("plots/ratios_v3")){
#     dir.create("plots/ratios_v3")
#   }
#   
#   # derived variables
#   cols = c("rna_psiTE" = "#B07AA1",
#            "rna_occupancy" = "#76B7B2",
#            "rna_tlr" = "#F28E2B")
#   
#   ratioTrajectories = abundDerLong %>% 
#     filter(locus_tag %in% protSet) %>% 
#     filter(libtype != "protein_lysate") %>% 
#     filter(libtype != "rna_ribofraction") %>% 
#     filter(libtype != "rna_total") %>%
#     ggplot(aes(x = timepoint,
#                y = log2(mean),
#                color = libtype,
#                group = libtype)) +
#     geom_line() +
#     geom_linerange(aes(ymin = mean-se, ymax = mean+se)) +
#     facet_wrap(~ locus_tag) +
#     #  ylim(c(-12, 5)) +
#     ylab("log2(Ratio)") +
#     xlab("Time point") +
#     scale_color_manual(name = "Ratio",
#                        values = cols,
#                        breaks = c("rna_occupancy", "rna_psiTE", "rna_tlr"),
#                        labels = c("TE", "TC", "TLR"))
#   
#   # saving
#   filename = paste0("plots/ratios_v3/72_", allLocusTags[i], "-", allLocusTags[f], "_ratioTrajectories_tpwise.png")
#   ggsave(filename = filename,
#          plot = ratioTrajectories,
#          units = "in",
#          width = 10,
#          height = 8)
#   
#   # setting new vars
#   i=f+1
#   f=f+30
# }
