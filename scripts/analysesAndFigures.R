# alorenzetti 202008

# description ####
# in this script we will have
# figures and additional analysis

# loading libs #####
source("scripts/loadingLibs.R")

# adding functional categorization to
# long version of tables
abundLongFuncat = left_join(abundLong,
                            dictFunCat,
                            by = c("locus_tag" = "pfeiLocusTag"))

abundNormLongFuncat = left_join(abundNormLong,
                                dictFunCat,
                                by = c("locus_tag" = "pfeiLocusTag"))

# Lsm trajectory
abundNormLongFuncat %>% 
  filter(locus_tag %in% "VNG_1496G") %>% 
  filter(libtype != "protein_ribo" & libtype != "rna_occupancy" & libtype != "rna_psiTE") %>% 
  ggplot(aes(x = timepoint,
             y = log10(mean),
             color = libtype,
             group = libtype)) +
  geom_line() +
  geom_linerange(aes(ymin = log10(mean-se), ymax = log10(mean+se))) +
  facet_wrap(~ locus_tag)

# ISH2 trajectory
abundNormLongFuncat %>% 
  filter(locus_tag %in% "VNG_0210H") %>% 
  filter(libtype != "protein_ribo" & libtype != "rna_occupancy" & libtype != "rna_psiTE") %>% 
  ggplot(aes(x = timepoint,
             y = log10(mean),
             color = libtype,
             group = libtype)) +
  geom_line() +
  geom_linerange(aes(ymin = log10(mean-se), ymax = log10(mean+se))) +
  facet_wrap(~ locus_tag)

# cluster proteinDown mRNA Up in the context of LSm protein
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(x=regRule, fill = lsmSense)) + 
  geom_bar(position = "fill") +
  coord_flip()

# cluster proteinDown mRNA Up in the context of total mRNA abundance
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(y = regRule, x = log10(totrnaBaseMean))) +
  geom_violin() +
  geom_jitter(alpha = 0.2)

# lsm in context of HL
abundNormLongFuncat %>% 
  filter(libtype == "rna_total") %>% 
  ggplot(aes(x = log10(mean), y = HL, color = lsmSense)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  xlab("log10(RNA Abundance)")

# regRule in context of HL
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(y = regRule, x = HL)) +
  geom_violin() +
  geom_jitter(alpha = 0.2)

# regRule in context of cai
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(y = regRule, x = cai)) +
  geom_violin() +
  geom_jitter(alpha = 0.2)
  
# cluster proteinDown mRNA Up trajectory
pDmUTP4 = joinedTibble[["TP4_vs_TP4"]] %>% 
  filter(regRule == "ProteinDown_mRNAUp") %>% 
  dplyr::select(locus_tag) %>% 
  unlist(use.names = F)

abundNormLongFuncat %>% 
  filter(locus_tag %in% pDmUTP4) %>% 
  filter(libtype != "protein_ribo" & libtype != "rna_occupancy" & libtype != "rna_psiTE") %>% 
  ggplot(aes(x = timepoint,
             y = log10(mean),
             color = libtype,
             group = libtype)) +
  geom_line() +
  geom_linerange(aes(ymin = log10(mean-se), ymax = log10(mean+se))) +
  facet_wrap(~ locus_tag)

# cluster mobilome in the context of LSm protein
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(x=arCOG, fill = lsmSense)) + 
  geom_bar(position = "fill") +
  coord_flip()

# cluster mobilome in the context of total mRNA abundance
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(y = arCOG, x = log10(totrnaBaseMean))) +
  geom_violin() +
  geom_jitter(alpha = 0.2)

# cluster mobilome in the context of asRNA
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(x=arCOG, fill = asRNA)) + 
  geom_bar(position = "fill") +
  coord_flip()

# cluster mobilome in the context of regRules
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(x=arCOG, fill = regRule)) + 
  geom_bar(position = "fill") +
  coord_flip()

# lsm in context of totalmRNA abundance
abundNormLongFuncat %>% 
  filter(libtype == "rna_total") %>% 
  ggplot(aes(x = log10(mean), color = lsmSense)) +
  geom_density() +
  facet_wrap(~ timepoint)

# cluster proteinDown mRNA Up trajectory
mobilomeTP4 = joinedTibble[["TP4_vs_TP4"]] %>% 
  filter(str_detect(arCOG, "^Mobilome")) %>% 
  dplyr::select(locus_tag) %>% 
  unlist(use.names = F)

abundNormLongFuncat %>% 
  filter(locus_tag %in% mobilomeTP4) %>% 
  filter(libtype != "protein_ribo" & libtype != "rna_occupancy" & libtype != "rna_psiTE") %>% 
  ggplot(aes(x = timepoint,
             y = log10(mean),
             color = libtype,
             group = libtype)) +
  geom_line() +
  geom_linerange(aes(ymin = log10(mean-se), ymax = log10(mean+se))) +
  facet_wrap(~ locus_tag)

# exploring abundance variables ####
# protein abundance in function of mRNA
# no colorspace
breaks = 10^(-10:10)
minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))

othercols = c(ggthemes_data$tableau$`color-palettes`$`ordered-sequential`$Gray$value[c(1,5,15,20) %>% rev()],
              ggthemes_data$tableau$`color-palettes`$`ordered-sequential`$Blue$value[c(1,5,15,20) %>% rev()],
              ggthemes_data$tableau$`color-palettes`$`ordered-sequential`$Green$value[c(1,20)] %>% rev())
reds = ggthemes_data$tableau$`color-palettes`$`ordered-sequential`$Red$value[c(1,5,15,20) %>% rev()]

cols = c(othercols, reds)
names(cols) = gvp1b

labs = paste0(gvp1b, " ", gvp1bnames)
names(labs) = gvp1b

abundNormLongFuncat %>% 
  dplyr::select(-se) %>% 
  filter(locus_tag %in% gvp1b) %>%
  filter(timepoint != "TP0") %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
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
#  facet_wrap(~ lsmSense) +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F),
                limits = (c(1,500000))) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F),
                limits = (c(1,500000))) +
  annotation_logticks() +
  scale_color_manual("Locus tag",
                     labels = labs,
                     values = cols)

# finding clusters using gaussian mixture model
model = Mclust(data = data.frame(mrnaTP4 = abundNorm$mean_abundance_rna_total_TP4,
                                 proteinTP4 = abundNorm$mean_abundance_protein_lysate_TP4) %>% 
                 drop_na(),
               G = 2)

data.frame(mrnaTP4 = abundNorm$mean_abundance_rna_total_TP4,
           proteinTP4 = abundNorm$mean_abundance_protein_lysate_TP4) %>% 
  as_tibble() %>% 
  drop_na() %>% 
  mutate(class = model$classification) %>% 
  ggplot(aes(x = mrnaTP4 %>% log10(),
             y = proteinTP4 %>% log10(), 
             color = class)) +
  geom_point(show.legend = F)

# protein abundance in function of mRNA
# no colorspace  
abundLongFuncat %>% 
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

# slide prot from TPi+1 mRNA from TPi
# TP2 vs TP1
slide1 = abundNorm %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP2,
             x = mean_abundance_rna_total_TP1)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks() +
  ggtitle("Protein TP2 vs. mRNA TP1")

# TP3 vs TP2
slide2 = abundNorm %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP3,
             x = mean_abundance_rna_total_TP2)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks() +
  ggtitle("Protein TP3 vs. mRNA TP2")

# TP4 vs TP3
slide3 = abundNorm %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP4,
             x = mean_abundance_rna_total_TP3)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks() +
  ggtitle("Protein TP4 vs. mRNA TP3")

# arranging plots
ggarrange(slide1, slide2, slide3,
          ncol = 3, nrow = 1)

# color == asRNA
abundNormLongFuncat %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log10(protein_lysate),
             x = log10(rna_total),
             color = asRNA)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log10(Protein Abundance)") +
  xlab("log10(mRNA Abundance)") 

# protein abundance in function of mRNA
# color == LSm
abundNormLongFuncat %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log10(protein_lysate),
             x = log10(rna_total),
             color = lsmSense)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log10(Protein Abundance)") +
  xlab("log10(mRNA Abundance)") 

# proteins in function rpfs
abundNormLongFuncat %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log10(protein_lysate),
             x = log10(rna_ribofraction))) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log10(Protein Abundance)") +
  xlab("log10(RPF Abundance)")

# what is ribosome occupancy? ####
# are they active translating ribosomes
# or stalled ribosomes?
# occupancy vs total
abundNormLongFuncat %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_occupancy),
             x = log10(rna_total),
             color = cai)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log2(Ribosome Occupancy)") +
  xlab("log10(mRNA Abundance)")  +
  scale_color_distiller(name = "cai",
                        type = "div",
                        palette = "RdYlGn")

# occupancy vs protein abundance
abundNormLongFuncat %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_occupancy),
             x = log10(protein_lysate))) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log2(Ribosome Occupancy)") +
  xlab("log10(Protein Abundance)")

# occupancy vs cai
abundNormLongFuncat %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_occupancy),
             x = cai)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log2(Ribosome Occupancy)") +
  xlab("CAI")

# what means translational efficiency? ####
# translational efficiency in function of total mRNA
abundNormLongFuncat %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_psiTE),
             x = log10(rna_total))) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log2(Protein / mRNA)") +
  xlab("log10(mRNA Abundance)")

# translational efficiency in function of protein abundance
abundNormLongFuncat %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_psiTE),
             x = log10(protein_lysate))) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log2(Protein / mRNA)") +
  xlab("log10(Protein Abundance)")

# translational efficiency in function of rpf abundance
abundNormLongFuncat %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_psiTE),
             x = log10(rna_ribofraction))) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log2(Protein / mRNA)") +
  xlab("log10(RPF Abundance)")

# translational efficiency in function of occupancy
abundNormLongFuncat %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_psiTE),
             x = log2(rna_occupancy))) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log2(Protein / mRNA)") +
  xlab("log2(Ribosome Occupancy)")
  
# checking if abundance and LSm are related
abundNormLongFuncat %>% 
  ggplot(aes(y=log10(mean), x=timepoint, fill = libtype)) +
  geom_violin() +
  facet_grid(~ lsmSense)

# checking if abundance and asRNA are related
abundNormLongFuncat %>% 
  ggplot(aes(y=log10(mean), x=timepoint, fill = libtype)) +
  geom_boxplot() +
  facet_grid(~ asRNA)

# combining both with facets
abundNormLongFuncat %>% 
  ggplot(aes(y=log10(mean), x=timepoint, fill = libtype)) +
  geom_boxplot() +
  facet_grid(lsmSense ~ asRNA)

# checking if abundance and cai are related
abundNormLongFuncat %>% 
  ggplot(aes(x = cai, y = log10(mean), color = timepoint)) +
  geom_point() +
  facet_grid(~ libtype)

# trajectories #####
# checking trajectories of ribosome proteins
riboProts = dictFunCat$pfeiLocusTag[dictFunCat$pfeiProduct %>% str_detect("ribosomal")]
abundNormLongFuncat %>% 
  filter(locus_tag %in% riboProts) %>% 
  ggplot(aes(x = timepoint,
             y = log10(mean),
             color = libtype,
             group = libtype)) +
  geom_line() +
  geom_linerange(aes(ymin = log10(mean-se), ymax = log10(mean+se))) +
  facet_wrap(~ locus_tag)

# checking trajectories of gas vesicles
gasVesProts = dictFunCat$pfeiLocusTag[dictFunCat$pfeiProduct %>% str_detect("vesicle")]
abundNormLongFuncat %>% 
  filter(locus_tag %in% gasVesProts) %>% 
  ggplot(aes(x = timepoint,
             y = log10(mean),
             color = libtype,
             group = libtype)) +
  geom_line() +
  geom_linerange(aes(ymin = log10(mean-se), ymax = log10(mean+se))) +
  facet_wrap(~ locus_tag)

# checking trajectories of mobilome proteins
abundNormLongFuncat %>% 
  filter(arCOG == "Mobilome: prophages, transposons") %>% 
  filter(libtype != "rna_psiTE" & libtype != "rna_occupancy" & libtype != "protein_ribo") %>% 
  ggplot(aes(x = timepoint,
             y = log10(mean),
             color = libtype,
             group = libtype)) +
  geom_line() +
  geom_linerange(aes(ymin = log10(mean-se), ymax = log10(mean+se))) +
  facet_wrap(~ locus_tag)

# showing that proteins with abundance < 10^2 are low abundance
# i.e. lower tail 
abundNormLongFuncat %>% 
  filter(libtype == "protein_lysate") %>% 
  filter(timepoint == "TP4") %>% 
  ggplot(aes(x = log10(mean))) + 
  geom_density() +
  geom_vline(xintercept = 2)

# showing that mRNAs with abundance > 10^4 are high abundance
# i.e. upper tail
abundNormLongFuncat %>% 
  filter(libtype == "rna_total") %>% 
  filter(timepoint == "TP4") %>% 
  ggplot(aes(x = log10(mean))) + 
  geom_density() +
  geom_vline(xintercept = 4)

# this function will take a list
# of gene names and will return a list
# of plots, one plot per gene
# it plots mRNA and protein levels
# for each time point over a distribution
# this make sense only for the quantile
# normalized dataset. ie. identical distributions
# despite the libtype
plotDistribution = function(df, geneNames){
  # creating an empty list
  plots = list()
  
  # setting up breaks for x axis
  breaks = 10^(-10:10)
  minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))
  
  # run this for each gene
  for(i in geneNames){
    dataDist = df %>%
      filter(libtype == "protein_lysate")
    
    dataAbundLevels = df %>% 
      filter(libtype == "protein_lysate" | libtype == "rna_total") %>% 
      filter(locus_tag == i)
    
    distPlot = dataDist %>% 
      ggplot(aes(x = mean)) +
      geom_density()
    
    distPlot = distPlot + 
      geom_vline(data = dataAbundLevels,
                 mapping = aes(xintercept = mean,
                               color = timepoint,
                               linetype = libtype)) +
      ylab("Density") +
      xlab("Abundance") +
      labs(color = "Time point",
           linetype = "Lib. Type") +
      scale_x_log10(breaks = breaks, minor_breaks = minor_breaks) +
      annotation_logticks(sides = "b") +
      scale_color_manual(values = c("TP1" = "#E15759",
                                    "TP2" = "#59A14F",
                                    "TP3" = "#4E79A7",
                                    "TP4" = "#B07AA1")) +
      ggtitle(i)
    
    plots[[i]] = distPlot
  }
  
  return(plots)
}

# using the function above to plot
# all mobile elements
mobilomeProts = abundNormLongFuncat %>% 
  filter(arCOG == "Mobilome: prophages, transposons") %>% 
  dplyr::select(locus_tag) %>% 
  unlist(use.names = F) %>% 
  sort() %>% 
  unique()

mobilDistPlots = plotDistribution(abundNormLongFuncat, mobilomeProts)
ggarrange(plotlist = mobilDistPlots,
          common.legend = T,
          legend = "right")

# checking relationship between occupancy and protein abundance
# seems like the higher the occupancy, the lower is the protein content
# is RPFs reflecting stale ribosomes?
# that makes sense, considering no drug was added to stop the ribosomes
abundNormFuncat %>% 
  dplyr::select(-starts_with("se")) %>% 
  ggplot(aes(x=log10(mean_abundance_protein_lysate_TP1), y=log2(mean_abundance_rna_occupancy_TP1))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ lsmSense + asRNA)

# higher occupancy also reflects in low codon adaptation index
# suggesting occupancy is actually a measure of ribosomes stalled
abundNormFuncat %>% 
  dplyr::select(-starts_with("se")) %>% 
  ggplot(aes(x=cai, y=log2(mean_abundance_rna_occupancy_TP4))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ lsmSense + asRNA)

# occupancy is also supposed to play
# a role in half lives
abundNormFuncat %>% 
  dplyr::select(-starts_with("se")) %>% 
  ggplot(aes(x=HL, y=log2(mean_abundance_rna_occupancy_TP1))) +
  geom_point() +
  geom_smooth(method = "lm")

# what about GC content?
abundNormFuncat %>% 
  dplyr::select(-starts_with("se")) %>% 
  ggplot(aes(x=GC, y=log2(mean_abundance_rna_occupancy_TP1))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ lsmSense)

# checking what are the most abundant proteins in the lysate fraction
# that is going to show us a reference
# set that can be used to compute the codon adaptation index
mostAbundProtsLysate = abundNorm %>% 
  dplyr::select(locus_tag, starts_with("mean_abundance_protein_lysate")) %>% 
  mutate(meanAcrossTP = rowMeans(.[,-1], na.rm = T)) %>% 
  dplyr::select(locus_tag, meanAcrossTP) %>% 
  mutate(meanAcrossTP = case_when(is.nan(meanAcrossTP) ~ NA_real_,
                                  TRUE ~ as.numeric(meanAcrossTP)))

fil95 = mostAbundProtsLysate$meanAcrossTP %>% quantile(probs = .95, na.rm = T)
mostAbundProtsLysate = mostAbundProtsLysate$locus_tag[mostAbundProtsLysate$meanAcrossTP > fil95]
mostAbundProtsLysate = mostAbundProtsLysate[!is.na(mostAbundProtsLysate)]
dictFunCat[dictFunCat$pfeiLocusTag %in% mostAbundProtsLysate,c("pfeiLocusTag", "pfeiProduct", "arCOGproduct")] %>% view

# what are the most abundant proteins in the ribosome fraction proteome?
mostAbundProtsRibo = abundNorm %>% 
  dplyr::select(locus_tag, starts_with("mean_abundance_protein_ribo")) %>% 
  mutate(meanAcrossTP = rowMeans(.[,-1], na.rm = T)) %>% 
  dplyr::select(locus_tag, meanAcrossTP) %>% 
  mutate(meanAcrossTP = case_when(is.nan(meanAcrossTP) ~ NA_real_,
                                  TRUE ~ as.numeric(meanAcrossTP)))

fil95 = mostAbundProtsRibo$meanAcrossTP %>% quantile(probs = .95, na.rm = T)
mostAbundProtsRibo = mostAbundProtsRibo$locus_tag[mostAbundProtsRibo$meanAcrossTP > fil95]
mostAbundProtsRibo = mostAbundProtsRibo[!is.na(mostAbundProtsRibo)]
dictFunCat[dictFunCat$pfeiLocusTag %in% mostAbundProtsRibo,c("pfeiLocusTag", "pfeiProduct", "arCOGproduct")] %>% view

# checking what are the shared proteins in
# the ribosome and lysate fractions
plot(euler(compact(list(lysate = mostAbundProtsLysate,
           ribo = mostAbundProtsRibo))),
     quantities = T)

# 18 proteins in ribo not contained in lysate
dictFunCat[dictFunCat$pfeiLocusTag %in% (mostAbundProtsRibo[!mostAbundProtsRibo %in% mostAbundProtsLysate]),]
# another question: from the set of most abundant in lysate
# how many are ribosome proteins?

# 22 proteins in lysate not contained in ribo
dictFunCat[dictFunCat$pfeiLocusTag %in% (mostAbundProtsLysate[!mostAbundProtsLysate %in% mostAbundProtsRibo]),]
# another question: from the set of most abundant in ribo
# how many are ribosome proteins?

# exploring proteome abundance in function of number of IRs ####
hmaFuncat %>%
  ggplot(aes(y=mean_abundance_protein_lysate_TP2 %>% log10,
             x=factor(x = IRs, levels = as.character(0:15)))) +
  geom_boxplot() +
  geom_quasirandom(size = .1, alpha = 0.5) +
  ylim(c(1,6)) +
  ylab("log10(Protein abundance) TP2") +
  xlab("Number of inverted repeats within sequence")
