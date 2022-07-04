# 20210813

# description ####
# this script will perform most of the descriptive analysis
# in the paper. also it is going to generate
# a couple of charts to show the correlation between 
# transcripts, RPF, and protein levels

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
xname = "mRNA Abundance (TPM+1)"

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
                color = "black") +
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
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = breaks,
                  minor_breaks = minor_breaks,
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks() +
    scale_color_manual(values = c("prot_top_mrna_bot" = "#F28E2B",
                                  "prot_bot_mrna_top" = "#59A14F",
                                  "regular" = "grey80")) +
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
             alpha = 0.75) +
  scale_color_viridis(name = "CAI") +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = trans_format("log10", math_format(10^.x))) +
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
       dpi = 600,
       width = 5.5,
       height = 10)

# big panel containing
# abund abundance and fc approaches to find
# post-transcriptionally regulated genes
abundFcProtRNApanel = ggarrange(plotlist = list(ggarrange(plotlist = list(protvsmrna$TP1,
                                                                          protvsmrna$TP2,
                                                                          protvsmrna$TP3,
                                                                          protvsmrna$TP4,
                                                                          slide1,
                                                                          slide2,
                                                                          slide3,
                                                                          p_prot_mRNA$TP2_vs_TP1,
                                                                          p_prot_mRNA$TP3_vs_TP2,
                                                                          p_prot_mRNA$TP4_vs_TP3,
                                                                          p_prot_mRNA$TP3_vs_TP1,
                                                                          p_prot_mRNA$TP4_vs_TP1),
                                                          ncol = 4, nrow = 3,
                                                          labels = LETTERS[1:12]),
                                                p_prot_mRNA$LEGEND),
                                ncol = 1, nrow = 2,
                                heights = c(15,1))

ggsave(filename = "plots/18_abund_fc_prot_mrna_panel.png",
       plot = abundFcProtRNApanel,
       units = "in",
       dpi = 600,
       width = 12,
       height = 8)

# finding significance of overlap SmAP1 BRs####
# computing significance of overlap for SmAP1
# biological replicates using the rationale
# available at
# http://nemates.org/MA/progs/representation.stats.html
# or in data folder overlap_stats.v0.011.tgz
x = overlaps[["biorep1_nr"]] %in% overlaps[["biorep2_nr"]] %>% sum()
n = overlaps[["biorep1_nr"]] %>% length()
D = overlaps[["biorep2_nr"]] %>% length()
N = dim(nrtx)[1]

p = D / N
q = 1-p

# computing representation factor
RF = x / ((n * D) / N) 

# computing probability
Prob = NULL

# it is hard to get the exact hypergeometric prob
# if the following conditions is satisfied
if((p + 2*sqrt(p*q/n) > 0) & (n * 10 < N) ){
  Z = abs(((x-.5) - n * p) / sqrt(n * p * q))
  
  Prob = 1 - (gsl::erf(Z / sqrt(2)) + 1) / 2
}else{
  # get the exact hypergeometric prob
  if(RF < 1){
    for(i in 0:as.integer(x)){
      Prob[i] = exp(lnchoose(D,i) + lnchoose(N-D, n-i) - lnchoose(N, n))
    }
    Prob = sum(Prob)
  }else{
    for(i in 0:(as.integer(x)-1)){
      Prob[i] = exp(lnchoose(D,i) + lnchoose(N-D, n-i) - lnchoose(N, n))
    }
    Prob = 1 - sum(Prob)
  }
}

# p for the overlap set
Prob

# zero due to lack of precision
# the webserver reported
# Representation Factor 5.3 and p < 3.171e-83 for
# Set1: 190 
# Set2: 348 
# Overlap: 134 
# Total number of genes: 2631 

# an alternative approach to compute
# the fold enrichment with an associated
# probability taken from
# https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html
smap1PerBrList = list(BR1 = overlaps$biorep1,
                      BR2 = overlaps$biorep2)

length.gene.sets = sapply(smap1PerBrList,length)

total = dim(nrtx)[1]

num.expcted.overlap = total * do.call(prod, as.list(length.gene.sets/total))

p = sapply(0:232,function(i) dpsets(i, length.gene.sets, n=total))

common.genes = intersect(smap1PerBrList$BR1, smap1PerBrList$BR2)

num.observed.overlap = length(common.genes)

# fold enrichment of the overlap over expected
FE = num.observed.overlap/num.expcted.overlap

# probability density of the observed intersection size
dpsets(num.observed.overlap, length.gene.sets, n=total)

# probability of observing an intersection of 157 or more genes
# using the cumulative probability function
cpsets(num.observed.overlap-1, length.gene.sets, n=total, lower.tail=FALSE)

# same operations done using a more succing approach
fit=MSET(smap1PerBrList, n=total, lower.tail=FALSE)

# fold enrichment
fit$FE

# probability
fit$p.value

# comparing pearson correlations ####
comparisons = list("RTP1 vs. RTP2" = c("TP1", "TP1", "TP2", "TP2"),
                   "RTP2 vs. RTP3" = c("TP2", "TP2", "TP3", "TP3"),
                   "RTP3 vs. RTP4" = c("TP3", "TP3", "TP4", "TP4"),
                   "RTP1 vs. RTP3" = c("TP1", "TP1", "TP3", "TP3"),
                   "RTP1 vs. RTP4" = c("TP1", "TP1", "TP4", "TP4"),
                   "RTP2 vs. RTP4" = c("TP2", "TP2", "TP4", "TP4"),
                   "RTP21 vs. RTP2" = c("TP2", "TP1", "TP2", "TP2"),
                   "RTP32 vs. RTP3" = c("TP3", "TP2", "TP3", "TP3"),
                   "RTP43 vs. RTP4" = c("TP4", "TP3", "TP4", "TP4"))

corcompRes = list()

for(i in names(comparisons)){
  form = paste0("~", "mean_abundance_protein_lysate_", comparisons[[i]][1], "+",
                "mean_abundance_rna_total_", comparisons[[i]][2], "|",
                "mean_abundance_protein_lysate_", comparisons[[i]][3], "+",
                "mean_abundance_rna_total_", comparisons[[i]][4])
  
  form = as.formula(form)
  
  # cocor does not work with tibbles
  corcompRes[[i]] = cocor(formula = form,
                          data = abund %>%
                            mutate(across(.cols = starts_with("mean_abundance"),
                                          .fns = ~ log10(.x))) %>% 
                            as.data.frame())
}

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
       dpi = 600,
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
       dpi = 600,
       width = 6.5,
       height = 6.5)

# alternatively to Venn, I am also building UpSet plots #####
# declaring list
lt = list(
  # "A: TP1" = ptgsAbund$TP1$prot_bot_prot_non_mrna_top,
  # "A: TP2" = ptgsAbund$TP2$prot_bot_prot_non_mrna_top,
  # "A: TP3" = ptgsAbund$TP3$prot_bot_prot_non_mrna_top,
  # "A: TP4" = ptgsAbund$TP4$prot_bot_prot_non_mrna_top,
  "TP2 vs. TP1" = ptgs$TP2_vs_TP1$Q4$locus_tag,
  "TP3 vs. TP2" = ptgs$TP3_vs_TP2$Q4$locus_tag,
  "TP4 vs. TP3" = ptgs$TP4_vs_TP3$Q4$locus_tag,
  "TP3 vs. TP1" = ptgs$TP3_vs_TP1$Q4$locus_tag,
  "TP4 vs. TP1" = ptgs$TP4_vs_TP1$Q4$locus_tag
)

# making the combination matrix
cm1 = make_comb_mat(lt, mode = "distinct")

# plotting
usplotfc = grid.grabExpr(draw(UpSet(m = cm1,
                                    set_order = names(lt))),
                         wrap = T)

ggsave(filename = "plots/upset_plot_fcbased.png",
       plot = usplotfc,
       dpi = 600,
       units = "in",
       height = 3.5,
       width = 5)

# plotting venn diagrams for post-transcriptional ####
# related features
# table containing general counts of asRNAs, TPS, SmAP1 binding, RNase differential expression
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

lt = list("SmAP1" = summaryTPelements$SmAP1,
          "asRNA" = summaryTPelements$asRNA,
          "TPS" = c(summaryTPelements$tps$tps1,
                    summaryTPelements$tps$tps2to5,
                    summaryTPelements$tps$tps5) %>% 
            unique(),
          "RNase" = c(summaryTPelements$`2099`$up,
                            summaryTPelements$`2099`$down) %>% 
            unique())

ptgsFeaturesVenn = plot(venn(combinations = lt),
                        fills = "white",
                        edges = list(lex = 2),
                        quantities = list(cex = 2),
                        labels = list(cex = 2))

ggsave(filename = "plots/17_ptgsFeatures_venn.png",
       plot = ptgsFeaturesVenn,
       units = "in",
       dpi = 600,
       width = 8.5,
       height = 7.5)

# genes being intersected by all those features
# making the combination matrix
cm9 = make_comb_mat(lt, mode = "distinct")
# dictFunCat %>% dplyr::select(pfeiLocusTag, product) %>% filter(pfeiLocusTag %in% extract_comb(m = cm9, comb_name = "1111")) %>% clipr::write_clip()

# checking enrichment of features (eg TPS) inside other features (eg SmAP1) ####
# of post transcriptional regulation

# first question: are SmAP1-bound genes enriched for TPS? Yes
enrichAnalysis(df = hmaFuncat, subvec = summaryTPelements$SmAP1,
               testcol = "tps", lev =  "yes")

# second question: are asRNA cognate genes enriched for TPS? Yes
enrichAnalysis(df = hmaFuncat, subvec = summaryTPelements$asRNA,
               testcol = "tps", lev =  "yes")

# third question: are delta2099 differentially expressed
# genes enriched for TPS? 
# upregulated: Yes
enrichAnalysis(df = hmaFuncat, subvec = summaryTPelements$`2099`$up,
               testcol = "tps", lev =  "yes")

# downregulated: No
enrichAnalysis(df = hmaFuncat, subvec = summaryTPelements$`2099`$down,
               testcol = "tps", lev =  "yes")

# both: Yes
enrichAnalysis(df = hmaFuncat, subvec = c(summaryTPelements$`2099`$up,
                                          summaryTPelements$`2099`$down),
               testcol = "tps", lev =  "yes")

# extra questions:
# are delta2099 differentially expressed genes
# enriched for SmAP1? yes
enrichAnalysis(df = hmaFuncat, subvec = summaryTPelements$`2099`$up,
               testcol = "lsmSense", lev =  "yes")

enrichAnalysis(df = hmaFuncat, subvec = summaryTPelements$`2099`$down,
               testcol = "lsmSense", lev =  "yes")

# enriched for asRNA? yes
enrichAnalysis(df = hmaFuncat, subvec = summaryTPelements$`2099`$up,
               testcol = "asRNA", lev =  "yes")

enrichAnalysis(df = hmaFuncat, subvec = summaryTPelements$`2099`$down,
               testcol = "asRNA", lev =  "yes")

# fourth question: does the intersection of SmAP1, asRNA,
# and delta2099 DE genes is enriched for TPS? No
# for that we are going to use a function
subvectIntersectionEnrich = function(list = list){
  subvectIntersection = list()
  subvectIntersection$list = list
  subvectIntersection$matrix = list_to_matrix(subvectIntersection$list)
  subvectIntersection$vector = subvectIntersection$matrix %>%
    rowSums()
  subvectIntersection$names = subvectIntersection$vector[subvectIntersection$vector == length(subvectIntersection$list)] %>% 
    names()
  
  return(enrichAnalysis(df = hmaFuncat,
                        subvec = subvectIntersection$names,
                        testcol = "tps",
                        lev =  "yes"))
}

subvectIntersectionEnrich(list = list("SmAP1" = summaryTPelements$SmAP1,
                                      "asRNA" = summaryTPelements$asRNA,
                                      "RNase" = c(summaryTPelements$`2099`$up,
                                                  summaryTPelements$`2099`$down) %>% 
                                        unique()))

# fifth question: does the intersection of SmAP1 and asRNA 
# is enriched for TPS? Yes
subvectIntersectionEnrich(list = list("SmAP1" = summaryTPelements$SmAP1,
                                      "asRNA" = summaryTPelements$asRNA))

# sixth question: does the intersection of SmAP1 and 2099 DE genes 
# is enriched for TPS? Yes
subvectIntersectionEnrich(list = list("SmAP1" = summaryTPelements$SmAP1,
                                      "RNase" = c(summaryTPelements$`2099`$up,
                                                  summaryTPelements$`2099`$down) %>% 
                                        unique()))

# seventh question: does the intersection of asRNA and 2099 DE genes 
# is enriched for TPS? Yes
subvectIntersectionEnrich(list = list("asRNA" = summaryTPelements$asRNA,
                                      "RNase" = c(summaryTPelements$`2099`$up,
                                                  summaryTPelements$`2099`$down) %>% 
                                        unique()))

# out of abund post-transcriptionally regulated genes (248)
# how many are within the set targeted by features
# smap1, asrna, tps
# testing for enrichment
myset = ptgsAbund$union$prot_bot_prot_non_mrna_top[
  ptgsAbund$union$prot_bot_prot_non_mrna_top %in%
    (lt %>% unlist() %>% unique())
  ]

enrichAnalysis(df = hmaFuncat, subvec = myset,
               testcol = "lsmSense", lev =  "yes")

enrichAnalysis(df = hmaFuncat, subvec = myset,
               testcol = "asRNA", lev =  "yes")

enrichAnalysis(df = hmaFuncat, subvec = myset,
               testcol = "tps", lev =  "yes")

# out of all poor correlation genes (269)
# how many are within the set targeted by features
# smap1, asrna, tps
# testing for enrichment
poorCorGenes = c(ptgs$union$Q4$locus_tag, ptgsAbund$union$prot_bot_prot_non_mrna_top) %>% unique()
myset2 = poorCorGenes[
  poorCorGenes %in%
    (lt %>% unlist() %>% unique())
]

enrichAnalysis(df = hmaFuncat, subvec = myset2,
               testcol = "lsmSense", lev =  "yes")

enrichAnalysis(df = hmaFuncat, subvec = myset2,
               testcol = "asRNA", lev =  "yes")

enrichAnalysis(df = hmaFuncat, subvec = myset2,
               testcol = "tps", lev =  "yes")

# testing enrichment of functions in abundance-based ####
# genes
ptgsAbundEnrich = tibble()
for(i in hmaFuncat$cog_category %>% unique()){
  
  curTib = tibble(cog_category = i,
                  pvalue = enrichAnalysis(df = hmaFuncat,
                                          subvec = ptgsAbund$union$prot_bot_prot_non_mrna_top,
                                          testcol = "cog_category",
                                          lev = i))
  
  ptgsAbundEnrich = bind_rows(ptgsAbundEnrich,
                              curTib)
}

# filtering results
ptgsAbundEnrich = ptgsAbundEnrich %>%
  filter(pvalue < 0.05)

# inspecting gvp1 trajectories ####
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
  mutate(missingprot = case_when(is.na(mean_abundance_protein_lysate_TP1) &
                                   is.na(mean_abundance_protein_lysate_TP2) &
                                   is.na(mean_abundance_protein_lysate_TP3) &
                                   is.na(mean_abundance_protein_lysate_TP4) ~ TRUE,
                                 TRUE ~ FALSE)) %>%
  pull(missingprot)

mp = abund$locus_tag[mp] %>% unique()

mp = mp[mp %in% gvp1a]

p_gvp_traj = abundLongFuncat %>% 
  dplyr::select(-se) %>% 
  filter(locus_tag %in% gvp1a) %>%
  filter(!locus_tag %in% mp) %>% 
  filter(timepoint != "TP0") %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  mutate(repellabs = case_when(locus_tag == "VNG_7024" |
                                 locus_tag == "VNG_7023" ~ str_replace(timepoint, "TP", ""),
                               TRUE ~ NA_character_)) %>% 
  mutate(genelab = case_when(timepoint == "TP4" & locus_tag == "VNG_7023" ~ str_match(product, "Gvp.*$"),
                             timepoint == "TP1" & locus_tag == "VNG_7024" ~ str_match(product, "Gvp.*$"),
                             TRUE ~ NA_character_)) %>% 
  ggplot(aes(y = protein_lysate,
             x = rna_total,
             group = locus_tag,
             color = locus_tag)) +
  geom_point(size = 1) +
  geom_path(arrow = arrow(ends = "last",
                          type = "closed",
                          length = unit(0.1, "inches"),
                          angle = 20),
            size = 1.25,
            alpha = 1,
            show.legend = T) +
  geom_label_repel(aes(label = repellabs),
                   show.legend = F,
                   size = 2.5,
                   label.padding = 0.125,
                   label.r = 0.18,
                   min.segment.length = 0,
                   segment.colour = "black",
                   color = "black") + 
  geom_label(aes(label = genelab),
             show.legend = F,
             size = 3,
             nudge_x = 1.5) +
  # geom_label_repel(aes(label = genelab),
  #                  show.legend = F,
  #                  size = 3,
  #                  label.padding = 0.125,
  #                  label.r = 0.18,
  #                  force_pull = 0.5) +
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
                     limits = (c(2^5,2^16))) +
  annotation_logticks(base = 2) +
  scale_color_manual("Locus tag",
                     labels = labs[!names(labs) %in% mp],
                     values = cols[!names(cols) %in% mp]) +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown())

# saving
ggsave(filename = "plots/gvp_traj.png",
       plot = p_gvp_traj,
       width = 5.5,
       height = 4,
       units = "in",
       dpi = 600)

# saving
ggsave(filename = "plots/gvp_traj.svg",
       plot = p_gvp_traj,
       width = 5.5,
       height = 4,
       units = "in",
       dpi = 600)

# out of gvp1 cluster genes
# testing for enrichment
# smap1, asrna, tps
enrichAnalysis(df = hmaFuncat, subvec = gvp1a,
               testcol = "lsmSense", lev =  "yes")

enrichAnalysis(df = hmaFuncat, subvec = gvp1a,
               testcol = "asRNA", lev =  "yes")

enrichAnalysis(df = hmaFuncat, subvec = gvp1a,
               testcol = "tps", lev =  "yes")


# trajectories of SmAP1 ####
# raw variables
cols = c("protein_lysate" = "#E15759",
         "rna_ribofraction" = "#59A14F",
         "rna_total" = "#4E79A7")
breaks = 10^(-10:10)
minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))

smap1Traj = abundLongFuncat %>%
  filter(locus_tag %in% "VNG_1496G") %>%
  filter(libtype != "protein_ribo" & libtype != "rna_occupancy" & libtype != "rna_psiTE") %>%
  ggplot(aes(x = timepoint,
             y = mean,
             color = libtype,
             group = libtype)) +
  geom_line() +
  geom_linerange(aes(ymin = mean-se, ymax = mean+se)) +
  facet_wrap(~ locus_tag) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(log10(x), scientific = F)) +
  ylab("Log<sub>10</sub> (Abundance)") +
  xlab("Time point") +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown()) +
  annotation_logticks(sides = "l") +
  scale_color_manual(name = "Lib. Type",
                     values = cols,
                     breaks = c("protein_lysate", "rna_ribofraction", "rna_total"),
                     labels = c("Protein", "RPFs", "mRNA"))

# saving
ggsave(filename = "plots/smap1_traj.png",
       plot = smap1Traj,
       width = 4,
       height = 3,
       units = "in",
       dpi = 600)

# comparing GC of SmAP1-bound transcripts ####
# getting GC distributions for those
# features binding to SmAP1 vs non binding
# wilcoxon unpaired two sample test
# alt name: Mann–Whitney U test
gcComparisonPlot = hmaFuncat %>% 
  ggplot(aes(y = GC,
             x = lsmSense)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_quasirandom(size = 0.5,
                   alpha = 0.2) +
  stat_compare_means(method = "wilcox.test",
                     mapping = aes(label = ..p.signif..),
                     label.x = 1.5) +
  scale_x_discrete(labels = c("yes" = "Yes", "no" = "No")) +
  xlab(label = "SmAP1 interaction")

ggsave(filename = "plots/smap1_gc_comparison.png",
       plot = gcComparisonPlot,
       height = 4,
       width = 3,
       unit = "in",
       dpi = 600)

# plotting enrichment of COG for smap1 ####
# binding transcripts
# testing each category for enrichment
# using the protein coding transcripts
# bound to SmAP1
proteinCodingDf = hmaFuncat %>%
  filter(protein_coding == "yes") 

categories = proteinCodingDf$cog_category %>% 
  unique()

smap1EnrichRes = tibble()
for(category in categories){
  p = enrichAnalysis(df = hmaFuncat,
                     subvec = hmaFuncat %>% 
                       filter(cog_category == category) %>% 
                       pull(locus_tag),
                     testcol = "lsmSense",
                     lev = "yes")
  line = tibble(cog_category = category, pval = p)
  
  smap1EnrichRes = bind_rows(smap1EnrichRes,
                             line)
}

# correcting pvalues using BH method
# filtering by pval
smap1EnrichRes$qval = p.adjust(smap1EnrichRes$pval, method = "BH")
smap1EnrichRes = smap1EnrichRes[smap1EnrichRes$pval < 0.05,]

# plotting counts of categories with enrichment score
cogEnrichPlot = hmaFuncat %>%
  filter(lsmSense == "yes") %>%
  group_by(cog_category) %>% 
  summarise(count = n()) %>%
  left_join(x = .,
            y = smap1EnrichRes,
            by = c("cog_category")) %>%
  mutate(enrichStatus = case_when(is.na(pval) ~ NA_character_,
                                  TRUE ~ "*"),
         cog_category = factor(cog_category, levels = rev(cog_category)),
         size = case_when(count < 25 ~ "< 25 genes",
                          TRUE ~ "> 25 genes")) %>% 
  mutate(size = factor(size, levels = size %>% unique())) %>% 
  ggplot(aes(x = count,
             y = cog_category)) +
  geom_col(fill = "white", col = "black", size = 0.25) +
  geom_text(aes(label = enrichStatus),
            vjust = 0.75,
            hjust = 1.25) +
  facet_wrap(~ size,
             scales = "free_x") + 
  xlab("Observations") + 
  ylab("COG category")

ggsave(filename = "plots/smap1_cogcategory_enrich.png",
       plot = cogEnrichPlot,
       width = 7.5,
       height = 4,
       unit = "in",
       dpi = 600)

# unifying the previous three plots ####
# to make a SmAP1 panel
smap1FeaturesGCcats = ggarrange(plotlist = list(ggarrange(plotlist = list(gcComparisonPlot,
                                                                          smap1Traj),
                                                          ncol = 2, nrow = 1,
                                                          labels = c("A", "B"),
                                                          widths = c(1,1.5)),
                                                cogEnrichPlot),
                                ncol = 1, nrow = 2,
                                labels = c(" ", "C"))

ggsave(filename = "plots/smap1_panel.png",
       plot = smap1FeaturesGCcats,
       width = 7.5,
       height = 7,
       unit = "in",
       dpi = 600)

# DNA-Seq mobilization plots ####
dfins = read_tsv(file = "data/dfins.txt")
dfdel = read_tsv(file = "data/dfdel.txt")

# getting aligned reads per lib
resCounts = read_tsv(file = "data/resCounts.tsv")

# adjusting ISnames and their levels levels of factor ISname
lvs = mixedsort(levels(as.factor(dfins$ISName)))
dfins$ISName = factor(dfins$ISName, levels=rev(lvs))

lvs = mixedsort(levels(as.factor(dfdel$ISName)))
dfdel$ISName = factor(dfdel$ISName, levels=rev(lvs))

# plots
# how many IS per barcode
# insertions
isCountPerLib = list()
isCountPerLib$ins = ggplot(dfins,
                           (aes(ISName, fill = ISFamily))) +
  geom_bar(show.legend = F) +
  ylim(c(0, 30)) +
  xlab(label = "") +
  ylab(label = "Observed events") +
  facet_wrap(strain ~ .) +
  coord_flip() +
  scale_fill_viridis(discrete = T,
                     name = "Family:") +
  theme(legend.position = "bottom",
        strip.text = element_markdown(),
        axis.text.y = element_markdown(),
        legend.text = element_markdown())

# deletions
isCountPerLib$del = ggplot(dfdel,
                           (aes(ISName, fill = ISFamily))) +
  geom_bar() +
  ylim(c(0, 3)) +
  facet_wrap(strain ~ .) +
  coord_flip() +
  ylab("Observed events") +
  xlab(label = "") +
  scale_fill_viridis(discrete = T,
                     name = "Family:",
                     labels = c("IS*4*" = "IS*4*",
                                "IS*H3*" = "IS*H3*",
                                "Outras famílias" = "Other families")) +
  theme(legend.position = "bottom",
        strip.text = element_markdown(),
        axis.text.y = element_markdown(),
        legend.text = element_markdown())

isCountPerLibPanel = ggarrange(plotlist = isCountPerLib,
                               labels = "AUTO",
                               ncol = 1,
                               nrow = 2,
                               heights = c(1.5, 1))

ggsave(filename = "plots/mob_observed_events_panel.png",
       plot = isCountPerLibPanel,
       dpi = 600,
       width = 6,
       height = 7.25)

# merging dfins and dfdel to observe frequency status
dfins$svType = "insertion"
dfdel$svType = "excision"
dfinsdel = bind_rows(dfins, dfdel)

dfinsdel = dfinsdel %>%
  mutate(status = case_when(status == "predominant" ~ "Predominant",
                            status == "common" ~ "Common",
                            status == "rare" ~ "Rare",
                            TRUE ~ as.character(status)),
         svType = case_when(svType == "insertion" ~ "Insertion",
                            svType == "excision" ~ "Excision",
                            TRUE ~ as.character(status)))

lvs = mixedsort(levels(as.factor(dfinsdel$ISName)), decreasing = T)
dfinsdel$ISName = factor(dfinsdel$ISName, levels=lvs)
dfinsdel$status = factor(dfinsdel$status, levels=c("Predominant", "Common", "Rare"))
dfinsdel$svType = factor(dfinsdel$svType, levels=c("Insertion", "Excision"))

isCountPerStatus = ggplot(dfinsdel, (aes(ISName, fill=ISFamily))) +
  geom_bar() +
  facet_grid(status ~ svType, scales = "free_x") +
  ylab("Observed events") +
  coord_flip() +
  xlab(label="") +
  scale_fill_viridis(discrete = T,
                     name = "Family:",
                     labels = c("IS*4*" = "IS*4*",
                                "IS*H3*" = "IS*H3*",
                                "Outras famílias" = "Other families")) +
  theme(legend.position = "bottom",
        strip.text = element_markdown(),
        axis.text.y = element_markdown(),
        legend.text = element_markdown())

# creating a tibble containing the sum of insertions and deletions
# and normalizing it by read depth
# all cases
insDelCounts = dfinsdel %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(count = n()) %>% 
  dplyr::mutate(readCount = resCounts$counts) %>% 
  dplyr::mutate(norm = (max(readCount)/readCount)*count) %>% 
  dplyr::mutate(strain = str_replace(strain, " .$", "")) %>% 
  dplyr::group_by(strain) %>% 
  dplyr::summarise(mean = mean(norm),
                   sd = sd(norm),
                   n = n())

insDelCounts$margin = 0.9945*(insDelCounts$sd/sqrt(insDelCounts$n))
insDelCounts$lower95ci = insDelCounts$mean - insDelCounts$margin
insDelCounts$upper95ci = insDelCounts$mean + insDelCounts$margin

# plotting
comparisonPlot = insDelCounts %>% 
  ggplot(aes(x = strain, y = mean)) +
  geom_point() +
  geom_pointrange(aes(ymin = lower95ci,
                      ymax = upper95ci)) +
  coord_flip() +
  xlab("Strain") +
  ylab("Average number of observed mobilizations") +
  ggtitle("Total") +
  theme(axis.text.y.left = ggtext::element_markdown(),
        plot.title = ggtext::element_markdown())

# separated by family
fams = dfinsdel$ISFamily %>% 
  unique() %>% 
  as.character() 

isfamdf = list()
isfamplot = list()
isfamplot[["Todas"]] = comparisonPlot
for(i in fams){
  if(i != "Outras famílias"){plottitle = i}else{plottitle = "Other families"}
  
  isfamdf[[i]] = dfinsdel %>%
    dplyr::filter(ISFamily == i) %>% 
    dplyr::group_by(strain) %>%
    dplyr::summarise(count = n()) %>% 
    dplyr::mutate(readCount = resCounts$counts) %>% 
    dplyr::mutate(norm = (max(readCount)/readCount)*count) %>% 
    dplyr::mutate(strain = str_replace(strain, " .$", "")) %>% 
    dplyr::group_by(strain) %>% 
    dplyr::summarise(mean = mean(norm),
                     sd = sd(norm),
                     n = n())
  
  isfamdf[[i]]$margin = 0.9945*(isfamdf[[i]]$sd/sqrt(isfamdf[[i]]$n))
  isfamdf[[i]]$lower95ci = isfamdf[[i]]$mean - isfamdf[[i]]$margin
  isfamdf[[i]]$upper95ci = isfamdf[[i]]$mean + isfamdf[[i]]$margin
  
  isfamplot[[i]] = isfamdf[[i]] %>% 
    ggplot(aes(x = strain, y = mean)) +
    geom_point() +
    geom_pointrange(aes(ymin = lower95ci,
                        ymax = upper95ci)) +
    coord_flip() +
    xlab("Strain") +
    ylab("Average number of observed mobilizations") +
    ggtitle(plottitle) +
    theme(axis.text.x = ggtext::element_markdown(),
          axis.text.y = ggtext::element_markdown(),
          plot.title = ggtext::element_markdown())
}

ggsave(filename = "plots/mobilizationComparisonPerFamily_en.png",
       plot = ggarrange(plotlist = isfamplot,
                        nrow = isfamplot %>% length(),
                        ncol = 1,
                        labels = "AUTO"),
       dpi = 600,
       width = 7,
       height = 8)

# growth curve of SmAP1 knockout strain ####
# loading data
dsmap1_I = c(0.05, 0.21, 0.29, 0.41, 0.54, 0.60, 0.62, 0.60, 0.60, 0.61, 0.64, 0.70, 0.73, 0.74, 0.76, 0.82, 0.83, 0.83, 0.82, 0.85, 0.87, 0.85, 0.83, 0.93, 0.91)
dsmap1_II = c(0.05, 0.21, 0.28, 0.39, 0.53, 0.58, 0.63, 0.58, 0.58, 0.58, 0.60, 0.71, 0.69, 0.70, 0.72, 0.78, 0.79, 0.80, 0.79, 0.78, 0.79, 0.81, 0.83, 0.89, 0.85)
dsmap1_III = c(0.05, 0.22, 0.32, 0.42, 0.53, 0.57, 0.59, 0.57, 0.58, 0.57, 0.61, 0.70, 0.71, 0.73, 0.75, 0.81, 0.83, 0.83, 0.84, 0.82, 0.85, 0.86, 0.86, 0.94, 0.91)
dura3_I = c(0.05, 0.26, 0.38, 0.51, 0.59, 0.60, 0.58, 0.58, 0.60, 0.61, 0.64, 0.76, 0.80, 0.80, 0.82, 0.86, 0.92, 0.92, 0.93, 0.93, 0.94, 0.95, 0.95, 1.04, 1.02)
dura3_II = c(0.05, 0.25, 0.33, 0.45, 0.55, 0.53, 0.58, 0.58, 0.57, 0.58, 0.60, 0.74, 0.77, 0.78, 0.80, 0.83, 0.91, 0.92, 0.92, 0.93, 0.95, 0.96, 0.96, 1.05, 1.05)
dura3_III = c(0.05, 0.28, 0.39, 0.51, 0.60, 0.61, 0.60, 0.61, 0.63, 0.64, 0.64, 0.81, 0.83, 0.86, 0.88, 0.93, 0.98, 0.99, 1.01, 0.99, 1.03, 1.02, 1.01, 1.10, 1.07)

time = c(0,13,17,20,23,27,37,40,43,45,48,61,64,67,69,70,83,86,89,91,96,112,136,155,159)

# creating a tibble
growthCurveTib = tibble(od = c(dsmap1_I,dsmap1_II,dsmap1_III,dura3_I,dura3_II,dura3_III),
                        Strain = c(rep("&Delta;*ura3* &Delta;*smap1*", length(dsmap1_I) * 3),
                                   rep("&Delta;*ura3*", length(dsmap1_I) * 3)),
                        Replicate = rep(c(rep("BR1", length(dsmap1_I)),
                                   rep("BR2", length(dsmap1_I)),
                                   rep("BR3", length(dsmap1_I))),
                                 2),
                        time = rep(time, 6))

# plotting a time series
growthCurvePlot = growthCurveTib %>% 
  ggplot(aes(x = time,
             y = od,
             color = Strain,
             shape = Strain)) +
  geom_point() +
  geom_line(aes(linetype = Replicate)) +
  scale_color_tableau() +
  ylab("OD<sub>600nm</sub>") +
  xlab("Time (hours)") +
  ylim(c(0, 1.2)) +
  theme(legend.text = element_markdown(),
        axis.title.y = element_markdown())

# saving
ggsave(filename = "plots/smap1_ura3_growth_curve.png",
       plot = growthCurvePlot,
       dpi = 600,
       width = 6,
       height = 4)

# Whitehead et al. 2006 pot. post-transcription regulated ####
# getting genes from table  supplemental table 2
# only the unshaded that should mRNA should not correlate with protein
whitehead2006unshaded = list()
whitehead2006unshaded$mrnaup = c("VNG0150H","VNG0152G","VNG0177G","VNG0324G","VNG0330G","VNG0557H","VNG0592G","VNG0769H",
                                 "VNG0865C","VNG0867G","VNG0906H","VNG1084G","VNG1103G","VNG1133G","VNG1138G","VNG1143G",
                                 "VNG1255C","VNG1291H","VNG1442G","VNG1451C","VNG1511C","VNG1538H","VNG1698G","VNG1703G",
                                 "VNG1709G","VNG1715G","VNG1766C","VNG1775C","VNG1809H","VNG2006C","VNG2010G","VNG2043G",
                                 "VNG2047G","VNG2048G","VNG2115H","VNG2153G","VNG2170H","VNG2196G","VNG2469G","VNG2471G",
                                 "VNG2577C","VNG2595G","VNG2622H","VNG2648G","VNG6188H","VNG5068G","VNG5146H")

whitehead2006unshaded$mrnadown = c("VNG0330G",
                                   "VNG0574C","VNG0689G","VNG0743H","VNG0771G","VNG0775G","VNG0794G","VNG0804C","VNG0837H",
                                   "VNG0963G","VNG0995H","VNG1009G","VNG1021C","VNG1121G","VNG1213C","VNG1332G","VNG1408G",
                                   "VNG1467G","VNG1524C","VNG1536C","VNG1551G","VNG1558H","VNG1559H","VNG1568G","VNG1603G",
                                   "VNG1624G","VNG1667G","VNG1676G","VNG1898C","VNG2436G","VNG2627C","VNG5033G")

gammaradptgs = nrtxsep %>%
  filter(locus_tag %in% whitehead2006unshaded$mrnaup) %>%
  pull(representative) %>%
  unique()

gammaradptgsNewPtgs = gammaradptgs[!gammaradptgs %in% (c(ptgs$union$Q4$locus_tag, ptgsAbund$union$prot_bot_prot_non_mrna_top) %>% unique())]

# plotting all genes for manual inspection: abundance #####
# allLocusTags = abundLongFuncat$locus_tag %>%
#   sort() %>%
#   unique()

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
