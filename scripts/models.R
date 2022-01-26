# alorenze 20200525

# description ####
# this script will generate
# linear models to explain
# proteins in function of mRNA and/or RPFs
# it will also plot several charts

# setting up thresholds and vars ####
log2fcthreshold = 1
borderlinezero = 0.25

# defining tab10 color scheme
# ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`
tab10col = c(
  "Q1"="#E15759", # red
  "Q2"="#F28E2B", # orange
  "Q3"="#4E79A7", # blue
  "Q4"="#59A14F", # green
  "BL12"="#FF9DA7", # pink
  "BL41"="#B07AA1", # purple
  "BL23"="#76B7B2", # light teal 
  "BL34"="#9C755F", # brown
  "core"="#EDC948", # yellow
  "center"="grey80") # grey

tab10labs = c("Protein Up & mRNA Up",
              "Protein Up & mRNA Down",
              "Protein Down & mRNA Down",
              "Protein Down & mRNA Up",
              "Protein Up & mRNA Flat",
              "Protein Flat & mRNA Up",
              "Protein Flat & mRNA Down",
              "Protein Down & mRNA Flat",
              "Protein Flat & mRNA Flat",
              "None of the previous categories")

names(tab10labs) = names(tab10col)

# mod to allow black lines using fill aes
tab10colmod = tab10col
tab10colmod["Q1"] = "black"
tab10colmod["Q3"] = "black"

# setting ggplot blank theme
theme_set(theme_bw())

# unifying lfc datasets ####
alldfs = list()
quantileNorm = list()
for(i in names(unifiedFin)){
  
  tphigh = str_replace(i, "(.*)_vs_.*", "\\1")
  tplow = str_replace(i, ".*_vs_(.*)", "\\1")
    
  varName[1] = paste0("sig_protein")
  varName[2] = paste0("borderlineZero_protein")
  varName[3] = paste0("sig_totrna")
  varName[4] = paste0("sig_ribo")
  varName[5] = paste0("borderlineZero_totrna")
  varName[6] = paste0("borderlineZero_ribo")
  
  alldfs[[i]] = inner_join(unifiedFinDEProt[[i]],
                           unifiedFin[[i]],
                           by = "locus_tag") %>% 
    dplyr::select(locus_tag,
                  matches("lfc"),
                  matches("sig"),
                  matches("border")) %>%
    rename("mRNA" = mean_lfc_rna_total,
           "RPF" = mean_lfc_rna_ribo,
           "protein" = mean_lfc_protein_lysate,
           "lfcse_mRNA" = se_lfc_rna_total,
           "lfcse_RPF" = se_lfc_rna_ribo,
           "lfcse_protein" = se_lfc_protein_lysate)
  
  # # I will apply quantile normalization to mRNA, RPF and protein
  # # and then I will replace the original values by the normalized ones
  # # the normalized values correlate quite well with the original ones
  # # so I don't see a reason to be concerned about adjusting significance
  # # normalizing lfc values
  # quantileNorm[[i]] = alldfs[[i]] %>%
  #   dplyr::select(mRNA, RPF, protein) %>%
  #   as.matrix() %>%
  #   normalize.quantiles() %>%
  #   as_tibble()
  # colnames(quantileNorm[[i]]) = c("mRNA", "RPF", "protein")
  # 
  # alldfs[[i]][,c("mRNA","RPF","protein")] = quantileNorm[[i]][,c("mRNA","RPF","protein")]
  # 
  # # normalizing standard errors
  # quantileNorm[[i]] = alldfs[[i]] %>%
  #   dplyr::select(lfcse_mRNA, lfcse_RPF, lfcse_protein) %>%
  #   as.matrix() %>%
  #   normalize.quantiles() %>%
  #   as_tibble()
  # colnames(quantileNorm[[i]]) = c("lfcse_mRNA", "lfcse_RPF", "lfcse_protein")
  # 
  # alldfs[[i]][,c("lfcse_mRNA","lfcse_RPF","lfcse_protein")] = quantileNorm[[i]][,c("lfcse_mRNA","lfcse_RPF","lfcse_protein")]
  
  # finding significant and borderline status
  for(type in c("totrna", "ribo", "protein")){
    var1 = paste0("sig_",type)
    var2 = paste0("borderlineZero_",type)
    
    if(type == "totrna"){var3 = "mRNA"}
    if(type == "ribo"){var3 = "RPF"}
    if(type == "protein"){var3 = "protein"}
    
    alldfs[[i]] = alldfs[[i]] %>% 
      mutate(!!var1 := case_when(abs(get(var3)) >= log2fcthreshold & get(var1) == "yes" ~ "yes",
                                 TRUE ~ "no")) %>% 
      mutate(!!var2 := case_when(abs(get(var3)) < borderlinezero ~ "yes",
                                 TRUE ~ "no"))
  }
  
  # creating new categorical variables
  alldfs[[i]] = alldfs[[i]] %>% 
    mutate(`RPF-mRNA-quad` = case_when(mRNA >= log2fcthreshold & RPF >= log2fcthreshold ~ "Q1",
                                       mRNA <= -log2fcthreshold & RPF >= log2fcthreshold ~ "Q2",
                                       mRNA <= -log2fcthreshold & RPF <= -log2fcthreshold ~ "Q3",
                                       mRNA >= log2fcthreshold & RPF <= -log2fcthreshold ~ "Q4",
                                       mRNA >= borderlinezero & abs(RPF) < borderlinezero ~ "BL41",
                                       mRNA <= -borderlinezero & abs(RPF) < borderlinezero ~ "BL23",
                                       RPF >= borderlinezero & abs(mRNA) < borderlinezero ~ "BL12",
                                       RPF <= -borderlinezero & abs(mRNA) < borderlinezero ~ "BL34",
                                       TRUE ~ "center")) %>% 
    
    mutate(`protein-mRNA-quad` = case_when(mRNA >= log2fcthreshold & protein >= log2fcthreshold ~ "Q1",
                                           mRNA <= -log2fcthreshold & protein >= log2fcthreshold ~ "Q2",
                                           mRNA <= -log2fcthreshold & protein <= -log2fcthreshold ~ "Q3",
                                           mRNA >= log2fcthreshold & protein <= -log2fcthreshold ~ "Q4",
                                           mRNA >= borderlinezero & abs(protein) < borderlinezero ~ "BL41",
                                           mRNA <= -borderlinezero & abs(protein) < borderlinezero ~ "BL23",
                                           protein >= borderlinezero & abs(mRNA) < borderlinezero ~ "BL12",
                                           protein <= -borderlinezero & abs(mRNA) < borderlinezero ~ "BL34",
                                           TRUE ~ "center")) %>%
    
    mutate(`protein-RPF-quad` = case_when(RPF >= log2fcthreshold & protein >= log2fcthreshold ~ "Q1",
                                          RPF <= -log2fcthreshold & protein >= log2fcthreshold ~ "Q2",
                                          RPF <= -log2fcthreshold & protein <= -log2fcthreshold ~ "Q3",
                                          RPF >= log2fcthreshold & protein <= -log2fcthreshold  ~ "Q4",
                                          RPF >= borderlinezero & abs(protein) < borderlinezero ~ "BL41",
                                          RPF <= -borderlinezero & abs(protein) < borderlinezero ~ "BL23",
                                          protein >= borderlinezero & abs(RPF) < borderlinezero ~ "BL12",
                                          protein <= -borderlinezero & abs(RPF) < borderlinezero ~ "BL34",
                                          TRUE ~ "center")) %>% 
    
    mutate(`RPF-mRNA-quad` = case_when(get(varName[4]) == "yes" & get(varName[3]) == "yes" ~ as.character(`RPF-mRNA-quad`),
                                       get(varName[4]) == "yes" & get(varName[3]) == "no" ~ as.character(`RPF-mRNA-quad`),
                                       get(varName[4]) == "no" & get(varName[3]) == "yes" ~ as.character(`RPF-mRNA-quad`),
                                       abs(RPF) < borderlinezero & abs(mRNA) < borderlinezero ~ "core",
                                       TRUE ~ "center"),
           `protein-mRNA-quad` = case_when(get(varName[1]) == "yes" & get(varName[3]) == "yes" ~ as.character(`protein-mRNA-quad`),
                                           get(varName[1]) == "yes" & get(varName[3]) == "no" ~ as.character(`protein-mRNA-quad`),
                                           get(varName[1]) == "no" & get(varName[3]) == "yes" ~ as.character(`protein-mRNA-quad`),
                                           abs(protein) < borderlinezero & abs(mRNA) < borderlinezero ~ "core",
                                           TRUE ~ "center"),
           `protein-RPF-quad` = case_when(get(varName[1]) == "yes" & get(varName[4]) == "yes" ~ as.character(`protein-RPF-quad`),
                                          get(varName[1]) == "yes" & get(varName[4]) == "no" ~ as.character(`protein-RPF-quad`),
                                          get(varName[1]) == "no" & get(varName[4]) == "yes" ~ as.character(`protein-RPF-quad`),
                                          abs(protein) < borderlinezero & abs(RPF) < borderlinezero ~ "core",
                                          TRUE ~ "center")) %>% 
    
    mutate(RPF_mRNA_sig = case_when(get(varName[4]) == "yes" & get(varName[3]) == "yes" ~ "bold",
                                    get(varName[4]) == "yes" & get(varName[3]) == "no" ~ "mid",
                                    get(varName[4]) == "no" & get(varName[3]) == "yes" ~ "mid",
                                    TRUE ~ "faint"),
           protein_mRNA_sig = case_when(get(varName[1]) == "yes" & get(varName[3]) == "yes" ~ "bold",
                                        get(varName[1]) == "yes" & get(varName[3]) == "no" ~ "mid",
                                        get(varName[1]) == "no" & get(varName[3]) == "yes" ~ "mid",
                                        TRUE ~ "faint"),
           protein_RPF_sig = case_when(get(varName[1]) == "yes" & get(varName[4]) == "yes" ~ "bold",
                                       get(varName[1]) == "yes" & get(varName[4]) == "no" ~ "mid",
                                       get(varName[1]) == "no" & get(varName[4]) == "yes" ~ "mid",
                                       TRUE ~ "faint")) %>% 
    
    mutate(contrast = i)
}

# functions to find axes lims
findXLims = function(df){
  mrnamax = df[,"mRNA"] %>% abs() %>% max() %>% RoundTo(x = ., multiple = 2, FUN = ceiling)
  rpfmax = df[,"RPF"] %>% abs() %>% max() %>% RoundTo(x = ., multiple = 2, FUN = ceiling)
  xmax = max(mrnamax,rpfmax)
  xmin = -xmax
  Xlims = c(xmin,xmax)
  return(Xlims)
}
findYLims = function(df){
  rpfmax = df[,"RPF"] %>% abs() %>% max() %>% RoundTo(x = ., multiple = 2, FUN = ceiling)
  protmax = df[,"protein"] %>% abs() %>% max() %>% RoundTo(x = ., multiple = 2, FUN = ceiling)
  ymax = max(rpfmax,protmax)
  ymin = -ymax
  Ylims = c(ymin,ymax)
  return(Ylims)
}

# getting models to explain relationship
# between y and x
modelGeneral = function(df, type){
  if(type == "RPF-mRNA"){
    y = "RPF"; x = "mRNA"}
  if(type == "protein-mRNA"){
    y = "protein"; x = "mRNA"}
  if(type == "protein-RPF"){
    y = "protein"; x = "RPF"}
  md = lm(formula=paste0(y,"~",x),data=df)
  return(md)
}

# running modelling
models = list()
types = c("RPF-mRNA", "protein-mRNA", "protein-RPF")
for(i in names(alldfs)){
  for(k in types){
    models[[i]][[k]] = modelGeneral(alldfs[[i]], k) %>% summary()
  }
}

# plotting general trends
plotGeneral = function(df, type, xLim, yLim){
  if(type == "RPF-mRNA"){
    y = "RPF"; x = "mRNA"; form = paste0(y,"~",x)}
  if(type == "protein-mRNA"){
    y = "protein"; x = "mRNA"; form = paste0(y,"~",x)}
  if(type == "protein-RPF"){
    y = "protein"; x = "RPF"; form = paste0(y,"~",x)}
  
  model = lm(formula = form, data = df)
  cortext = paste0("*R<sup>2</sup>* = ", summary(model)$r.squared %>% round(digits = 3))
  p = pf(summary(model)$fstatistic[1],
         summary(model)$fstatistic[2],
         summary(model)$fstatistic[3],
         lower.tail = FALSE) %>% 
    formatC(format = "e", digits = 2)
  pvaltext = paste0("*p* = ", p)
  
  grob = grobTree(richtext_grob(text = paste0(cortext, "; ", pvaltext),
                           gp=gpar(fontsize=8),
                           x=0.075, y=0.95, hjust=0))
  
  alp = 0.05
  sz = 2
  
  plot = ggplot(data = df) +
    geom_point(aes(x=get(x), y=get(y)), alpha = alp, stroke = 0) +
    xlim(xLim) + ylim(yLim) + xlab(x) + ylab(y) +
#    geom_vline(xintercept = 0, size = 0.3) +
#    geom_hline(yintercept = 0, size = 0.3) +
    geom_point(inherit.aes = F,
                data = df,
                aes(x=get(x), y=get(y)),
                alpha = alp,
                show.legend = F,
                stroke = 0) #+
    # geom_smooth(inherit.aes = F,
    #             data = df,
    #             aes(x=get(x), y=get(y)),
    #             method="lm",
    #             alpha = alp,
    #             color = "blue",
    #             show.legend = F,
    #             size = 0.5) +
    # annotation_custom(grob = grob)
  
  return(plot)
}

# generating general plots
genplots = list()
for(i in names(alldfs)){
  xLim = findXLims(alldfs[[i]])
  yLim = findYLims(alldfs[[i]])
  for(k in types){
    genplots[[i]][[k]] = plotGeneral(alldfs[[i]], k, xLim, yLim)}
}

# arranging plots (full)
ggarrange(plotlist = c(genplots$TP2_vs_TP1,
                       genplots$TP3_vs_TP1,
                       genplots$TP4_vs_TP1,
                       genplots$TP3_vs_TP2,
                       genplots$TP4_vs_TP3),
          nrow = 5, ncol = 3)

# plot functions with sig dataset
plotLM = function(df, type, xLim, yLim){
  if(type == "RPF-mRNA"){
    y = "RPF"; x = "mRNA"; form = paste0(y,"~",x); color = paste0(type, "-quad"); alpha = paste0(y,"_",x,"_","sig")}
  if(type == "protein-mRNA"){
    y = "protein"; x = "mRNA"; form = paste0(y,"~",x); color = paste0(type, "-quad"); alpha = paste0(y,"_",x,"_","sig")}
  if(type == "protein-RPF"){
    y = "protein"; x = "RPF"; form = paste0(y,"~",x); color = paste0(type, "-quad"); alpha = paste0(y,"_",x,"_","sig")}

  if(y == "RPF"){var1 = "sig_ribo"; ylabtxt = "Log<sub>2</sub> Fold Change RPF"}
  if(y == "protein"){var1 = "sig_protein"; ylabtxt = "Log<sub>2</sub> Fold Change Protein"}
  if(x == "mRNA"){var2 = "sig_totrna"; xlabtxt = "Log<sub>2</sub> Fold Change mRNA"}
  if(x == "RPF"){var2 = "sig_ribo"; xlabtxt = "Log<sub>2</sub> Fold Change RPF"}
  
  varName1 = df %>%
    dplyr::select(starts_with(var1)) %>% names() %>% as.symbol()
  varName2 = df %>%
    dplyr::select(starts_with(var2)) %>% names() %>% as.symbol()
  colsymb = color %>% as.symbol()
  dfSig = df %>% filter((!!varName1 == "yes" & !!varName2 == "yes") & (!!colsymb == "Q1" | !!colsymb == "Q3"))
#  dfNonSig =  df %>% filter(!!varName2 == "no")
  
#  alp = 0.25
  sz = 1.5
  
  var = color %>% as.symbol()
  model = list() ; cortext = list()
  p = list() ; pvaltext = list()
  grob = list() ; grobY = 0.92
  
  quadsItSummy = dfSig %>% dplyr::select(!!var) %>% table()
  quadsIt = names(quadsItSummy[quadsItSummy > 1])
  quadsIt = quadsIt[quadsIt %in% paste0("Q", c(1,3))]
  
  for(i in quadsIt){
    model[[i]] = lm(formula = form, data = dfSig %>% filter(!!var == i))
    cortext[[i]] = paste0("*R<sup>2</sup>* = ", summary(model[[i]])$r.squared %>% round(digits = 3))
    p[[i]] = pf(summary(model[[i]])$fstatistic[1],
                summary(model[[i]])$fstatistic[2],
                summary(model[[i]])$fstatistic[3],
                lower.tail = FALSE) %>% 
      formatC(format = "e", digits = 2)
    pvaltext[[i]] = paste0("*p* = ", p[[i]])
    if(i == "Q1"){col = tab10col["Q1"]}
#    if(i == "Q2"){col = tab10col["Q2"]}
    if(i == "Q3"){col = tab10col["Q3"]}
#    if(i == "Q4"){col = tab10col["Q4"]}
    grob[[i]] = grobTree(richtext_grob(text = paste0(cortext[[i]], "; ", pvaltext[[i]]),
                                  gp=gpar(fontsize=8, col=col),
                                  x=0.075, y=grobY, hjust=0))
    grobY = grobY - 0.05
  }
  
  plot = ggplot(data = df, aes(x=get(x), y=get(y), colour=get(color))) +
    geom_point(aes(alpha=get(alpha)), size = sz, show.legend = T, stroke = 0) +
    xlim(xLim) + ylim(yLim) + xlab(xlabtxt) + ylab(ylabtxt) +
    # geom_vline(xintercept = 0, size = 0.3) +
    # geom_hline(yintercept = 0, size = 0.3) +
    # geom_point(inherit.aes = F,
    #             data = dfSig,
    #             aes(x=get(x), y=get(y), color=get(color)),
    #             alpha = alp, size = sz,
    #             show.legend = F,
    #             stroke = 0) +
    # geom_smooth(inherit.aes = F,
    #             data = dfSig,
    #             aes(x=get(x), y=get(y), fill=get(color)),
    #             colour="black",
    #             size=0.5,
    #             method="lm",
    #             show.legend = F) +
    scale_colour_manual(name = NULL,
                        values = tab10col,
                        labels = tab10labs,
                        guide = guide_legend(override.aes = list(size = 3))) +
    scale_alpha_manual(guide = "none",
                       values = c("bold"=0.8, "mid"= 0.3, "faint"=0.1)) +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown())
#    scale_shape_manual(values = c("bold"=4, "mid"=5, "faint"=1))
  
  # for(i in quadsIt){
  #   plot = plot +
  #     annotation_custom(grob = grob[[i]])
  # }
  
  return(plot)
}

# types of plot we should get
types = c("RPF-mRNA", "protein-mRNA", "protein-RPF")

######## linear models and plots for full dataset
lmplots = list()
for(i in names(alldfs)){
  xLim = findXLims(alldfs[[i]])
  yLim = findYLims(alldfs[[i]])
  for(k in types){
    # with computed lims
    #lmplots[[i]][[k]] = plotLM(alldfs[[i]], k, xLim, yLim)
    
    # with preset lims
    lmplots[[i]][[k]] = plotLM(alldfs[[i]], k, c(-10,10), c(-10,10)) +
      scale_x_continuous(limits = c(-10,10), breaks = seq(-10,10,2)) +
      scale_y_continuous(limits = c(-10,10), breaks = seq(-10,10,2))
  }
}

# arranging plots (full panel)
# ggarrange(plotlist = c(lmplots$TP2_vs_TP1,
#                        lmplots$TP3_vs_TP1,
#                        lmplots$TP4_vs_TP1,
#                        lmplots$TP3_vs_TP2,
#                        lmplots$TP4_vs_TP3),
#           nrow = 5, ncol = 3)

# arranging plots (only protein vs mRNA)
#creating list of plots
p_prot_mRNA = list("TP2_vs_TP1" = lmplots$TP2_vs_TP1$`protein-mRNA` + ggtitle("TP2_vs_TP1") + theme(legend.position="none"),
                   "TP3_vs_TP2" = lmplots$TP3_vs_TP2$`protein-mRNA` + ggtitle("TP3_vs_TP2") + theme(legend.position="none"),
                   "TP4_vs_TP3" = lmplots$TP4_vs_TP3$`protein-mRNA` + ggtitle("TP4_vs_TP3") + theme(legend.position="none"),
                   "TP3_vs_TP1" = lmplots$TP3_vs_TP1$`protein-mRNA` + ggtitle("TP3_vs_TP1") + theme(legend.position="none"),
                   "TP4_vs_TP1" = lmplots$TP4_vs_TP1$`protein-mRNA` + ggtitle("TP4_vs_TP1") + theme(legend.position="none"),
                   "LEGEND" = get_legend(lmplots$TP4_vs_TP1$`protein-mRNA`, position = "bottom"))

prot_mRNA_panel = ggarrange(plotlist = p_prot_mRNA,
                            nrow = 2, ncol = 3,
                            labels = c(LETTERS[1:5]))

ggsave(filename = "plots/ptgs_lfc_panel.png",
       plot = prot_mRNA_panel,
       width = 8,
       height = 5.5,
       units = "in",
       dpi = 300)

# plotting elements of the panel
# individually
for(i in names(lmplots)){
  ggsave(filename = paste0("plots/ptgs", i, ".png"),
         plot = lmplots[[i]]$`protein-mRNA` +
           ggtitle(gsub(i, pattern = "_vs_", replacement = " vs. ")) +
           theme(legend.position = "none"),
         width = 2.5,
         height = 2.5,
         units = "in",
         dpi = 300)
}

# getting properties of potentially
# post-transcriptionally regulated genes
ptgs = list()
for(i in names(alldfs)){
  # Q1 (red)
  ptgs[[i]][["Q1"]]$locus_tag = alldfs[[i]] %>%
    filter(sig_protein == "yes" & sig_totrna == "yes" & `protein-mRNA-quad` == "Q1") %>%
    pull(locus_tag)
  
  ptgs[[i]][["Q1"]]$tbl = dictFunCat %>% 
    filter(pfeiLocusTag %in% ptgs[[i]][["Q1"]]$locus_tag) %>% 
    select(-matches("ChIP"))
  
  # Q2 (orange)
  ptgs[[i]][["Q2"]]$locus_tag = alldfs[[i]] %>%
    filter(sig_protein == "yes" & sig_totrna == "yes" & `protein-mRNA-quad` == "Q2") %>%
    pull(locus_tag)
  
  ptgs[[i]][["Q2"]]$tbl = dictFunCat %>% 
    filter(pfeiLocusTag %in% ptgs[[i]][["Q2"]]$locus_tag) %>% 
    select(-matches("ChIP"))
  
  # Q3 (red)
  ptgs[[i]][["Q3"]]$locus_tag = alldfs[[i]] %>%
    filter(sig_protein == "yes" & sig_totrna == "yes" & `protein-mRNA-quad` == "Q3") %>%
    pull(locus_tag)
  
  ptgs[[i]][["Q3"]]$tbl = dictFunCat %>% 
    filter(pfeiLocusTag %in% ptgs[[i]][["Q3"]]$locus_tag) %>% 
    select(-matches("ChIP"))
  
  # Q4 (green)
  ptgs[[i]][["Q4"]]$locus_tag = alldfs[[i]] %>%
    filter(sig_protein == "yes" & sig_totrna == "yes" & `protein-mRNA-quad` == "Q4") %>%
    pull(locus_tag)
  
  ptgs[[i]][["Q4"]]$tbl = dictFunCat %>% 
    filter(pfeiLocusTag %in% ptgs[[i]][["Q4"]]$locus_tag) %>% 
    select(-matches("ChIP"))
  
  # BL12 (pink)
  ptgs[[i]][["BL12"]]$locus_tag = alldfs[[i]] %>%
    filter(sig_protein == "yes" & `protein-mRNA-quad` == "BL12") %>%
    pull(locus_tag)
  
  ptgs[[i]][["BL12"]]$tbl = dictFunCat %>% 
    filter(pfeiLocusTag %in% ptgs[[i]][["BL12"]]$locus_tag) %>% 
    select(-matches("ChIP"))
  
  # BL23 (light teal)
  ptgs[[i]][["BL23"]]$locus_tag = alldfs[[i]] %>%
    filter(sig_totrna == "yes" & `protein-mRNA-quad` == "BL23") %>%
    pull(locus_tag)
  
  ptgs[[i]][["BL23"]]$tbl = dictFunCat %>% 
    filter(pfeiLocusTag %in% ptgs[[i]][["BL23"]]$locus_tag) %>% 
    select(-matches("ChIP"))
  
  # BL34 (brown)
  ptgs[[i]][["BL34"]]$locus_tag = alldfs[[i]] %>%
    filter(sig_protein == "yes" & `protein-mRNA-quad` == "BL34") %>%
    pull(locus_tag)
  
  ptgs[[i]][["BL34"]]$tbl = dictFunCat %>% 
    filter(pfeiLocusTag %in% ptgs[[i]][["BL34"]]$locus_tag) %>% 
    select(-matches("ChIP"))
  
  # BL41 (purple)
  ptgs[[i]][["BL41"]]$locus_tag = alldfs[[i]] %>%
    filter(sig_totrna == "yes" & `protein-mRNA-quad` == "BL41") %>%
    pull(locus_tag)
  
  ptgs[[i]][["BL41"]]$tbl = dictFunCat %>% 
    filter(pfeiLocusTag %in% ptgs[[i]][["BL41"]]$locus_tag) %>% 
    select(-matches("ChIP"))
}

# getting union of ptgs locus tags
for(i in names(ptgs$TP2_vs_TP1)){
  ptgs$union[[i]][["locus_tag"]] = c(ptgs$TP2_vs_TP1[[i]]$locus_tag,
                                     ptgs$TP3_vs_TP2[[i]]$locus_tag,
                                     ptgs$TP4_vs_TP3[[i]]$locus_tag,
                                     ptgs$TP3_vs_TP1[[i]]$locus_tag,
                                     ptgs$TP4_vs_TP1[[i]]$locus_tag) %>% 
    unique()
}
