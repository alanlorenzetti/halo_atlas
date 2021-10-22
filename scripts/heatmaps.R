# alorenzetti 20200625

# description ####
# this script will plot a heatmap
# based on fold change data
# for proteins, mRNAs and RPFs
# functional categorization
# will also be incorporated into
# the generated heatmaps

# loading libs #####
source("scripts/loadingLibs.R")

# disabling warnings about rstudio ide
# thrown by ComplexHeatmap
ht_opt$message = FALSE

# heatmaps
joinedTibble = list()
M = list()
M2 = list()
M3 = list()
htlist = list()
heatLegs = list()
row_ha = list()
p = list()
colors = colorRamp2(c(-5, 0, 5), c("#4E79A7", "white", "#E15759"))

# integrating costFunc_ProtmRNA_RPFmRNA_AllTP into heatmap
for(i in names(regRules)){
  
  # setting higher and lower time point contrasts
  ihigh = sub("^(.*)_vs_.*$","\\1",i)
  ilow = sub("^.*_vs_(.*$)","\\1",i)

  # creating and wrangling a dataframe based on the regRules dataset
  joinedTibble[[i]] = regRules[[i]][["Prot_mRNA"]] %>%
    bind_rows() %>% 
    distinct() %>% 
    left_join(., dictFunCat, by = c("locus_tag" = "pfeiLocusTag")) %>%
    mutate(lsmSense = case_when(is.na(lsmSense) ~ "no",
                                TRUE ~ as.character(lsmSense)),
           lsmAntiSense = case_when(is.na(lsmAntiSense) ~ "no",
                                    TRUE ~ as.character(lsmAntiSense)),
           arCOG = case_when(is.na(arCOG) ~ "Function unknown",
                             arCOG == "NA" ~ "Function unknown",
                             TRUE ~ as.character(arCOG)),
           utrSize = case_when(is.na(utrSize) ~ "Unknown",
                               TRUE ~ as.character(utrSize)),
           asRNA = case_when(is.na(asRNA) ~ "no",
                             TRUE ~ as.character(asRNA))
           )
  
  # creating the matrix of foldchanges
  M[[i]] = joinedTibble[[i]] %>%
    dplyr::select(matches("protein.*high") |
                    matches("mRNA.*low") |
                    matches("RPF.*low")) %>% 
    dplyr::select(-matches("quad"),
                  -matches("sig"),
                  -matches("lfcse")) %>% 
    as.matrix()
  
  rownames(M[[i]]) = paste0(joinedTibble[[i]]$locus_tag, "|", joinedTibble[[i]]$pfeiProduct) %>%
    unlist() %>%
    unname()
  colnames(M[[i]]) = c(paste0("Protein"," ", ihigh),
                       paste0("mRNA"," ", ilow),
                       paste0("RPF"," ", ilow))
  
  # creating the matrix of baseMean
  baseMean = tibble(locus_tag = sub("\\|.*$", "", unifiedFin$TP4$locus_tag),
                    totrnaBaseMean = unifiedFin$TP4$baseMean.x,
                    ribornaBaseMean = unifiedFin$TP4$baseMean.y)
  
  joinedTibble[[i]] = left_join(joinedTibble[[i]], baseMean, by = "locus_tag")
  
  M2[[i]] = joinedTibble[[i]] %>% 
    dplyr::select(totrnaBaseMean) %>% 
    as.matrix()
  
  rownames(M2[[i]]) = joinedTibble[[i]]$locus_tag %>% unlist() %>% unname()
  colnames(M2[[i]]) = "mRNA (Base Mean)"
  
  M3[[i]] = joinedTibble[[i]] %>% 
    dplyr::select(ribornaBaseMean) %>% 
    as.matrix()
  
  rownames(M3[[i]]) = joinedTibble[[i]]$locus_tag %>% unlist() %>% unname()
  colnames(M3[[i]]) = "RPF (Base Mean)"
  
  # setting colors for abundance matrices
  colors2 = colorRamp2(c(0, M2[[i]][,1] %>% quantile(0.95)), c("white", "#59A14F"))
  colors3 = colorRamp2(c(0, M3[[i]][,1] %>%  quantile(0.95)), c("white", "#59A14F"))
  
  #narcogs = joinedTibble[[i]]$arCOG %>% unique() %>% length(.)
  #arCOGcols = circlize::rand_color(n=narcogs, luminosity="bright")
  # setting up manually
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
  arCOGClasses = joinedTibble[[i]]$arCOG %>% sort() %>% unique()
  arCOGcols = arCOGcols[names(arCOGcols) %in% arCOGClasses]
  
  # defining colors for the heatmap and annots
  heatCols = list(
    lsmCol = c("no" = "white",
               "yes" = "#E15759"),
    regGroupCol =   c("ProteinUp_mRNAUp" = "#E15759",
                      "ProteinUp_mRNAZero" = "#FF9DA7",
                      "ProteinUp_mRNADown" = "#F28E2B",
                      "ProteinZero_mRNAUp" = "#B07AA1",
                      "ProteinZero_mRNAZero" = "#EDC948",
                      "ProteinZero_mRNADown" = "#76B7B2",
                      "ProteinDown_mRNAUp" = "#59A14F",
                      "ProteinDown_mRNAZero" = "#9C755F",
                      "ProteinDown_mRNADown" = "#4E79A7"),
    pmrmColFun = c("no" = "white",
                   "yes" = "#59A14F"),
    arCOGCol = arCOGcols,
    utrSizeCol = c("Unknown" = "grey70",
                   "leaderless" = "#4E79A7",
                   "short" = "#EDC948",
                   "mid" = "#F28E2B",
                   "long" = "#E15759"),
    asRNACol = c("no" = "white",
                 "yes" = "#B07AA1"),
    mfeCol = colorRamp2(breaks = c(0,joinedTibble[[i]]$utrmfe %>% min(na.rm = T) %>% floor()),
                        colors = c("white", "#4E79A7")),
    HLCol = colorRamp2(breaks = c(0,joinedTibble[[i]]$HL %>% max(na.rm = T) %>% ceiling()),
                       colors = c("white", "#4E79A7")),
    caiCol = colorRamp2(breaks = c(0.5,joinedTibble[[i]]$cai %>% max(na.rm = T)),
                        colors = c("white", "#4E79A7"))
  )
  
  # defining colors and values for legends
  heatLegs[[i]] = list(
    regGroup =  Legend(title = "Regulation Group",
                       at = heatCols$regGroupCol %>% names() %>% sub("_"," ", .),
                       legend_gp  = gpar(fill = heatCols$regGroupCol %>% unname()),
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
    #    pmrm = Legend(title = "PMRM",
    #                  at = c("No","Yes"),
    #                  legend_gp = gpar(fill = heatCols$pmrmColFun %>% unname()),
    #                  border="black"),
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
    arCOG = Legend(title = "arCOG",
                   at = names(heatCols$arCOGCol),
                   legend_gp  = gpar(fill = heatCols$arCOGCol %>% unname()),
                   border="black")
  )
  
  # defining annotations
  row_ha[[i]] = HeatmapAnnotation(which = "row",
                                  regGroup = anno_simple(joinedTibble[[i]]$regRule,
                                                         col = heatCols$regGroupCol,
                                                         border = T),
#                                 operon = anno_simple(joinedTibble[[i]]$operon,
#                                                       border = T),
                                  arCOG = anno_simple(joinedTibble[[i]]$arCOG,
                                                      border = T,
                                                      col = heatCols$arCOGCol),
                                  lsmSense = anno_simple(joinedTibble[[i]]$lsmSense,
                                                         col = heatCols$lsmCol,
                                                         border = T),
                                  lsmAntiSense = anno_simple(joinedTibble[[i]]$lsmAntiSense,
                                                             col = heatCols$lsmCol,
                                                             border = T),
                                  asRNA = anno_simple(joinedTibble[[i]]$asRNA,
                                                      col = heatCols$asRNACol,
                                                      border = T),
#                                  pmrm = anno_simple(joinedTibble[[i]]$pmrm,
#                                                     col = heatCols$pmrmCol,
#                                                     border = T),
                                  utrSize = anno_simple(joinedTibble[[i]]$utrSize,
                                                        col = heatCols$utrSizeCol,
                                                        border = T),
                                  mfe = anno_simple(joinedTibble[[i]]$utrmfe,
                                                    col = heatCols$mfeCol,
                                                    border = T),
                                  HL = anno_simple(joinedTibble[[i]]$HL,
                                                   col = heatCols$HLCol,
                                                   border = T),
                                  cai = anno_simple(joinedTibble[[i]]$cai,
                                                    col = heatCols$caiCol,
                                                    border = T),
                                  annotation_label = c("Regulation Group",
#                                                       "Operon",
                                                       "arCOG",
                                                       "LSm Sense",
                                                       "LSm AntiSense",
#                                                       "PMRM",
                                                       "asRNA",
                                                       "5' UTR Length",
                                                       "5' UTR MFE",
                                                       "Half life",
                                                       "CAI"))
  
  # creating heatmaps
  p[[i]][["lfc"]] = Heatmap(M[[i]],
                            name = "LFC",
                            col = colors,
                            show_heatmap_legend = T,
                            row_names_side = "right",
                            show_row_names = T,
                            column_order = colnames(M[[i]]),
                            left_annotation = row_ha[[i]],
                            row_split = factor(joinedTibble[[i]]$regRule,
                                               levels=joinedTibble[[i]]$regRule %>% unique()),
#                                              row_split = factor(joinedTibble[[i]]$utrSize),
#                                              row_split = factor(joinedTibble[[i]]$lsmSense),
#                                              row_split = factor(joinedTibble[[i]]$arCOG),
#                                              row_split = factor(joinedTibble[[i]]$asRNA),
                            cluster_row_slices = FALSE,
                            row_title = NULL,
                            border = T,
                            heatmap_legend_param = list(
                              title = "LFC",
                              at = c(-5, 0, 5),
                              border = T))
  
  p[[i]][["totrnaAbundance"]] = Heatmap(M2[[i]],
                                  name = "mRNA",
                                  col = colors2,
                                  show_heatmap_legend = T,
                                  row_names_side = "right",
                                  show_row_names = F,
                                  row_names_gp = gpar(fontsize = 5),
                                  column_order = colnames(M2[[i]]),
                                  cluster_row_slices = FALSE,
                                  row_title = NULL,
                                  border = T,
                                  heatmap_legend_param = list(
                                    title = "mRNA",
                                    border = T))
  
  p[[i]][["riboAbundance"]] = Heatmap(M3[[i]],
                                  name = "RPF",
                                  col = colors3,
                                  show_heatmap_legend = T,
                                  row_names_side = "right",
                                  show_row_names = F,
                                  row_names_gp = gpar(fontsize = 5),
                                  column_order = colnames(M3[[i]]),
                                  cluster_row_slices = FALSE,
                                  row_title = NULL,
                                  border = T,
                                  heatmap_legend_param = list(
                                    title = "RPF",
                                    border = T))
  
  # creating list of heatmaps
  htlist[[i]] = p[[i]][["lfc"]] + p[[i]][["totrnaAbundance"]] + p[[i]][["riboAbundance"]]
  
  # drawing heatmaps using
  # complex heatmap shiny app
  draw(htlist[[i]], annotation_legend_list = heatLegs[[i]])
}

#ht_shiny(htlist[[i]])