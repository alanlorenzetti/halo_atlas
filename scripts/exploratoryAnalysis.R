# alorenzetti 20210119

# description ####
# this script will perform PCA
# and plot the results
# using protein raw data
# rna-seq tpms and
# ribo-seq tpms

# setting up ####
# creating dir to store plots 
if(!dir.exists("plots/exploratory_analysis_v3")){dir.create("plots/exploratory_analysis_v3")}

# exploring proteome dataset ####
# generating a heat map to check
# sample grouping
# defining heatmap col funct
colfunct = circlize::colorRamp2(breaks = c(0,3,6), colors = viridis(3))

# creating matrix
M = spectro[ ,-1] %>%
  select(-contains("ribo")) %>% 
  drop_na() %>%
  as.matrix() %>%
  log10()

# defining annotation colors
legendCols = list(
  timepoints = c("TP1" = "#E15759",
                 "TP2" = "#F28E2B",
                 "TP3" = "#4E79A7",
                 "TP4" = "#59A14F"),
  bioReps = c("BR1" = adjustcolor(col = "#76B7B2", alpha.f = 0.2),
              "BR2" = adjustcolor(col = "#76B7B2", alpha.f = 0.5),
              "BR3" = adjustcolor(col = "#76B7B2", alpha.f = 1)),
  runs = c("r01" = adjustcolor(col = "#EDC948", alpha.f = 0.2),
           "r02" = adjustcolor(col = "#EDC948", alpha.f = 0.5),
           "r03" = adjustcolor(col = "#EDC948", alpha.f = 1))
)

# defining annotations
tpannots = colnames(M) %>% str_replace(".*(TP[1-4]).*", "\\1")
brannots = colnames(M) %>% str_replace(".*(BR[1-3]).*", "\\1")
runannots = colnames(M) %>% str_replace(".*(r0[1-3])$", "\\1")

annot = HeatmapAnnotation(which = "col",
                          timepoints = anno_simple(tpannots,
                                                   col = legendCols$timepoints,
                                                   border = T),
                          bioReps = anno_simple(brannots,
                                                col = legendCols$bioReps,
                                                border = T),
                          runs = anno_simple(runannots,
                                             col = legendCols$runs,
                                             border = T),
                          annotation_label = c("Time Point",
                                               "Bio. Replicate",
                                               "Run")
                          # annotation_label = c("Ponto de Coleta",
                          #                      "Réplica Biológica",
                          #                      "Corrida")
)

# defining legends
legs = list(
  timepoints = Legend(title = "Time Point",
#  timepoints = Legend(title = "Ponto de Coleta",
                      at = legendCols$timepoints %>% names(),
                      legend_gp = gpar(fill = legendCols$timepoints %>% unname()),
                      border = "black"),
  bioReps = Legend(title = "Bio. Replicate",
#  bioReps = Legend(title = "Réplica Biológica",
                   at = legendCols$bioReps %>% names(),
                   legend_gp = gpar(fill = legendCols$bioReps %>% unname()),
                   border = "black"),
  runs = Legend(title = "Run",
#  runs = Legend(title = "Corrida",
                at = legendCols$runs %>% names() %>% str_replace("r0", "R"),
                legend_gp = gpar(fill = legendCols$runs %>% unname()),
                border = "black")
)

# plotting heat map
ht = Heatmap(M, col = colfunct,
             top_annotation = annot,
             border = T,
             show_column_names = F,
             heatmap_legend_param = list(
               title = expression(Log[10](Abundância)),
               border = T,
               direction = "horizontal"))

#svglite("./plots/exploratory_analysis_v3/spectronautHeatmap.svg", width = 3, height = 4.5)
drawnht = grid.grabExpr(draw(ht,
               heatmap_legend_side = "bottom",
               annotation_legend_list = legs))
#dev.off()

# PCA
spectroPCA = spectro %>% 
  dplyr::select(-starts_with("ribo")) %>% 
  drop_na() %>% 
  transpose_tibble(locus_tag, id_col = "sample")

# adapted from:
# https://tbradley1013.github.io/2018/02/01/pca-in-a-tidy-verse-framework/
# creating a dataframe containing all pca info we need
dfpcaComplete = spectroPCA %>% 
  tidyr::nest(data=everything()) %>% 
  dplyr::mutate(pca = purrr::map(data, ~ prcomp(.x %>% select(-sample), 
                                  center = T, scale = T)),
         pca_aug = purrr::map2(pca, data, ~ broom::augment(.x, data = .y)))

# computing the explained variance for
# each of the principal components
expVar = dfpcaComplete %>% 
  unnest(pca_aug) %>% 
  summarize_at(.vars = vars(contains(".fittedPC")), .funs = list(~var(.))) %>% 
  gather(key = pc, value = variance) %>% 
  mutate(var_exp = variance/sum(variance),
         cum_var_exp = cumsum(var_exp),
         pc = str_replace(pc, ".fitted", ""))

# preparing df to plot PCA
pcaAugDf = dfpcaComplete %>%
  unnest(pca_aug) %>%
  mutate(tp = str_replace(sample,".*(TP[1-4]).*", "\\1"),
         br = str_replace(sample,".*(BR[1-3]).*", "\\1"),
         run = str_replace(sample,".*r0([1-3])$", "R\\1"))

# finding top three loadings in PC1
#dfpcaComplete$pca[[1]]$rotation[dfpcaComplete$pca[[1]]$rotation %>% .[,1] %>% order(decreasing = T),1] %>% head(3)

# finding bottom three loadings in PC1
#dfpcaComplete$pca[[1]]$rotation[dfpcaComplete$pca[[1]]$rotation %>% .[,1] %>% order(decreasing = F),1] %>% head(3)

# plotting 2D scatter (alternative version to the 3D scatter)
twodscatterplot = pcaAugDf %>% 
  ggplot(aes(x=.fittedPC1, y=.fittedPC2, color = tp, shape = br)) +
  geom_point(alpha = 0.75) +
  xlab(paste0("PC1", " (", round(expVar$var_exp[1] * 100, digits = 2), "%)")) +
  ylab(paste0("PC2", " (", round(expVar$var_exp[2] * 100, digits = 2), "%)")) +
  scale_color_manual(values = c("TP1" = "#E15759",
                                "TP2" = "#F28E2B",
                                "TP3" = "#4E79A7",
                                "TP4" = "#59A14F"),
                     # name="Time Point") +
                     name="Ponto de Coleta") +
  scale_shape_discrete(name="Bio. Replicate") +
#  scale_shape_discrete(name="Réplica Biológica") +
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2))

# saving previous plot
ggsave("./plots/exploratory_analysis_v3/spectronautPCA.svg",
       plot=twodscatterplot,
       width = 6, height = 3.5)

# arranging plots
fullpanel = ggarrange(plotlist = list(twodscatterplot, drawnht),
                      labels = "AUTO",
                      widths = c(1.5,1))

ggsave("./plots/exploratory_analysis_v3/spectronautPanelExploratory.png",
       plot=fullpanel,
       width = 7, height = 6,
       dpi = 300)

# plotting explained variance
cumuvarplot = expVar %>%
  pivot_longer(cols = contains("exp"),
               values_to = "var",
               names_to = "key") %>%
  mutate(key = case_when(key == "cum_var_exp" ~ "Cumulative Explained Variance",
                         key == "var_exp" ~ "Explained Variance",
                         TRUE ~ as.character(key))) %>% 
  ggplot(aes(x = factor(pc,levels = pc %>% unique() %>% str_sort(numeric = T)),
             y = var,
             group = key)) +
  geom_point() +
  geom_line() +
  facet_wrap(~key, scales = "free_y") +
  lims(y = c(0, 1)) +
  labs(y = "Variance",
       x = "Principal Components") +
  theme(axis.text.x = element_text(angle = 90))

# saving previous plot
ggsave("./plots/exploratory_analysis_v3/spectronautCumuvarplot.svg",
       plot=cumuvarplot,
       width = 8.5, height = 5)

ggsave("./plots/exploratory_analysis_v3/spectronautCumuvarplot.png",
       plot=cumuvarplot,
       width = 8.5, height = 5)

# exploring rna-seq dataset ####
# adding pseudocount
totrnatpmraw[,-1] = totrnatpmraw[,-1] + 1

# creating matrix
M = totrnatpmraw[ ,-1] %>% 
  drop_na() %>%
  as.matrix() %>%
  log10()

# defining annotation colors
legendCols = list(
  timepoints = c("TP1" = "#E15759",
                 "TP2" = "#F28E2B",
                 "TP3" = "#4E79A7",
                 "TP4" = "#59A14F"),
  bioReps = c("BR1" = adjustcolor(col = "#76B7B2", alpha.f = 0.2),
              "BR2" = adjustcolor(col = "#76B7B2", alpha.f = 0.5),
              "BR3" = adjustcolor(col = "#76B7B2", alpha.f = 1))
)

# defining annotations
tpannots = colnames(M) %>% str_replace(".*(TP[1-4]).*", "\\1")
brannots = colnames(M) %>% str_replace(".*(BR[1-3]).*", "\\1")

annot = HeatmapAnnotation(which = "col",
                          timepoints = anno_simple(tpannots,
                                                   col = legendCols$timepoints,
                                                   border = T),
                          bioReps = anno_simple(brannots,
                                                col = legendCols$bioReps,
                                                border = T),
                          annotation_label = c("Timepoint",
                                               "Bio. Replicate")
                          # annotation_label = c("Ponto de Coleta",
                          #                      "Réplica Biológica")
                          
)

# defining legs
legs = list(
  timepoints = Legend(title = "Time Point",
#   timepoints = Legend(title = "Ponto de Coleta",
                     at = legendCols$timepoints %>% names(),
                      legend_gp = gpar(fill = legendCols$timepoints %>% unname()),
                      border = "black"),
   bioReps = Legend(title = "Bio. Replicate",
#   bioReps = Legend(title = "Réplica Biológica",
                   at = legendCols$bioReps %>% names(),
                   legend_gp = gpar(fill = legendCols$bioReps %>% unname()),
                   border = "black")
)

# plotting heat map
ht = Heatmap(M, col = colfunct,
             top_annotation = annot,
             border = T,
             show_column_names = F,
             heatmap_legend_param = list(
               title = expression(Log[10](TPM+1)),
               border = T,
               direction = "horizontal"))

#svglite("./plots/exploratory_analysis_v3/totalrnaHeatmap.svg", width = 3, height = 4.5)
drawnht = grid.grabExpr(draw(ht,
                             heatmap_legend_side = "bottom",
                             annotation_legend_list = legs))
#dev.off()

# PCA
totrnatpmrawPCA = totrnatpmraw %>% 
  drop_na() %>% 
  transpose_tibble(target_id, id_col = "sample")

# removing genes with var == 0
filout = which(apply(X = totrnatpmrawPCA[,-1], MARGIN = 2, FUN = function(x) var(x)) == 0 ) %>% names()
totrnatpmrawPCA = totrnatpmrawPCA %>% 
  select(-all_of(filout))

# adapted from:
# https://tbradley1013.github.io/2018/02/01/pca-in-a-tidy-verse-framework/
# creating a dataframe containing all pca info we need
dfpcaComplete = totrnatpmrawPCA %>% 
  tidyr::nest(data=everything()) %>% 
  dplyr::mutate(pca = purrr::map(data, ~ prcomp(.x %>% select(-sample), 
                                                center = T, scale = T)),
                pca_aug = purrr::map2(pca, data, ~ broom::augment(.x, data = .y)))

# computing the explained variance for
# each of the principal components
expVar = dfpcaComplete %>% 
  unnest(pca_aug) %>% 
  summarize_at(.vars = vars(contains(".fittedPC")), .funs = list(~var(.))) %>% 
  gather(key = pc, value = variance) %>% 
  mutate(var_exp = variance/sum(variance),
         cum_var_exp = cumsum(var_exp),
         pc = str_replace(pc, ".fitted", ""))

# preparing df to plot PCA
pcaAugDf = dfpcaComplete %>%
  unnest(pca_aug) %>%
  mutate(tp = str_replace(sample,".*(TP[1-4]).*", "\\1"),
         br = str_replace(sample,".*(BR[1-3]).*", "\\1"))

# finding top three loadings in PC1
#dfpcaComplete$pca[[1]]$rotation[dfpcaComplete$pca[[1]]$rotation %>% .[,1] %>% order(decreasing = T),1] %>% head(3)

# finding bottom three loadings in PC1
#dfpcaComplete$pca[[1]]$rotation[dfpcaComplete$pca[[1]]$rotation %>% .[,1] %>% order(decreasing = F),1] %>% head(3)

# plotting 2D scatter (alternative version to the 3D scatter)
twodscatterplot = pcaAugDf %>% 
  ggplot(aes(x=.fittedPC1, y=.fittedPC2, color = tp, shape = br)) +
  geom_point(alpha = 0.75) +
  xlab(paste0("PC1", " (", round(expVar$var_exp[1] * 100, digits = 2), "%)")) +
  ylab(paste0("PC2", " (", round(expVar$var_exp[2] * 100, digits = 2), "%)")) +
  scale_color_manual(values = c("TP1" = "#E15759",
                                "TP2" = "#F28E2B",
                                "TP3" = "#4E79A7",
                                "TP4" = "#59A14F"),
                     name="Time Point") +
#                     name="Ponto de Coleta") +
  scale_shape_discrete(name="Bio. Replicate") +
#  scale_shape_discrete(name="Réplica Biológica") +
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2))

# saving previous plot
ggsave("./plots/exploratory_analysis_v3/totalrnaPCA.svg",
       plot=twodscatterplot,
       width = 6, height = 3.5)

# arranging plots
fullpanel = ggarrange(plotlist = list(twodscatterplot, drawnht),
                      labels = "AUTO",
                      widths = c(1.5,1))

ggsave("./plots/exploratory_analysis_v3/totalrnaPanelExploratory.png",
       plot=fullpanel,
       width = 7, height = 6,
       dpi = 300)

# plotting explained variance
cumuvarplot = expVar %>%
  pivot_longer(cols = contains("exp"),
               values_to = "var",
               names_to = "key") %>%
  mutate(key = case_when(key == "cum_var_exp" ~ "Cumulative Explained Variance",
                         key == "var_exp" ~ "Explained Variance",
                         TRUE ~ as.character(key))) %>% 
  ggplot(aes(x = factor(pc,levels = pc %>% unique() %>% str_sort(numeric = T)),
             y = var,
             group = key)) +
  geom_point() +
  geom_line() +
  facet_wrap(~key, scales = "free_y") +
  lims(y = c(0, 1)) +
  labs(y = "Variance",
       x = "Principal Components") +
  theme(axis.text.x = element_text(angle = 90))

# saving previous plot
ggsave("./plots/exploratory_analysis_v3/totalrnaCumuvarplot.svg",
       plot=cumuvarplot,
       width = 8.5, height = 5)

ggsave("./plots/exploratory_analysis_v3/totalrnaCumuvarplot.png",
       plot=cumuvarplot,
       width = 8.5, height = 5)

# exploring ribo-seq dataset ####
# adding pseudocount
ribornatpmraw[,-1] = ribornatpmraw[,-1] + 1

# creating matrix
M = ribornatpmraw[ ,-1] %>% 
  drop_na() %>%
  as.matrix() %>%
  log10()

# defining annotation colors
legendCols = list(
  timepoints = c("TP1" = "#E15759",
                 "TP2" = "#F28E2B",
                 "TP3" = "#4E79A7",
                 "TP4" = "#59A14F"),
  bioReps = c("BR1" = adjustcolor(col = "#76B7B2", alpha.f = 0.2),
              "BR2" = adjustcolor(col = "#76B7B2", alpha.f = 0.5),
              "BR3" = adjustcolor(col = "#76B7B2", alpha.f = 1))
)

# defining annotations
tpannots = colnames(M) %>% str_replace(".*(TP[1-4]).*", "\\1")
brannots = colnames(M) %>% str_replace(".*(BR[1-3]).*", "\\1")


annot = HeatmapAnnotation(which = "col",
                          timepoints = anno_simple(tpannots,
                                                   col = legendCols$timepoints,
                                                   border = T),
                          bioReps = anno_simple(brannots,
                                                col = legendCols$bioReps,
                                                border = T),
                          annotation_label = c("Timepoint",
                                               "Bio. Replicate")
)

# defining legs
legs = list(
  timepoints = Legend(title = "Time Point",
                      at = legendCols$timepoints %>% names(),
                      legend_gp = gpar(fill = legendCols$timepoints %>% unname()),
                      border = "black"),
  bioReps = Legend(title = "Bio. Replicate",
                   at = legendCols$bioReps %>% names(),
                   legend_gp = gpar(fill = legendCols$bioReps %>% unname()),
                   border = "black")
)

# plotting heat map
ht = Heatmap(M, col = colfunct,
             top_annotation = annot,
             border = T,
             show_column_names = F,
             heatmap_legend_param = list(
               title = expression(Log[10](Abundance)),
               border = T,
               direction = "horizontal"))

svglite("./plots/exploratory_analysis_v3/ribornaHeatmap.svg", width = 3, height = 4.5)
draw(ht,
     heatmap_legend_side = "bottom",
     annotation_legend_list = legs)
dev.off()

# PCA
ribornatpmrawPCA = ribornatpmraw %>% 
  drop_na() %>% 
  transpose_tibble(target_id, id_col = "sample")

# removing genes with var == 0
filout = which(apply(X = ribornatpmrawPCA[,-1], MARGIN = 2, FUN = function(x) var(x)) == 0 ) %>% names()
ribornatpmrawPCA = ribornatpmrawPCA %>% 
  select(-all_of(filout))

# adapted from:
# https://tbradley1013.github.io/2018/02/01/pca-in-a-tidy-verse-framework/
# creating a dataframe containing all pca info we need
dfpcaComplete = ribornatpmrawPCA %>% 
  tidyr::nest(data=everything()) %>% 
  dplyr::mutate(pca = purrr::map(data, ~ prcomp(.x %>% select(-sample), 
                                                center = T, scale = T)),
                pca_aug = purrr::map2(pca, data, ~ broom::augment(.x, data = .y)))

# computing the explained variance for
# each of the principal components
expVar = dfpcaComplete %>% 
  unnest(pca_aug) %>% 
  summarize_at(.vars = vars(contains(".fittedPC")), .funs = list(~var(.))) %>% 
  gather(key = pc, value = variance) %>% 
  mutate(var_exp = variance/sum(variance),
         cum_var_exp = cumsum(var_exp),
         pc = str_replace(pc, ".fitted", ""))

# preparing df to plot PCA
pcaAugDf = dfpcaComplete %>%
  unnest(pca_aug) %>%
  mutate(tp = str_replace(sample,".*(TP[1-4]).*", "\\1"),
         br = str_replace(sample,".*(BR[1-3]).*", "\\1"))

# finding top three loadings in PC1
dfpcaComplete$pca[[1]]$rotation[dfpcaComplete$pca[[1]]$rotation %>% .[,1] %>% order(decreasing = T),1] %>% head(3)

# finding bottom three loadings in PC1
dfpcaComplete$pca[[1]]$rotation[dfpcaComplete$pca[[1]]$rotation %>% .[,1] %>% order(decreasing = F),1] %>% head(3)

# plotting 2D scatter (alternative version to the 3D scatter)
twodscatterplot = pcaAugDf %>% 
  ggplot(aes(x=.fittedPC1, y=.fittedPC2, color = tp, shape = br)) +
  geom_point(alpha = 0.75) +
  xlab(paste0("PC1", " (", round(expVar$var_exp[1] * 100, digits = 2), "%)")) +
  ylab(paste0("PC2", " (", round(expVar$var_exp[2] * 100, digits = 2), "%)")) +
  scale_color_manual(values = c("TP1" = "#E15759",
                                "TP2" = "#F28E2B",
                                "TP3" = "#4E79A7",
                                "TP4" = "#59A14F"),
                     name="Time Point") +
  scale_shape_discrete(name="Bio. Replicate") +
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2))

# saving previous plot
ggsave("./plots/exploratory_analysis_v3/ribornaPCA.svg",
       plot=twodscatterplot,
       width = 6, height = 3.5)

# plotting explained variance
cumuvarplot = expVar %>%
  pivot_longer(cols = contains("exp"),
               values_to = "var",
               names_to = "key") %>%
  mutate(key = case_when(key == "cum_var_exp" ~ "Cumulative Explained Variance",
                         key == "var_exp" ~ "Explained Variance",
                         TRUE ~ as.character(key))) %>% 
  ggplot(aes(x = factor(pc,levels = pc %>% unique() %>% str_sort(numeric = T)),
             y = var,
             group = key)) +
  geom_point() +
  geom_line() +
  facet_wrap(~key, scales = "free_y") +
  lims(y = c(0, 1)) +
  labs(y = "Variance",
       x = "Principal Components") +
  theme(axis.text.x = element_text(angle = 90))

# saving previous plot
ggsave("./plots/exploratory_analysis_v3/ribornaCumuvarplot.svg",
       plot=cumuvarplot,
       width = 8.5, height = 5)

ggsave("./plots/exploratory_analysis_v3/ribornaCumuvarplot.png",
       plot=cumuvarplot,
       width = 8.5, height = 5)