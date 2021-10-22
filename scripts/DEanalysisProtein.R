# alorenzetti 20210818

# description ####
# this script will take protein groups intensity
# matrix and compute differential expression
# using proDA from Bioconductor

# setting up thresholds ####
padjthreshold = 0.05
log2fcthreshold = 1
borderlinezero = 0.25

# loading required SN file file ####
spectro = read_tsv("data/Halo_WholeLysate_HaloIonLibraryVNG_merged_final_2018-04-07_SW100_unnormalized_Protein_Report_C4.xls")

# repeating protein group filtering 
# from parseProteomicsSpectronaut script
# adjusting colnames
colnames(spectro)[-1] = sub("^\\[.* ", "", colnames(spectro)[-1])
colnames(spectro)[-1] = sub("^.*_?(BR[1-3])(TP[1-4]).*_.*_.*_(r0[1-3]).*$",
                            "\\2_\\1_\\3",
                            colnames(spectro)[-1],
                            perl = T)
colnames(spectro)[-1] = paste0(rep("lysate_", 36),
                               colnames(spectro)[-1])
colnames(spectro)[1] = "locus_tag"

# removing first entry line
spectro = spectro[-1, ]

# warnings will be thrown due to
# missing entries ("Filtered" value)
spectro = spectro %>%
  mutate_at(.vars = vars(contains("BR")),
            .funs = as.numeric)

# selection of entries ####
# we have to filter the table due to the existence
# of peptides matching more than one protein
# our simplification is:
# the members of a protein group can exist in
# the form of an individual observation.
# if only one member of a protein group exists individually,
# protein group values and individual protein values
# are tallied up;
# if more than two members of a protein group exist individually
# the protein group is discarded
# if no member of a protein group exists individually
# the protein group is kept (first locus_tag is preserved)

# finding out solos and protein groups
solo = spectro %>% 
  filter(!str_detect(locus_tag, ","))
pg = spectro %>% 
  filter(str_detect(locus_tag, ","))

ex = NULL
pgvspg = NULL
for(i in 1:dim(pg)[1]){
  curline = pg$locus_tag[i] %>% 
    str_split(pattern = ",") %>% 
    unlist()
  for(j in 1:length(curline)){
    ex[j] = curline[j] %in% solo$locus_tag
    names(ex)[j] = curline[j]
  }
  if(sum(ex) == 0){pass = "pass0"; repres = names(ex)[1]}
  if(sum(ex) == 1){pass = "pass1"; repres = names(which(ex == TRUE))}
  if(sum(ex) > 1){pass = "fail"; repres = NA_character_} 
  
  y = 1:dim(pg)[1]
  for(k in y[-i]){
    curpg = pg$locus_tag[k] %>% 
      str_split(pattern = ",") %>% 
      unlist()
    pgvspg = c(pgvspg, sum(names(ex) %in% curpg))
  }
  
  npg = sum(pgvspg > 0)
  if(sum(pgvspg) == 0){passpg = "pass"}
  if(sum(pgvspg) > 0){passpg = paste0("fail", npg)}
  
  if((pass == "pass0" | pass == "pass1") & passpg == "pass"){
    pg$locus_tag[i] = repres
  }
  
  ex = NULL
  pgvspg = NULL
}

# filtering out those protein groups
# that failed our simplification
pg = pg %>% 
  filter(!str_detect(locus_tag, ","))

# VNG5030G (GvpA) and VNG5033G (GvpN) were manually included
# since they are protein groups contained within
# another protein group but satisfy the nonambiguity criterion
# e.g. the protein groups come from isoforms of GvpA1a, GvpA1b, and GvpA2
gvpAN = spectro %>% 
  filter(str_detect(string = locus_tag, pattern = "^VNG5030G|^VNG5033G")) %>% 
  mutate(locus_tag = str_replace(string = locus_tag, pattern = ",.*", replacement = ""))

# unifying individual proteins to
# simplified protein groups
soloPlusSPG = bind_rows(solo, pg, gvpAN) %>% 
  arrange(locus_tag) %>% 
  group_by(locus_tag) %>%
  dplyr::summarise(across(.cols = starts_with("lysate"),
                          .fns = ~ sum(.x, na.rm = T))) %>% 
  ungroup() %>% 
  dplyr::mutate(across(.cols = starts_with("lysate"),
                       .fns = ~ case_when(.x == 0 ~ NA_real_,
                                          TRUE ~ as.numeric(.x))))

spectro = soloPlusSPG
snv2 = spectro

# following with proDA vignette
# https://www.bioconductor.org/packages/devel/bioc/vignettes/proDA/inst/doc/Introduction.html
snv2M = snv2[,-1] %>%
  as.matrix()
rownames(snv2M) = snv2$locus_tag

# according to the workflow
# raw data present linear mean-variance relationship
# and should be transformed by log function
# I would rather check first if our dataset
# follow this trend
snv2M = log2(snv2M)

# # first QC step to check missing values in samples
# barplot(colSums(is.na(snv2M)),
#         ylab = "# missing values",
#         xlab = "Sample 1 to 36")
# 
# boxplot(snv2M,
#         ylab = "Intensity Distribution",
#         xlab = "Sample 1 to 36")

# normalizing using median normalization
normalized_abundance_matrix = median_normalization(snv2M)

# exploratory analysis ####
# creating distance matrix
# da = dist_approx(normalized_abundance_matrix)
# 
# # plotting matrix
# sel = 1:36
# plot_mat = as.matrix(da$mean)[sel, sel]
# 
# # Remove diagonal elements, so that the colorscale is not distorted
# plot_mat[diag(36) == 1] = NA
# 
# # 95% conf interval is approx `sd * 1.96`
# uncertainty = matrix(paste0(" ± ",round(as.matrix(da$sd * 1.96)[sel, sel], 1)), nrow=36)
# pheatmap::pheatmap(plot_mat, 
#                    cluster_rows = FALSE, cluster_cols = FALSE,
#                    display_numbers= uncertainty,
#                    number_color = "black")

# The best way to create this data.frame depends on the column naming scheme
sample_info_df = data.frame(name = colnames(normalized_abundance_matrix),
                             stringsAsFactors = FALSE)
sample_info_df$timepoint = str_replace(colnames(normalized_abundance_matrix),
                                       pattern = ".*_(.*)_.*_.*",
                                       replacement = "\\1")
sample_info_df$replicate = c(rep(1:3,4), rep(4:6,4), rep(7:9,4))

# fit = proDA(normalized_abundance_matrix, design = ~ timepoint, 
#             col_data = sample_info_df, reference_level = "TP1")
# 
# # getting differential abundance
# test_res = test_diff(fit = fit, "conditionTP4")

# result tables for contrasts ####
resultsDEProt = list()
for(i in 2:4){
  # fitting model
  resultsDEProt[[paste0("protein_TP", i, "_vs_TP", 1)]][["model"]] = proDA(normalized_abundance_matrix,
                                                                design = ~ timepoint,
                                                                col_data = sample_info_df,
                                                                reference_level = paste0("TP", 1))
  
  # getting differential abundance
  resultsDEProt[[paste0("protein_TP", i, "_vs_TP", 1)]][["results"]] = test_diff(fit = resultsDEProt[[paste0("protein_TP", i, "_vs_TP", 1)]][["model"]],
                                                                                 paste0("timepointTP", i)) %>% 
    as_tibble()
}

for(i in 3:4){
  # fitting model
  resultsDEProt[[paste0("protein_TP", i, "_vs_TP", i-1)]][["model"]] = proDA(normalized_abundance_matrix,
                                                                           design = ~ timepoint,
                                                                           col_data = sample_info_df,
                                                                           reference_level = paste0("TP", i-1))
  
  # getting differential abundance
  resultsDEProt[[paste0("protein_TP", i, "_vs_TP", i-1)]][["results"]] = test_diff(fit = resultsDEProt[[paste0("protein_TP", i, "_vs_TP", i-1)]][["model"]],
                                                                                 paste0("timepointTP", i)) %>% 
    as_tibble()
}

# classification of genes regarding DE status ####
resultsFinDEProt = list()
for(i in 2:4){
  varName = paste0("sig_protein","_TP",i,"_vs_TP", 1)
  varName2 = paste0("borderlineZero_protein","_TP",i,"_vs_TP", 1)
  resultsFinDEProt[[paste0("protein","_TP", i, "_vs_TP", 1)]] = resultsDEProt[[paste0("protein_TP", i, "_vs_TP", 1)]][["results"]] %>% 
    mutate(!!varName := case_when(abs(diff) >= log2fcthreshold & adj_pval < padjthreshold ~ "yes",
                                  TRUE ~ "no")) %>% 
    mutate(!!varName2 := case_when(abs(diff) < borderlinezero ~ "yes", 
                                   TRUE ~ "no"))
}

for(i in 3:4){
  varName = paste0("sig_protein","_TP",i,"_vs_TP", i-1)
  varName2 = paste0("borderlineZero_protein","_TP",i,"_vs_TP", i-1)
  resultsFinDEProt[[paste0("protein","_TP", i, "_vs_TP", i-1)]] = resultsDEProt[[paste0("protein_TP", i, "_vs_TP", i-1)]][["results"]] %>% 
    mutate(!!varName := case_when(abs(diff) >= log2fcthreshold & adj_pval < padjthreshold ~ "yes",
                                  TRUE ~ "no")) %>% 
    mutate(!!varName2 := case_when(abs(diff) < borderlinezero ~ "yes", 
                                   TRUE ~ "no"))
}

# parsing final results
unifiedFinDEProt = list()
for(i in 2:4){
  unifiedFinDEProt[[paste0("TP", i, "_vs_TP", 1)]] = resultsFinDEProt[[paste0("protein","_TP", i, "_vs_TP", 1)]] %>% 
    rename("id" = name) %>% 
    left_join(x = .,
              y = nrtxsep %>% select(-product),
              by = c("id" = "locus_tag")) %>% 
    filter(!is.na(representative)) %>% 
    relocate(representative) %>% 
    select(-id) %>% 
    rename_with(.fn = ~ str_replace(string = .x,
                                    pattern = "_TP.*$",
                                    replacement = ""),
                .cols = matches("TP")) %>% 
    rename("mean_lfc_protein_lysate" = diff,
           "se_lfc_protein_lysate" = se,
           "locus_tag" = representative)
}
for(i in 3:4){
  unifiedFinDEProt[[paste0("TP", i, "_vs_TP", i-1)]] = resultsFinDEProt[[paste0("protein","_TP", i, "_vs_TP", i-1)]] %>% 
    rename("id" = name) %>% 
    left_join(x = .,
              y = nrtxsep %>% select(-product),
              by = c("id" = "locus_tag")) %>% 
    filter(!is.na(representative)) %>% 
    relocate(representative) %>% 
    select(-id) %>% 
    rename_with(.fn = ~ str_replace(string = .x,
                                    pattern = "_TP.*$",
                                    replacement = ""),
                .cols = matches("TP")) %>% 
    rename("mean_lfc_protein_lysate" = diff,
           "se_lfc_protein_lysate" = se,
           "locus_tag" = representative)
}
