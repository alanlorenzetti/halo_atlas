# 20210308 alorenzetti

# description ####
# this script will take the non normalized
# abundance values and compute time point-wise
# fold changes
# then the genes are going to be grouped
# according to the direction of change

# vars ####
# setting log2 FC minimum difference
# to change status
fcthr = 1.5
lfcthr = log2(1.5)

# if user wants NCBI locus_tags
# as primary, change this var to
# y
# locus_tag NCBI
ncbiLocusTag = "n"

# saving NCBI locus_tag dictionary
ncbiDict = nrtx %>% 
  separate_rows(locus_tag, sep = ",") %>% 
  dplyr::filter(str_detect(locus_tag, "VNG_RS")) %>% 
  group_by(representative) %>% 
  summarise(locus_tag = paste0(locus_tag, collapse = ",")) %>% 
  ungroup() %>% 
  mutate(locus_tag = str_replace(locus_tag, ",.*$", ""))

# starting processing ####
# 20210308 ####
# creating a new object to compute fold changes
# and changing
fcs = abund %>% 
  mutate(mean_lfc_rna_total_TP2vsTP1 = log2(mean_abundance_rna_total_TP2 / mean_abundance_rna_total_TP1),
         mean_lfc_rna_total_TP3vsTP2 = log2(mean_abundance_rna_total_TP3 / mean_abundance_rna_total_TP2),
         mean_lfc_rna_total_TP4vsTP3 = log2(mean_abundance_rna_total_TP4 / mean_abundance_rna_total_TP3),
         
         mean_lfc_rna_ribofraction_TP2vsTP1 = log2(mean_abundance_rna_ribofraction_TP2 / mean_abundance_rna_ribofraction_TP1),
         mean_lfc_rna_ribofraction_TP3vsTP2 = log2(mean_abundance_rna_ribofraction_TP3 / mean_abundance_rna_ribofraction_TP2),
         mean_lfc_rna_ribofraction_TP4vsTP3 = log2(mean_abundance_rna_ribofraction_TP4 / mean_abundance_rna_ribofraction_TP3),
         
         mean_lfc_protein_lysate_TP2vsTP1 = log2(mean_abundance_protein_lysate_TP2 / mean_abundance_protein_lysate_TP1),
         mean_lfc_protein_lysate_TP3vsTP2 = log2(mean_abundance_protein_lysate_TP3 / mean_abundance_protein_lysate_TP2),
         mean_lfc_protein_lysate_TP4vsTP3 = log2(mean_abundance_protein_lysate_TP4 / mean_abundance_protein_lysate_TP3)) %>% 
  select(c(locus_tag, contains("lfc")))

# checking distributions
# for alternative version of fcs
fcsalt = abund %>% 
  mutate(mean_lfc_rna_total_TP2 = log2(mean_abundance_rna_total_TP2 / mean_abundance_rna_total_TP1),
         mean_lfc_rna_total_TP3 = log2(mean_abundance_rna_total_TP3 / mean_abundance_rna_total_TP1),
         mean_lfc_rna_total_TP4 = log2(mean_abundance_rna_total_TP4 / mean_abundance_rna_total_TP1),
         
         mean_lfc_rna_ribofraction_TP2 = log2(mean_abundance_rna_ribofraction_TP2 / mean_abundance_rna_ribofraction_TP1),
         mean_lfc_rna_ribofraction_TP3 = log2(mean_abundance_rna_ribofraction_TP3 / mean_abundance_rna_ribofraction_TP1),
         mean_lfc_rna_ribofraction_TP4 = log2(mean_abundance_rna_ribofraction_TP4 / mean_abundance_rna_ribofraction_TP1),
         
         mean_lfc_protein_lysate_TP2 = log2(mean_abundance_protein_lysate_TP2 / mean_abundance_protein_lysate_TP1),
         mean_lfc_protein_lysate_TP3 = log2(mean_abundance_protein_lysate_TP3 / mean_abundance_protein_lysate_TP1),
         mean_lfc_protein_lysate_TP4 = log2(mean_abundance_protein_lysate_TP4 / mean_abundance_protein_lysate_TP1)) %>% 
  select(c(locus_tag, contains("lfc")))

fcsalt %>% 
  pivot_longer(cols = starts_with("mean"),
               names_to = c("libtype", "timepoint"),
               names_pattern = ".*_(.*_.*)_(.*)",
               values_to = "value") %>% 
  filter(libtype == "protein_lysate") %>% 
  ggplot(aes(x=value, color = timepoint)) +
  geom_density() +
  xlim(c(-6, 6)) +
  xlab("log2(Fold Change)")

# checking correlation of foldchanges between in house
# vs oneomics
oosnplots = list()
# setting correlation method
mt = "pearson"
  
# TP2 vs TP1
oosnplots$TP2 = oneomicsWide %>% 
  mutate(locus_tag = str_replace(locus_tag, "VNG", "VNG")) %>% 
  left_join(x = .,
            y = fcsalt,
            by = "locus_tag") %>% 
  select(c(locus_tag, contains("lfc"))) %>% 
  ggplot(aes(x = mean_lfc_protein_lysate_TP2.x,
             y = mean_lfc_protein_lysate_TP2.y)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  xlab("OneOmics") +
  ylab("Spectronaut")

# TP3 vs TP1
oosnplots$TP3 = oneomicsWide %>% 
  mutate(locus_tag = str_replace(locus_tag, "VNG", "VNG")) %>% 
  left_join(x = .,
            y = fcsalt,
            by = "locus_tag") %>% 
  select(c(locus_tag, contains("lfc"))) %>% 
  ggplot(aes(x = mean_lfc_protein_lysate_TP3.x,
             y = mean_lfc_protein_lysate_TP3.y)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  xlab("OneOmics") +
  ylab("Spectronaut")

# TP4 vs TP1
oosnplots$TP4 = oneomicsWide %>% 
  mutate(locus_tag = str_replace(locus_tag, "VNG", "VNG")) %>% 
  left_join(x = .,
            y = fcsalt,
            by = "locus_tag") %>% 
  select(c(locus_tag, contains("lfc"))) %>% 
  ggplot(aes(x = mean_lfc_protein_lysate_TP4.x,
             y = mean_lfc_protein_lysate_TP4.y)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  xlab("OneOmics") +
  ylab("Spectronaut")

ggarrange(plotlist = oosnplots,
          labels = c("TP2 vs. TP1",
                     "TP3 vs. TP1",
                     "TP4 vs. TP1"),
          nrow = 1,
          ncol = 3)

# checking correlation of fold changes between protein vs mrna ####
mrna_prot = list()

# TP2 vs TP1
mrna_prot$TP2 = fcsalt %>% 
  ggplot(aes(x = mean_lfc_rna_total_TP2,
             y = mean_lfc_protein_lysate_TP2)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  xlab("mRNA") +
  ylab("Protein")

# TP3 vs TP1
mrna_prot$TP3 = fcsalt %>% 
  ggplot(aes(x = mean_lfc_rna_total_TP3,
             y = mean_lfc_protein_lysate_TP3)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  xlab("mRNA") +
  ylab("Protein")

# TP4 vs TP1
mrna_prot$TP4 = fcsalt %>% 
  ggplot(aes(x = mean_lfc_rna_total_TP4,
             y = mean_lfc_protein_lysate_TP4)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  xlab("mRNA") +
  ylab("Protein")

ggarrange(plotlist = mrna_prot,
          labels = c("TP2 vs. TP1",
                     "TP3 vs. TP1",
                     "TP4 vs. TP1"),
          nrow = 1,
          ncol = 3)

# replacing locus_tags if requested by the user
if(ncbiLocusTag == "y"){
  fcs = fcs %>% 
    left_join(x = ., y = ncbiDict,
              by = c("locus_tag" = "representative")) %>% 
    dplyr::mutate(locus_tag = locus_tag.y) %>% 
    dplyr::select(-locus_tag.y) %>% 
    filter(!is.na(locus_tag))
}

# creating new columns to store change status
fcs = fcs %>% 
  mutate(change_rna_total_TP2vsTP1 = case_when(mean_lfc_rna_total_TP2vsTP1 > lfcthr ~ "Up",
                                               mean_lfc_rna_total_TP2vsTP1 < -lfcthr ~ "Down",
                                               abs(mean_lfc_rna_total_TP2vsTP1) < lfcthr ~ "Flat",
                                               TRUE ~ "?"),
         change_rna_total_TP3vsTP2 = case_when(mean_lfc_rna_total_TP3vsTP2 > lfcthr ~ "Up",
                                               mean_lfc_rna_total_TP3vsTP2 < -lfcthr ~ "Down",
                                               abs(mean_lfc_rna_total_TP3vsTP2) < lfcthr ~ "Flat",
                                               TRUE ~ "?"),
         change_rna_total_TP4vsTP3 = case_when(mean_lfc_rna_total_TP4vsTP3 > lfcthr ~ "Up",
                                               mean_lfc_rna_total_TP4vsTP3 < -lfcthr ~ "Down",
                                               abs(mean_lfc_rna_total_TP4vsTP3) < lfcthr ~ "Flat",
                                               TRUE ~ "?")) %>% 
  
  mutate(change_rna_ribofraction_TP2vsTP1 = case_when(mean_lfc_rna_ribofraction_TP2vsTP1 > lfcthr ~ "Up",
                                                      mean_lfc_rna_ribofraction_TP2vsTP1 < -lfcthr ~ "Down",
                                                      abs(mean_lfc_rna_ribofraction_TP2vsTP1) < lfcthr ~ "Flat",
                                                      TRUE ~ "?"),
         change_rna_ribofraction_TP3vsTP2 = case_when(mean_lfc_rna_ribofraction_TP3vsTP2 > lfcthr ~ "Up",
                                                      mean_lfc_rna_ribofraction_TP3vsTP2 < -lfcthr ~ "Down",
                                                      abs(mean_lfc_rna_ribofraction_TP3vsTP2) < lfcthr ~ "Flat",
                                                      TRUE ~ "?"),
         change_rna_ribofraction_TP4vsTP3 = case_when(mean_lfc_rna_ribofraction_TP4vsTP3 > lfcthr ~ "Up",
                                                      mean_lfc_rna_ribofraction_TP4vsTP3 < -lfcthr ~ "Down",
                                                      abs(mean_lfc_rna_ribofraction_TP4vsTP3) < lfcthr ~ "Flat",
                                                      TRUE ~ "?")) %>% 
  
  mutate(change_protein_lysate_TP2vsTP1 = case_when(mean_lfc_protein_lysate_TP2vsTP1 > lfcthr ~ "Up",
                                                    mean_lfc_protein_lysate_TP2vsTP1 < -lfcthr ~ "Down",
                                                    abs(mean_lfc_protein_lysate_TP2vsTP1) < lfcthr ~ "Flat",
                                                    TRUE ~ "?"),
         change_protein_lysate_TP3vsTP2 = case_when(mean_lfc_protein_lysate_TP3vsTP2 > lfcthr ~ "Up",
                                                    mean_lfc_protein_lysate_TP3vsTP2 < -lfcthr ~ "Down",
                                                    abs(mean_lfc_protein_lysate_TP3vsTP2) < lfcthr ~ "Flat",
                                                    TRUE ~ "?"),
         change_protein_lysate_TP4vsTP3 = case_when(mean_lfc_protein_lysate_TP4vsTP3 > lfcthr ~ "Up",
                                                    mean_lfc_protein_lysate_TP4vsTP3 < -lfcthr ~ "Down",
                                                    abs(mean_lfc_protein_lysate_TP4vsTP3) < lfcthr ~ "Flat",
                                                    TRUE ~ "?")) %>% 
  
  # creating a code for the time course changes
  mutate(timecourse_code_rna_total = paste(change_rna_total_TP2vsTP1,
                                           change_rna_total_TP3vsTP2,
                                           change_rna_total_TP4vsTP3,
                                           sep = "_"),
         timecourse_code_rna_ribofraction = paste(change_rna_ribofraction_TP2vsTP1,
                                                  change_rna_ribofraction_TP3vsTP2,
                                                  change_rna_ribofraction_TP4vsTP3,
                                                  sep = "_"),
         timecourse_code_protein_lysate = paste(change_protein_lysate_TP2vsTP1,
                                                change_protein_lysate_TP3vsTP2,
                                                change_protein_lysate_TP4vsTP3,
                                                sep = "_"))

# organizing data structure to get insights
# about the genes following each one of the 
# change codes
tcCodes = list()

# everything
tcCodes[["all"]] = fcs %>%
  select(c(locus_tag,
           starts_with("timecourse_code")))

# subsetting and grouping rna, ribo, and protein
# rna
tcCodes[["rna_total"]] = tcCodes[["all"]] %>%
  select(c(locus_tag,
           contains("rna_total"))) %>% 
  group_by(timecourse_code_rna_total) %>% 
  summarise(entries_rna_total = paste0(locus_tag, collapse = ","))

# ribo 
tcCodes[["rna_ribofraction"]] = tcCodes[["all"]] %>%
  select(c(locus_tag,
           contains("rna_ribofraction"))) %>% 
  group_by(timecourse_code_rna_ribofraction) %>% 
  summarise(entries_rna_ribofraction = paste0(locus_tag, collapse = ","))

# protein
tcCodes[["protein_lysate"]] = tcCodes[["all"]] %>%
  select(c(locus_tag,
           contains("protein_lysate"))) %>% 
  group_by(timecourse_code_protein_lysate) %>% 
  summarise(entries_protein_lysate = paste0(locus_tag, collapse = ","))

# unifying datasets 
tcCodes[["unified_written"]] = full_join(x = tcCodes[["rna_total"]],
                                         y = tcCodes[["rna_ribofraction"]],
                                         by = c("timecourse_code_rna_total" = "timecourse_code_rna_ribofraction"))

tcCodes[["unified_written"]] = full_join(x = tcCodes[["unified_written"]],
                                         y = tcCodes[["protein_lysate"]],
                                         by = c("timecourse_code_rna_total" = "timecourse_code_protein_lysate")) %>% 
  dplyr::rename(timecourse_code = "timecourse_code_rna_total")

# ordering data set
tcCodes[["unified_written"]] = tcCodes[["unified_written"]][mixedorder(tcCodes[["unified_written"]]$timecourse_code,
                                                                       decreasing = T),]
# creating cols to store intersections
# and then finding them
tcCodes$unified_written$intersection_total_prot = NA_character_
tcCodes$unified_written$intersection_total_ribo = NA_character_
tcCodes$unified_written$intersection_prot_ribo = NA_character_
tcCodes$unified_written$intersection_all = NA_character_

for(i in 1:dim(tcCodes$unified_written)[1]){
  
  rna = tcCodes$unified_written$entries_rna_total[i] %>%
    str_split(pattern = ",", simplify = F) %>%
    unlist()
  
  ribo = tcCodes$unified_written$entries_rna_ribofraction[i] %>%
    str_split(pattern = ",", simplify = F) %>%
    unlist()
  
  prot = tcCodes$unified_written$entries_protein_lysate[i] %>%
    str_split(pattern = ",", simplify = F) %>%
    unlist()
  
  rnaprot = base::intersect(rna, prot) ; if(identical(rnaprot, character(0))){rnaprot = NA_character_}
  rnaribo = base::intersect(rna, ribo) ; if(identical(rnaribo, character(0))){rnaribo = NA_character_}
  protribo = base::intersect(prot, ribo) ; if(identical(protribo, character(0))){protribo = NA_character_}
  all = base::intersect(base::intersect(rna, prot), ribo) ; if(identical(all, character(0))){all = NA_character_}
  
  tcCodes$unified_written$intersection_total_prot[i] = paste0(rnaprot, collapse = ",")
  tcCodes$unified_written$intersection_total_ribo[i] = paste0(rnaribo, collapse = ",")
  tcCodes$unified_written$intersection_prot_ribo[i] = paste0(protribo, collapse = ",")
  tcCodes$unified_written$intersection_all[i] = paste0(all, collapse = ",")
}

# adjusting NAs of intersection cols
# and headers
tcCodes[["unified_written"]] = tcCodes[["unified_written"]] %>% 
  mutate(across(.cols = starts_with("intersection"),
                .fns = ~ case_when(.x  == "NA" ~ NA_character_,
                                   TRUE ~ as.character(.x))))

# renaming cols
colnames(tcCodes[["unified_written"]]) = c("timecourse_code",
                                           "Entries mRNA",
                                           "Entries RiboSeq",
                                           "Entries Protein",
                                           "Intersection mRNA Protein",
                                           "Intersection mRNA RiboSeq",
                                           "Intersection Protein RiboSeq",
                                           "Intersection mRNA Protein RiboSeq")

# creating a count matrix based in the above matrix
tcCodes[["unified_counts"]] = tcCodes[["unified_written"]] %>% 
  mutate(across(.cols = -timecourse_code,
                .fns = ~ str_count(string = .x, pattern = "VNG_")))

tcCodes[["webObj"]] = map2_df(tcCodes[["unified_counts"]],
                              tcCodes[["unified_written"]],
                              ~ paste(.x, .y)) %>% 
  mutate(timecourse_code = str_replace(timecourse_code, " .*$", "")) %>% 
  mutate(across(.cols = -timecourse_code,
                .fns = ~ case_when(str_detect(.x, "NA NA") ~ NA_character_,
                                   TRUE ~ str_replace(string = .x,
                                                      pattern = '^(.*) (.*)$',
                                                      replacement = '<a href="#" onclick="alert(\'\\2\');">\\1</a>'))))  %>% 
  mutate(across(.cols = -timecourse_code,
                .fns = ~ case_when(is.na(.x) ~ NA_character_,
                                   TRUE ~ str_replace_all(string = .x,
                                                          pattern = ',',
                                                          replacement = ' '))))


tcCodes[["unified_written"]] %>% 
  write_tsv(x = .,
            file = "results/20210311_globalFCanalysis_entries_table.tsv")

tcCodes[["unified_counts"]] %>% 
  write_tsv(x = .,
            file = "results/20210311_globalFCanalysis_counts_table.tsv")

# generating and saving datatable obj
tcCodes[["webObj"]] %>% 
  datatable(data = ., escape = F) %>% 
  DT::saveWidget(widget = .,
                 file = "20210311_globalFCanalysis_int_table.html",
                 selfcontained = T)

# creating and generating histograms for lfcs tpi vs. tp1
if(!dir.exists("plots/20210312_lfc_hists")){dir.create("plots/20210312_lfc_hists")}

# new tidy lfc object
fcs2 = abund %>% 
  mutate(mean_lfc_rna_total_TP2vsTP1 = log2(mean_abundance_rna_total_TP2 / mean_abundance_rna_total_TP1),
         mean_lfc_rna_total_TP3vsTP1 = log2(mean_abundance_rna_total_TP3 / mean_abundance_rna_total_TP1),
         mean_lfc_rna_total_TP4vsTP1 = log2(mean_abundance_rna_total_TP4 / mean_abundance_rna_total_TP1),
         
         mean_lfc_rna_ribofraction_TP2vsTP1 = log2(mean_abundance_rna_ribofraction_TP2 / mean_abundance_rna_ribofraction_TP1),
         mean_lfc_rna_ribofraction_TP3vsTP1 = log2(mean_abundance_rna_ribofraction_TP3 / mean_abundance_rna_ribofraction_TP1),
         mean_lfc_rna_ribofraction_TP4vsTP1 = log2(mean_abundance_rna_ribofraction_TP4 / mean_abundance_rna_ribofraction_TP1),
         
         mean_lfc_protein_lysate_TP2vsTP1 = log2(mean_abundance_protein_lysate_TP2 / mean_abundance_protein_lysate_TP1),
         mean_lfc_protein_lysate_TP3vsTP1 = log2(mean_abundance_protein_lysate_TP3 / mean_abundance_protein_lysate_TP1),
         mean_lfc_protein_lysate_TP4vsTP1 = log2(mean_abundance_protein_lysate_TP4 / mean_abundance_protein_lysate_TP1)) %>% 
  select(c(locus_tag, contains("lfc"))) %>% 
  pivot_longer(cols = starts_with("mean"),
               names_to = c("libtype", "timepointContrast"),
               names_pattern = "mean_(.*_.*)_(TP.*)",
               values_to = "lfc")

# plotting and saving
plfc = fcs2 %>%
  ggplot(aes(x = lfc)) +
  geom_histogram(color = "black", fill = "white") +
  facet_wrap(timepointContrast ~ libtype )
ggsave(plot = plfc,
       filename = "plots/20210312_lfc_hists/lfc_hists.png",
       dpi = 300,
       unit = "in",
       width = 10,
       height = 7)

# 20210314 ####
# THIS IS NOT AN AUTOMATIC STEP YET
# save the 27x27 matrix with locus_tags from 
# 20210314_compendium_globalFCanalysis_counts_table_RV.xlsx to
# 20210314_compendium_globalFCanalysis_counts_table_RV.csv to
# and then proceed
# reading csv
df = read_csv(file = "data/20210314_compendium_globalFCanalysis_counts_table_RV.csv",
              col_names = F)

# getting a vector of genes interacting with lsm
lsm_lt = lsmGenes %>%
  filter(lsmSense == "yes") %>%
  select(name) %>%
  unlist(use.names = F)

# creating output matrix
outmatrix = matrix(nrow = 27,
                   ncol = 27)

# iterating over original tibble
for(i in 1:dim(df)[1]){
  for(j in 1:dim(df)[2]){
    if(!is.na(df[i,j] %>% unlist(use.names = F))){
      cell = df[i,j] %>% 
        unlist(use.names = F) %>%
        str_split(pattern = " ") %>%
        unlist()
      outmatrix[i,j] = sum(cell %in% lsm_lt)
    }else{
      outmatrix[i,j] = 0
    }
  }
}

# plots requested by rvencio
# plotting ####
breaks = 10^(-10:10)
minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))

# setting correlation method to compare
mt = "pearson"

# protein abundance in function of mRNA ####
protvsmrna = abundLong %>% 
  filter(timepoint != "TP0") %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = protein_lysate,
             x = rna_total)) +
  #    geom_density2d() +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  facet_wrap(~ timepoint) +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance") +
  # ylab("Abundância de Proteínas") +
  # xlab("Abundância de mRNAs") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  annotation_logticks()

# slide prot from TPi+1 mRNA from TPi
# TP2 vs TP1
slide1 = abund %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP2,
             x = mean_abundance_rna_total_TP1)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  annotation_logticks() +
  ggtitle("Protein TP2 vs. mRNA TP1")

# TP3 vs TP2
slide2 = abund %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP3,
             x = mean_abundance_rna_total_TP2)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  annotation_logticks() +
  ggtitle("Protein TP3 vs. mRNA TP2")

# TP4 vs TP3
slide3 = abund %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP4,
             x = mean_abundance_rna_total_TP3)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  annotation_logticks() +
  ggtitle("Protein TP4 vs. mRNA TP3")

# arranging plots
protvsmrnaslide = ggarrange(slide1, slide2, slide3,
                            ncol = 3, nrow = 1)

# 20210323 ####
# plots requested by rvencio
# plotting ####
breaks = 10^(-10:10)
minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))

# setting correlation method to compare
mt = "pearson"

# protein abundance in function of mRNA
protvsmrna = abundLong %>% 
  filter(timepoint != "TP0") %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = protein_lysate,
             x = rna_total)) +
  geom_density2d(color = "black") +
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  facet_wrap(~ timepoint) +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance") +
  # ylab("Abundância de Proteínas") +
  # xlab("Abundância de mRNAs") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  annotation_logticks()

# slide prot from TPi+1 mRNA from TPi
# TP2 vs TP1
slide1 = abund %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP2,
             x = mean_abundance_rna_total_TP1)) +
  geom_density2d(color = "black") + 
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  annotation_logticks() +
  ggtitle("Protein TP2 vs. mRNA TP1")

# TP3 vs TP2
slide2 = abund %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP3,
             x = mean_abundance_rna_total_TP2)) +
  geom_density2d(color = "black") + 
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  annotation_logticks() +
  ggtitle("Protein TP3 vs. mRNA TP2")

# TP4 vs TP3
slide3 = abund %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP4,
             x = mean_abundance_rna_total_TP3)) +
  geom_density2d(color = "black") +  
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  annotation_logticks() +
  ggtitle("Protein TP4 vs. mRNA TP3")

# arranging plots
protvsmrnaslide = ggarrange(slide1, slide2, slide3,
                            ncol = 3, nrow = 1)

# "NB" plots####
# setting correlation method to compare
mt = "pearson"

tmpdf = abund %>% 
  dplyr::select(!starts_with("se")) %>% 
  mutate(mean_ratio_translational_efficiency_TP1 = mean_abundance_protein_lysate_TP1 / mean_abundance_rna_total_TP1,
         mean_ratio_translational_efficiency_TP2 = mean_abundance_protein_lysate_TP2 / mean_abundance_rna_total_TP2,
         mean_ratio_translational_efficiency_TP3 = mean_abundance_protein_lysate_TP3 / mean_abundance_rna_total_TP3,
         mean_ratio_translational_efficiency_TP4 = mean_abundance_protein_lysate_TP4 / mean_abundance_rna_total_TP4)

tmpdflong = tmpdf %>% 
  pivot_longer(cols = contains(c("abundance", "ratio")),
               names_pattern = "^(.*)_.*_(.*_.*)_(.*)$",
               names_to = c("measure", "libtype", "timepoint"),
               values_to = "abundance") %>% 
  pivot_wider(names_from = c("measure"),
              values_from = "abundance")

# protein abundance in function of mRNA
protvsmrna = tmpdflong %>% 
  filter(timepoint != "TP0") %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(translational_efficiency) ,
             x = log2(protein_lysate)/2 + log2(rna_total)/2)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  facet_wrap(~ timepoint) +
  ylab("log2(Protein/mRNA)") +
  xlab("log2(Protein)/2 + log2(mRNA)/2")
  # ylab("Abundância de Proteínas") +
  # xlab("Abundância de mRNAs") +

# slide prot from TPi+1 mRNA from TPi
# TP2 vs TP1
slide1 = abund %>% 
  ggplot(aes(y = log2(mean_abundance_protein_lysate_TP2/mean_abundance_rna_total_TP1),
             x = log2(mean_abundance_protein_lysate_TP2)/2 + log2(mean_abundance_rna_total_TP1)/2)) +
  geom_point(alpha = 0.25) +  
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  ylab("log2(Protein/mRNA)") +
  xlab("log2(Protein)/2 + log2(mRNA)/2")  +
  ggtitle("Protein TP2 vs. mRNA TP1")

# TP3 vs TP2
slide2 = abund %>% 
  ggplot(aes(y = log2(mean_abundance_protein_lysate_TP3/mean_abundance_rna_total_TP2),
             x = log2(mean_abundance_protein_lysate_TP3)/2 + log2(mean_abundance_rna_total_TP1)/2)) +
  geom_point(alpha = 0.25) +  
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  ylab("log2(Protein/mRNA)") +
  xlab("log2(Protein)/2 + log2(mRNA)/2")  +
  ggtitle("Protein TP3 vs. mRNA TP2")

# TP4 vs TP3
slide3 = abund %>% 
  ggplot(aes(y = log2(mean_abundance_protein_lysate_TP4/mean_abundance_rna_total_TP3),
             x = log2(mean_abundance_protein_lysate_TP4)/2 + log2(mean_abundance_rna_total_TP3)/2)) +
  geom_point(alpha = 0.25) +  
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = mt) +
  ylab("log2(Protein/mRNA)") +
  xlab("log2(Protein)/2 + log2(mRNA)/2")  +
  ggtitle("Protein TP4 vs. mRNA TP3")

# arranging plots
protvsmrnaslide = ggarrange(slide1, slide2, slide3,
                            ncol = 3, nrow = 1)


