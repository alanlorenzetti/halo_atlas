# alorenzetti 20210616

# description ####
# this script will take data
# from Galal et al. 2021 (Table S1) 
# on Transcript Processing Sites
# and parse it accordingly

# loading original file and parsing ####
# Table_S1___TPS_to_cds_Hsal.tsv

# getting locus_tag of genes with associated TPS
tpsoritbl = read_tsv(file = "data/Table_S1___TPS_to_cds_Hsal.tsv")
tpsncbilt = tpsoritbl %>%
  group_by(Common.Name) %>%
  summarise(count = n())

# creating tibble to store data
tpsncbi = tibble(ncbi_locus_tag = tpsncbilt$Common.Name,
                 tps = tpsncbilt$count)

# translating locus_tag
tpscount = left_join(x = nrtx %>%
                       separate_rows(sep = ",", locus_tag),
                     y = tpsncbi,
                     by = c("locus_tag" = "ncbi_locus_tag")) %>% 
  mutate(tps = case_when(is.na(tps) ~ as.integer(0),
                         TRUE ~ tps)) %>% 
  group_by(representative) %>% 
  summarise(tps = max(tps))

tps = tpscount %>% 
  mutate(tps = case_when(tps > 0 ~ "yes",
                         TRUE ~ "no"))

# creating a gff files for genome browser
tpsgff = tibble(acc = tpsoritbl$Replicon,
                src = "Ibrahim_et_al_2021",
                genebiotype = "misc_feature",
                start = tpsoritbl$Position,
                end = tpsoritbl$Position,
                score = ".",
                strand = tpsoritbl$Strand,
                misc8 = ".",
                att = paste0("ID=", tpsoritbl$ID, ";",
                             "Name=", tpsoritbl$ID, ";",
                             "Note=pos_uncertainty:", tpsoritbl$Positional.Uncertainty, ",",
                             "common_name:", tpsoritbl$Common.Name, ",",
                             "old_name:", tpsoritbl$Old.Name, ",",
                             "pvalue:", tpsoritbl$p.Value)) %>% 
  mutate(acc = case_when(acc == "chr" ~ "NC_002607.1",
                         acc == "plasmid_pNRC100" ~ "NC_001869.1",
                         acc == "plasmid_pNRC200" ~ "NC_002608.1"))

# writing positive strand file
write_tsv(x = tpsgff %>% filter(strand == "+"),
          col_names = F,
          file = "data/TPS_gff_fwd.gff3")

# writing negative strand file
write_tsv(x = tpsgff %>% filter(strand == "-"),
          col_names = F,
          file = "data/TPS_gff_rev.gff3")
