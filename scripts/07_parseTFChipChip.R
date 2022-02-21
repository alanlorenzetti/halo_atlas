# alorenzetti 20201104

# description ####
# this script will parse
# tf binding site data
# coming from ChIP-chip.ggb file
# provided by rvencio

# loading and parsing tf file #####
tfggb = read_delim(file = "data/ChIP-chip.ggb",
                   delim = "\t")

# adjusting replicon names
tfggb = tfggb %>%
  mutate(Sequence = case_when(Sequence == "chr" ~ "NC_002607.1",
                              Sequence == "plasmid_pNRC100" ~ "NC_001869.1",
                              Sequence == "plasmid_pNRC200" ~ "NC_002608.1",
                              TRUE ~ as.character(Sequence)))

# filtering out non tbp and tfb
tfggb = tfggb %>% 
  filter(str_detect(Name, pattern = "^tb|^tf"))

# creating granges object
tfgr = GRanges(seqnames = tfggb$Sequence,
               ranges = IRanges(start = tfggb$Start,
                                end = tfggb$End),
               strand = tfggb$Strand)

names(tfgr) = tfggb$Name %>%
  str_replace(., "_ChIP-chip", "")
tfgr$name = names(tfgr)

# finding overlaps #####
# using promoter regions created by
# parseTFChipSeq.R script to 
# find overlapping features
ovlpsRes = GenomicRanges::findOverlaps(pfeiCDSpromoterGr, tfgr,
                                       ignore.strand = T) %>% 
  as_tibble()

ovlps = tibble(locus_tag = names(pfeiCDSpromoterGr)[ovlpsRes$queryHits],
               tf = names(tfgr)[ovlpsRes$subjectHits])

# pivoting dataframe
ovlps = ovlps %>%
  unique() %>% 
  pivot_wider(names_from = tf,
              values_from = tf) %>% 
  mutate(across(-locus_tag, .fns = ~ case_when(is.na(.x) ~ "no",
                                               TRUE ~ "yes")))

ovlps = ovlps %>% 
  dplyr::select(c("locus_tag", colnames(ovlps)[-1] %>% sort()))

colnames(ovlps)[-1] = paste0("ChIPChip_", colnames(ovlps)[-1])

# copying object
chipChipTFs = ovlps

# if at least one of the transcript groups
# has a ChIP-Seq binding site
# it is going to be == yes
chipChipTFs = left_join(x = nrtxsep,
                        y = chipChipTFs,
                        by = "locus_tag") %>% 
  select(-product) %>% 
  mutate(across(.cols = starts_with("Ch"),
                .fns = ~ case_when(.x == "yes" ~ TRUE,
                                   TRUE ~ FALSE))) %>% 
  group_by(representative) %>% 
  summarise(across(.cols = starts_with("Ch"),
                   .fns = ~ sum(.x))) %>% 
  ungroup() %>% 
  mutate(across(.cols = starts_with("Ch"),
                .fns = ~ case_when(.x >= 1 ~ "yes",
                                   TRUE ~ "no")))
