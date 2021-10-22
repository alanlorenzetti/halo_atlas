# alorenzetti 20201103

# description ####
# this script will parse
# tf binding site data
# coming from ChIP-seq.ggb file
# provided by rvencio

# setting up promoter regions from CDS start #####
upstream = 150
downstream = 150

# loading and parsing tf file #####
tfggb = read_delim(file = "data/ChIP-seq.ggb",
                   delim = "\t")

# adjusting replicon names
tfggb = tfggb %>%
  mutate(Sequence = case_when(Sequence == "chr" ~ "NC_002607.1",
                              Sequence == "plasmid_pNRC100" ~ "NC_001869.1",
                              Sequence == "plasmid_pNRC200" ~ "NC_002608.1",
                              TRUE ~ as.character(Sequence)))

# creating granges object
tfgr = GRanges(seqnames = tfggb$Sequence,
               ranges = IRanges(start = tfggb$Start,
                                end = tfggb$End),
               strand = tfggb$Strand)

names(tfgr) = tfggb$Name %>%
  str_replace(., "_Chip-seq", "")
tfgr$name = names(tfgr)

# creating CDS annotation object from pfeiffer et al 2019 annot CDS ####
# using pfeiffer et al 2019 annotation
pfeiAnnotCDS = rtracklayer::import("data/Hsalinarum-gene-annotation-pfeiffer2019.gff3") %>% 
  as_tibble() %>% 
  filter(type == "CDS")

# creating granges obj
pfeiCDSgr = GRanges(seqnames = pfeiAnnotCDS$seqnames,
                    IRanges(start = pfeiAnnotCDS$start,
                            end = pfeiAnnotCDS$end),
                    strand = pfeiAnnotCDS$strand)

names(pfeiCDSgr) = pfeiAnnotCDS$locus_tag
pfeiCDSgr$locus_tag = pfeiAnnotCDS$locus_tag

# extending start position x bp upstream and y bp downstream
# see first section for hard coded variables upstream and downstream
pfeiCDSpromoterGr = GenomicRanges::promoters(pfeiCDSgr,
                                             upstream = upstream,
                                             downstream = downstream)

# finding overlaps ####
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

colnames(ovlps)[-1] = paste0("ChIPSeq_", colnames(ovlps)[-1])

# copying object
chipSeqTFs = ovlps

# if at least one of the transcript groups
# has a ChIP-Seq binding site
# it is going to be == yes
chipSeqTFs = left_join(x = nrtxsep,
                       y = chipSeqTFs,
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
