# alorenzetti 202007

# description #####
# this script will take deAlmeida et al. 2019
# annotation of antisense RNAs (parsed from Table S4)
# and find which genes of pfeiffer et al. 2019 intersect
# with them. this is part of functional
# categorization step

# reading antisense gff
asrnas = rtracklayer::import("data/Hsalinarum-846asRNAs-deAlmeida2019.gff3")
genes = rtracklayer::import("data/Hsalinarum-gene-annotation-pfeiffer2019.gff3")
genes = genes[genes$type == "gene",]

# inverting asrnas strand to find intersections
strand(asrnas) = invertStrand(strand(asrnas))

# checking intersections
geneswasrnas = subsetByOverlaps(genes, asrnas, ignore.strand=F) %>%
  as_tibble() %>%
  dplyr::select(locus_tag) %>%
  mutate(asRNA = "yes")

# if at least one of the transcript groups
# has an antisense, it is going to be asRNA == yes
geneswasrnas = left_join(x = nrtxsep,
                         y = geneswasrnas,
                         by = "locus_tag") %>% 
  select(-product) %>% 
  mutate(asRNA = case_when(asRNA == "yes" ~ TRUE,
                           TRUE ~ FALSE)) %>% 
  group_by(representative) %>% 
  summarise(asRNA = sum(asRNA)) %>% 
  ungroup() %>% 
  mutate(asRNA = case_when(asRNA >= 1 ~ "yes",
                           TRUE ~ "no"))
