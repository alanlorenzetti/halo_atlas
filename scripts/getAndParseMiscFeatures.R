# alorenzetti 202008

# description ####
# this script will get/load
# a few additional misc features
# LSm binding sites
# UTR info
# GC content

# processing starts ####
# getting and parsing LSm binding info
# here I will load the table, remove the asRNAs
# and lookup into genes and 5UTRs
# if LSm binding sites are detected
# in genes or 5UTRs, sense or antisense,
# that will be taken into consideration
lsmGenes = read_tsv(file = "data/interactionListNRTX.tsv") %>% 
  select(name = representative,
         lsmSense = LSmInteraction,
         lsmAntiSense = LSmInteractionAS) %>% 
  mutate(lsmSense = case_when(lsmSense == "Sim" ~ "yes",
                              TRUE ~ "no"),
         lsmAntiSense = case_when(lsmAntiSense == "Sim" ~ "yes",
                                  TRUE ~ "no"))

# # getting and parsing info about UTRs
# utr = read_delim(file = "data/5UTR.txt",
#                  delim = "\t",
#                  col_names = T) %>% 
#   dplyr::select(-length)
# 
# # getting and parsing info about UTR mfe
# utrmfe = read_delim(file = "data/5UTRplus18ntdownstream_mfe.txt",
#                     delim = "\t",
#                     col_names = T) %>% 
#   dplyr::select(-seq,-dotbracket)

# getting GC content info for each gene
GCcontent = tibble(locus_tag = names(nrtxseqs) %>%
                     sub("\\|.*$", "", .),
                   GC = nrtxseqs %>%
                     letterFrequency(letters = "GC", as.prob = T) %>%
                     as.numeric(),
                   GCdev = GC - mean(GC))

# parsing insertion sequence information
isinfo = read_tsv(file = "https://alanlorenzetti.github.io/halo_nr_tx/data/is_info.tsv") %>% 
  select(representative,
         ISFamily)
