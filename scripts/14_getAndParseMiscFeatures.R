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
# lsmGenes = read_tsv(file = "data/interactionListNRTX.tsv") %>% 
#   select(name = representative,
#          lsmSense = LSmInteraction,
#          lsmAntiSense = LSmInteractionAS) %>% 
#   mutate(lsmSense = case_when(lsmSense == "Sim" ~ "yes",
#                               TRUE ~ "no"),
#          lsmAntiSense = case_when(lsmAntiSense == "Sim" ~ "yes",
#                                   TRUE ~ "no"))

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

# getting Tfb labels info
# based on Baliga Lab
# halo annotation resource
gtfInfo = list()

gtfInfo$tib = tibble(locus_tag = c("VNG5039G","VNG5142G",
                                   "VNG5163G","VNG2243G",
                                   "VNG5052G","VNG0315G",
                                   "VNG6438G","VNG0734G",
                                   "VNG6389G","VNG0254G",
                                   "VNG6351G","VNG0869G",
                                   "VNG2184G"),
                     label = c("TbpA","TbpC",
                               "TbpD","TbpE",
                               "TbpB","TfbF",
                               "TbpF","TfbB",
                               "TfbE","TfbG",
                               "TfbC","TfbD",
                               "TfbA")) %>% 
  left_join(x = .,
            y = nrtxsep,
            by = "locus_tag") %>% 
  arrange(label)

