# alorenzetti 20220303

# description ####
# this script will parse
# membrane proteins detected
# by Goo et al 2003 (table I - membrane fraction)
# and Klein et al 2005
# (tables I and II; only normal and reliable confidence)
# in order to eliminate false
# negatives not detected
# by our proteomics survey

# also we will parse proteins
# not included in the SWATH lib
# according to Kusebauch et al. (in prep)

# also it is going to
# parse TOPCONS webserver result
# generated for CDS (pfei_2019_cds.fa)
# annotated by Pfeiffer et al 2019.

# getting started ####
# membrane proteins ####
# loading data
membProtsGoo = read_csv(file = "data/tabula-goo2003onlytables.csv") %>% 
  select(-c(`Functional Category`, Probability, Protein, `Putative Function`)) %>% 
  rename(locus_tag = "Gene ID")

membProtsKlein = read_csv(file = "data/tabula-klein2005onlytables.csv") %>% 
  filter(identification_level != "q") %>% 
  mutate(locus_tag = case_when(!str_detect(locus_tag, "^VNG") ~ locus_tag_r1,
                               TRUE ~ locus_tag)) %>% 
  relocate(locus_tag)

membProtsLT = c(membProtsGoo$locus_tag, membProtsKlein$locus_tag) %>% unique()

membProtsFinal = nrtxsep %>% 
  filter(locus_tag %in% membProtsLT) %>% 
  pull(representative) %>% 
  unique()

# loading topcons output
topsconRes = read_tsv(file = "data/pfei_2019_cds_topcons_finished_seqs.txt",
                      col_names = c("No.",
                                    "Length",
                                    "numTM",
                                    "SignalPeptide",
                                    "Source",
                                    "Runtime",
                                    "locus_tag",
                                    "FinishDate")) %>% 
  select(locus_tag,
         Length,
         numTM,
         SignalPeptide)

# selecting proteins with at least two TM
# domain since the ones with just one TM
# does not show satisfying delta G
# (0 would be ideal for TM regions).
# we also observed that proteins with
# just one TM lack confidence in predictions
# (confidence < 0.9).
topsconMemb = list()
# topsconMemb$one = topsconRes %>%
#   filter(numTM > 0) %>%
#   pull(locus_tag)

topsconMemb$two = topsconRes %>%
  filter(numTM > 1) %>% 
  pull(locus_tag)

# since they were likely false positives
# for membrane proteins according to our
# manual inspection
membFP = c("VNG_0551G",
           "VNG_1801G",
           "VNG_0052H",
           "VNG_0059H",
           "VNG_0420H",
           "VNG_0767H",
           "VNG_1653H",
           "VNG_2089H")

memb_to_remove = c(membProtsFinal, topsconMemb$two)
memb_to_remove = memb_to_remove[!memb_to_remove %in% membFP]

# proteins not detectable by SWATH ####
# we need to check what was included
# in SWATH library to learn what 
# was not represented because of technical
# reasons and what was not annotated
# by the time of SWATH-MS procedure

# reading protein search db for swath
# swath_search_db = readAAStringSet(filepath = "data/Halobacterium-20080205_VNG_cRAP_TargDecoy_plusRT.fasta")
# swath_search_db_names = names(swath_search_db)
# swath_search_db_names = swath_search_db_names %>% 
#   str_detect(pattern = "^VNG") %>% 
#   swath_search_db_names[.] %>% 
#   str_replace(pattern = "^(.*?) .*$", "\\1")

# these 83 proteins were not represented
# in the SWATH assay library so
# they cannot be detected
unrepresentedTib = tibble(locus_tag = c("VNG0076H","VNG0205H","VNG0261H","VNG0323H","VNG0420H",
                                         "VNG0506H","VNG0518H","VNG0548C","VNG0568C","VNG0613H",
                                         "VNG0641C","VNG0642C","VNG0652H","VNG0656H","VNG0659H",
                                         "VNG0741H","VNG0762H","VNG0772H","VNG0863H",
                                         "VNG0934H","VNG0945H","VNG0988H","VNG0990H",
                                         "VNG1200H","VNG1229H","VNG1263C","VNG1365C",
                                         "VNG1402H","VNG1447H","VNG1479H","VNG1619H",
                                         "VNG1626C","VNG1630H","VNG1664H","VNG1678H",
                                         "VNG1801G","VNG1886a","VNG1904H","VNG1919H",
                                         "VNG1921H","VNG1927H","VNG1943H","VNG2081H",
                                         "VNG2091H","VNG2098H","VNG2129H","VNG2137G",
                                         "VNG2189H","VNG2236H","VNG2244H","VNG2264C",
                                         "VNG2290G","VNG2314H","VNG2353H","VNG2359G",
                                         "VNG2456Cm","VNG2477H","VNG2510H","VNG2556H",
                                         "VNG2581H","VNG2594C","VNG2673H","VNG2674H",
                                         "VNG5054H","VNG5073H","VNG5083H","VNG5112H",
                                         "VNG5175H","VNG5214H","VNG5224H","VNG5243H",
                                         "VNG6052H","VNG6070H","VNG6080H","VNG6109H",
                                         "VNG6159H","VNG6209H","VNG6276H","VNG6286H",
                                         "VNG6330H","VNG6446H","VNG6456H","VNG6474H"))

unrepresented = left_join(x = unrepresentedTib,
                          y = nrtxsep,
                          by = "locus_tag") %>% 
  filter(!is.na(representative)) %>% 
  group_by(representative) %>% 
  summarise(locus_tag = paste0(locus_tag, collapse = ","),
            product = paste0(product, collapse = ","))

# checking what is in nrtx but not in SWATH
nrtx_not_in_swath_assay = nrtx %>%
  filter(str_detect(string = representative, pattern = "^VNG_[0-9]{4}")) %>% 
  filter(!str_detect(string = locus_tag, pattern = "VNG[0-9]{4}[A-Z]{1,2}")) %>% 
  pull(representative)

nrtx_not_in_swath_assayTib = nrtx %>% 
  filter(representative %in% c(nrtx_not_in_swath_assay))

# finding proteins not producing peptides
# >= 7 aa and <= 30 aa

# for that
# we used rapid peptides generator
# rpg -i data/pfei_2019_cds.fa -c 4 -e 42 -f tsv > results/pfei_2019_cds_rpg.tsv
rpg_results = read_delim(file = "results/pfei_2019_cds_rpg.tsv") %>% 
  filter(Peptide_size >= 7 & Peptide_size <= 30)

rpg_valid_prots = left_join(x = nrtxsep,
                            y = rpg_results,
                            by = c("locus_tag" = "Original_header")) %>% 
  drop_na() %>% 
  pull(representative) %>% 
  unique()

rpg_not_valid_prots = nrtx %>%
  filter(str_detect(representative, "VNG_[0-9]{4}")) %>%
  filter(!representative %in% rpg_valid_prots) %>% 
  pull(representative)

# unifying all proteins that we can't know if are
# really undetected due to many reasons
unrepresentedFP = c(memb_to_remove,
                    unrepresented$representative,
                    nrtx_not_in_swath_assayTib$representative,
                    rpg_not_valid_prots) %>% 
  unique()
