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

# also it is going to
# parse TOPCONS webserver result
# generated for CDS (pfei_2019_cds.fa)
# annotated by Pfeiffer et al 2019.

# getting started ####
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

# checking the proteins
# nrtx %>%
#   filter(representative %in%
#            ptgsAbund$union$prot_non_mrna_top[ptgsAbund$union$prot_non_mrna_top %in%
#                                                membProtsFinal]) %>% 
#   pull(representative)

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

# checking the proteins
undetectedRemovDF = nrtx %>%
  filter(representative %in%
           ptgsAbund$union$prot_non_mrna_top[ptgsAbund$union$prot_non_mrna_top %in%
                                               c(membProtsFinal, topsconMemb$two)])

# the following proteins should be removed
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
           
# updating our final dataframe
undetectedRemovDF = undetectedRemovDF %>% 
  filter(!representative %in% membFP)
