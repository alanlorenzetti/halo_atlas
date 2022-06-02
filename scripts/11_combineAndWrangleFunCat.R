# alorenzetti 202103

# description ####
# this script will take parsed dataframes
# generated previously and unify
# in a comprenhensive dataframe

# starting processing ####
# COG categorization from now on will
# be collected from a companion repository:
# the non-redundant transcriptome of Hsalinarum NRC-1.
# https://alanlorenzetti.github.io/halo_nr_tx/
COG = read_tsv(file = "https://alanlorenzetti.github.io/halo_nr_tx/data/cog.tsv")
  
# there are just thirteen genes having more than one distinct COG
# COG %>%
#   filter(str_detect(string = cog_category,
#                     pattern = "\\|"))

# I am arbitrarily accepting only the first category for each one of those
COG = COG %>%
  mutate(cog_category = str_replace(string = cog_category,
                                    pattern = "\\|.*$",
                                    replacement = ""))

# operons
# uniArCOGKEGGOper = left_join(uniArCOGKEGG, operons, by = "locus_tag")

# adding coginfo and alternative locus_tags
dictFunCat = left_join(nrtx, COG, by = "representative") %>% 
  rename(pfeiLocusTag = "representative")

# adding lsm data
dictFunCat = left_join(dictFunCat, smap1InteractionDF, by = c("pfeiLocusTag" = "representative"))

# adding UTR info
#dictFunCat = left_join(dictFunCat, utr, by = c("locus_tag" = "locus_tag"))

# adding asRNA info
dictFunCat = left_join(dictFunCat, geneswasrnas, by = c("pfeiLocusTag" = "representative"))

# adding utr mfe info
# wasn't created using nrtx
#dictFunCat = left_join(dictFunCat, utrmfe, by = c("locus_tag" = "locus_tag")) %>% 
#  rename(utrGC = GC,
#         utrmfe = mfe)

# adding GC content info
dictFunCat = left_join(dictFunCat, GCcontent, by = c("pfeiLocusTag" = "locus_tag"))

# adding halflives from Hundt et al. 2007
dictFunCat = left_join(dictFunCat, halfLives, by = c("pfeiLocusTag" = "locus_tag"))

# adding codon index usage
dictFunCat = left_join(dictFunCat, cai, by = c("pfeiLocusTag" = "locus_tag"))

# # adding transcription factor information obtained
# # from chip-seq
# dictFunCat = left_join(dictFunCat, chipSeqTFs, by = c("pfeiLocusTag" = "representative"))
# 
# # chip-chip experiments
# dictFunCat = left_join(dictFunCat, chipChipTFs, by = c("pfeiLocusTag" = "representative"))

# adding location from where tx come from (chromosome, plasmids)
dictFunCat = left_join(dictFunCat, location, by = c("pfeiLocusTag" = "representative"))

# adding TPS info
dictFunCat = left_join(dictFunCat, tps, by = c("pfeiLocusTag" = "representative"))

# # adding IRs info
# dictFunCat = left_join(dictFunCat, irs, by = c("pfeiLocusTag" = "representative")) %>% 
#   mutate(IRs = case_when(is.na(IRs) ~ 0,
#                          TRUE ~ as.numeric(IRs)))

# adding delta2647 fold changes
# dictFunCat = left_join(dictFunCat, lfc2647, by = c("pfeiLocusTag" = "representative"))

# adding delta1503 fold changes 
# dictFunCat = left_join(dictFunCat, lfc1503, by = c("pfeiLocusTag" = "representative"))

# adding delta2099 fold changes 
dictFunCat = left_join(dictFunCat, lfc2099, by = c("pfeiLocusTag" = "representative"))

# adjusting missing values for delta 2647, 1503, and 2099 fold changes
dictFunCat = dictFunCat %>% 
  mutate(
    # lfc2647 = case_when(is.na(lfc2647) ~ 0,
    #                     TRUE ~ lfc2647),
    # lfc1503 = case_when(is.na(lfc1503) ~ 0,
    #                     TRUE ~ lfc1503),
    lfc2099 = case_when(is.na(lfc2099) ~ 0,
                        TRUE ~ lfc2099))

# adding IS family information for representative CDS loci
# within insertion sequences
dictFunCat = left_join(dictFunCat, isinfo, by = c("pfeiLocusTag" = "representative"))

# adjusting missing values for a few vars
dictFunCat = dictFunCat %>% 
  mutate(cog_category = case_when(is.na(cog_category) ~ "Function unknown",
                                  TRUE ~ as.character(cog_category)))

# ISH2
# also, all transposases according to pfeiProduct will receive mobilome class
dictFunCat = dictFunCat %>%
  mutate(cog_category = case_when(str_detect(product, pattern = "ISH2|transposase") &
                                    str_detect(product, pattern = "nonfunc", negate = T) ~ "Mobilome: prophages, transposons",
                                  TRUE ~ as.character(cog_category)))

# adding information if a gene
# encodes a protein or not
dictFunCat = dictFunCat %>% 
  mutate(protein_coding = case_when(str_detect(pfeiLocusTag,
                                               pattern = "VNG_r|VNG_t|VNG_s",
                                               negate = T) ~ "yes",
                                    TRUE ~ "no"))
  

########## VISUALIZATION ###############
# funcat visualization
# lsm per arcog
# the results show lsm is proportionally
# more represented in mobilome and defense mechanism
# category; we have to be cautious, because there
# are repetitive elements within this categories

# # lsm per arcog
# ggplot(data = dictFunCat) +
#   geom_bar(aes(x=arCOG, fill=lsmSense), position="fill") +
#   coord_flip() +
#   ggtitle("LSm bound RNAs and arCOGs")
# 
# # antisense per arcog
# ggplot(data = dictFunCat) +
#   geom_bar(aes(x=arCOG, fill=asRNA), position="fill") +
#   coord_flip() +
#   ggtitle("Antisense RNAs and arCOGs")
# 
# # utrSize per arcog
# ggplot(data = dictFunCat) +
#   geom_bar(aes(x=arCOG, fill=utrSize), position="fill") +
#   coord_flip() +
#   ggtitle("5' UTR length and arCOGs")
# 
# # # utr mfe vs lsm
# # ggplot(data = dictFunCat, aes(y=mfe, x=lsmSense)) +
# #   geom_violin() +
# #   geom_beeswarm()
# 
# ggplot(data = dictFunCat, aes(x=mfe, color=lsmSense)) +
#   geom_density() +
#   ggtitle("5' UTR MFE distributions")
# 
# # # utr mfe vs GC
# # ggplot(data = dictFunCat, aes(y=GC, x=lsmSense)) +
# #   geom_boxplot() +
# #   geom_beeswarm()
# 
# ggplot(data = dictFunCat, aes(x=GC, color=lsmSense)) +
#   geom_density() +
#   ggtitle("5' UTR GC content distributions")
# 
# # checking half lives per arcog
# ggplot(data = dictFunCat, aes(y=HL, x=arCOG)) +
#   geom_boxplot() +
# #  geom_beeswarm() + 
#   coord_flip() +
#   ggtitle("Half Life distribution per category")
# 
# # comparing if distributions of GC are
# # different for lsm bound and unbound UTRs
# # we are going to compute empirical cumulative
# # distribution function for both of them and
# # check if they are similar
# # lsm bound
# lsmbound = dictFunCat[dictFunCat$lsmSense == "no","GC"] %>% drop_na() %>% unlist() %>% unname()
# lsmboundFun = ecdf(lsmbound)
# 
# lsmunbound = dictFunCat[dictFunCat$lsmSense == "yes","GC"] %>% drop_na() %>% unlist() %>% unname()
# lsmunboundFun = ecdf(lsmunbound)
# 
# # plotting
# plot(ecdf(lsmbound), xlim=range(c(lsmbound, lsmunbound)))
# plot(ecdf(lsmunbound), add=TRUE, lty="dashed")
# 
# # performing kolmogorov-smirnoff to compare
# # both functions
# # and we reject the null hypothesis
# ks.test(x = lsmbound, y = lsmunboundFun)
