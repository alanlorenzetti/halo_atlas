# alorenzetti 202007

# description ######

# this script will take data for Hundt et al. 2007
# on the half lives of halo transcripts

# in Hundt et al. 2007 the strain was Hsal NRC1
# but the loci names are from Hsal R1

# # source: https://www.omicsdi.org/dataset/arrayexpress-repository/E-MEXP-1088
# processed data

# importing data ####
# it is gonna throw a warning since there are columns with
# repetitive names
hldf = read_delim(file="data/E-MEXP-1088-processed-data-1369116981.txt",
                  delim ="\t", skip = 1)

# seems like there are two fields corresponding to hl
# parsing them, ND to NA and negative to NA
hl = hldf %>%
  dplyr::select(`Reporter REF`, `GenePix:half-life (min)...5`, `GenePix:half-life (min)...10`)
names(hl) = c("ref", "HL1", "HL2")

hl = hl %>%
  mutate(HL1 = na_if(HL1, "ND"), HL2 = na_if(HL2, "ND")) %>% 
  mutate(HL1 = as.numeric(HL1), HL2 = as.numeric(HL2))
hl$HL1 = replace(hl$HL1, which(hl$HL1 < 0), NA)
hl$HL2 = replace(hl$HL2, which(hl$HL2 < 0), NA)

hl = hl %>% drop_na()

# hl1 and hl2 should be a mean of both values
hl$mean = rowMeans(hl[,2:3])

# I will drop every entry holding a suffix after OEXXXX[FR] since
# I don't know what they mean
# I am also computing a mean of repetitive refs
hl = hl %>%
  filter(str_detect(ref, "^OE")) %>% 
  filter(str_detect(ref, "for|_e|orf|fehlt", negate = T)) %>% 
  dplyr::select(ref, mean) %>% 
  group_by(ref) %>% 
  summarise(mean = mean(mean))

# we need an approach to remove outliers
# this is not optimal but
# I will remove all not included within
# two standard deviations
m = hl$mean %>% mean()
sd = sd(hl$mean)
hl = hl %>%
  filter(mean >= m - (2*sd) & mean <= m + (2*sd))

# using a dictionary of locus tags to unify R1 (Hundt et al 2007)
# locus tags with NRC1 data
nrtxsep = read_tsv("https://alanlorenzetti.github.io/halo_nr_tx/data/dictionary.tsv") %>% 
  separate_rows(locus_tag, sep = ",")

# some of the genes in the Hundt et al. 2007 have
# different locus_tags but the same sequence
# that is common in halo
# therefore, two different OE locus_tags can have
# a same VNG locus_tag in our redundant transcriptome
# for those cases, I will take the mean half life
halfLives = left_join(nrtxsep, hl, by = c("locus_tag" = "ref")) %>% 
  filter(!is.na(mean)) %>% 
  dplyr::select(locus_tag = representative, mean) %>% 
  group_by(locus_tag) %>% 
  summarise(HL = mean(mean)) %>% 
  ungroup()

# I will saturate values higher than 95 percentile
# which is ~ 24 min
# halfLives %>% ggplot(aes(x=HL)) + geom_histogram()
sat95 = quantile(halfLives$HL, probs = 0.95)
halfLives = halfLives %>% 
  mutate(HL = case_when(HL > sat95 ~ sat95,
                   TRUE ~ as.numeric(HL)))
