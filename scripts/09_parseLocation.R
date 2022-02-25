# alorenzetti 202108

# description ####
# this script will parse
# gene location considering
# chr and plasmids
# for that I am going to rely
# on an object created by
# parseTFChipSeq script

# getting started ####
loc = pfeiAnnotCDS %>% 
  mutate(seqnames = case_when(seqnames == "NC_002607.1" ~ "Chromosome",
                              seqnames == "NC_001869.1" ~ "pNRC100",
                              TRUE ~ "pNRC200")) %>% 
  select(seqnames,
         locus_tag)

# joining with nrtx dict
location = left_join(x = nrtxsep,
                     y = loc,
                     by = "locus_tag") %>% 
  filter(!is.na(seqnames)) %>% 
  select(-product) %>% 
  group_by(representative) %>% 
  summarise(seqnames = c(seqnames) %>% 
              unique() %>% 
              sort() %>% 
              paste0(collapse = ",")) %>% 
  separate_rows(seqnames) %>% 
  mutate(values = "yes") %>% 
  pivot_wider(names_from = seqnames,
              values_from = values,
              values_fill = "no")
  
