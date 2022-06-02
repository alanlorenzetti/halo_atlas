# alorenzetti 20210211

# description ####
# this script will load files,
# get smap1 interaction results

# start analysis ####

# loading and parsing halo annotation
pfeiAnnot = rtracklayer::import("data/Hsalinarum-gene-annotation-pfeiffer2019.gff3")
pfeiGenes = subset(pfeiAnnot, type == "gene" | type == "pseudogene")
names(pfeiGenes) = pfeiGenes$locus_tag

# loading interaction tracks
int = list()

# loading interaction regions generated
# using biorep1 vs control
int[["biorep1"]] = c(
  rtracklayer::import("data/igv/br1biorep1-interaction-regions-entire-genome-fwd.gff3"),
  rtracklayer::import("data/igv/br1biorep1-interaction-regions-entire-genome-rev.gff3")
)

# using biorep2 vs control
int[["biorep2"]] = c(
  rtracklayer::import("data/igv/br2biorep2-interaction-regions-entire-genome-fwd.gff3"),
  rtracklayer::import("data/igv/br2biorep2-interaction-regions-entire-genome-rev.gff3")
) 

# unifying br1 and br2 tracks for fwd and rev
# that is going to be used in our browser app


# getting what genes intersect with interaction regions
overlaps = list()

# interaction regions of biorep1
overlaps[["biorep1"]] = pfeiGenes[
  findOverlaps(query = pfeiGenes, subject = int[["biorep1"]]) %>%
    as_tibble() %>% 
    select(queryHits) %>% 
    unique() %>% 
    unlist(use.names = F)
] %>% names()

# interaction regions of biorep2
overlaps[["biorep2"]] = pfeiGenes[
  findOverlaps(query = pfeiGenes, subject = int[["biorep2"]]) %>%
    as_tibble() %>% 
    select(queryHits) %>% 
    unique() %>% 
    unlist(use.names = F)
] %>% names()

# finding the union list of locus_tags
overlaps[["union"]] = c(overlaps[["biorep1"]],
                        overlaps[["biorep2"]]) %>% 
  unique()

# doing the same but trying to find
# interaction regions on the opposite strand
# interaction regions of biorep1 antisense
overlaps[["biorep1_as"]] = pfeiGenes[
  findOverlaps(query = pfeiGenes %>% invertStrand(), subject = int[["biorep1"]]) %>%
    as_tibble() %>% 
    select(queryHits) %>% 
    unique() %>% 
    unlist(use.names = F)
] %>% names()

# interaction regions of biorep2 antisense
overlaps[["biorep2_as"]] = pfeiGenes[
  findOverlaps(query = pfeiGenes %>% invertStrand(), subject = int[["biorep2"]]) %>%
    as_tibble() %>% 
    select(queryHits) %>% 
    unique() %>% 
    unlist(use.names = F)
] %>% names()

# finding the union list of locus_tags antisense
overlaps[["union_as"]] = c(overlaps[["biorep1_as"]],
                           overlaps[["biorep2_as"]]) %>% 
  unique()

# creating an interaction tibble and
# getting the representative locus_tag
# for each one of the locus_tags within the union list
df = list()
df[["union"]] = tibble(locus_tag = overlaps[["union"]],
                       LSmInteraction = "Sim")

df[["redundant_df"]] = nrtx %>%
  separate_rows(locus_tag, sep = ",") %>% 
  left_join(x = ., y = df[["union"]], by = "locus_tag") %>% 
  drop_na()

df[["only_representative"]] = df[["redundant_df"]] %>% 
  dplyr::select(-locus_tag) %>% 
  unique()

df[["locus_tags_representative"]] = df[["only_representative"]] %>%
  select(representative) %>%
  unlist(use.names = F)

# doing the same for antisense
df[["union_as"]] = tibble(locus_tag = overlaps[["union_as"]],
                          LSmInteractionAS = "Sim")

df[["redundant_df_as"]] = nrtx %>%
  separate_rows(locus_tag, sep = ",") %>% 
  left_join(x = ., y = df[["union_as"]], by = "locus_tag") %>% 
  drop_na()

df[["only_representative_as"]] = df[["redundant_df_as"]] %>% 
  dplyr::select(-locus_tag) %>% 
  unique()

df[["locus_tags_representative_as"]] = df[["only_representative_as"]] %>%
  select(representative) %>%
  unlist(use.names = F)

# nonredundant list per biorep
overlaps[["biorep1_nr"]] = tibble(locus_tag = overlaps[["biorep1"]],
                                  smap1_interaction = "yes") %>% 
  left_join(x = nrtxsep,
            y = .,
            by = "locus_tag") %>% 
  drop_na() %>% 
  select(-locus_tag) %>% 
  distinct() %>% 
  pull(representative)

overlaps[["biorep2_nr"]] = tibble(locus_tag = overlaps[["biorep2"]],
                                  smap1_interaction = "yes") %>% 
  left_join(x = nrtxsep,
            y = .,
            by = "locus_tag") %>% 
  drop_na() %>% 
  select(-locus_tag) %>% 
  distinct() %>% 
  pull(representative)

# venn diagram of bioreps gene-wise
# plot(eulerr::venn(combinations = list("BR1" = overlaps[["biorep1_nr"]],
#                                       "BR2" = overlaps[["biorep2_nr"]])))

# getting an alias for df$representative
smap1InteractionDF = nrtx %>% 
  left_join(x = .,
            y = df[["only_representative"]] %>% select(-product),
            by = "representative") %>% 
  mutate(LSmInteraction = case_when(is.na(LSmInteraction) ~ "Não",
                                    TRUE ~ as.character(LSmInteraction))) %>% 
  left_join(x = .,
            y = df[["only_representative_as"]] %>% select(-product),
            by = "representative") %>% 
  mutate(LSmInteractionAS = case_when(is.na(LSmInteractionAS) ~ "Não",
                                      TRUE ~ as.character(LSmInteractionAS)),
         geneClass = case_when(str_detect(product, "ISH|transposase") & 
                                 str_detect(product, "nonfunc", negate = T) ~ "Transposase",
                               TRUE ~ "Outra")) %>% 
  select(representative,
         lsmSense = LSmInteraction,
         lsmAntiSense = LSmInteractionAS) %>% 
  mutate(lsmSense = case_when(lsmSense == "Não" ~ "no",
                              lsmSense == "Sim" ~ "yes",
                              TRUE ~ lsmSense),
         lsmAntiSense = case_when(lsmAntiSense == "Não" ~ "no",
                                  lsmAntiSense == "Sim" ~ "yes",
                                  TRUE ~ lsmAntiSense))
