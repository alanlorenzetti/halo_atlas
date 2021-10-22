# alorenzetti 202008

# description ######
# this script will take the proteomics
# files output by OneOmics software
# and perform wrangling and filtering
# in this version
# we are using OneOmics output
# combining all the biological
# replicates, therefore, only one
# input table

# thresholds ####
padj = 0.05
pthr = padj
log2fcthreshold = 0.75 # Ulli said this is recommended by Sciex
borderlinezero = 0.25

# loading files #####
oneomics = read_xlsx("data/20200729_HL_VNGcplte_3BRs_RTCal_heatmap_0conf_0reproducibility_all2070proteins_s.xlsx") %>% 
  dplyr::select(-secondary_name)

# parsing colnames
colnames(oneomics) = colnames(oneomics)
colnames(oneomics) = gsub(" ", "_", colnames(oneomics))
colnames(oneomics) = gsub("log2_sfc", "lfc", colnames(oneomics))
colnames(oneomics) = sub("_vs._TP1", "", colnames(oneomics))
colnames(oneomics)[-1] = sub("(.*)_(.*)", "\\2_\\1", colnames(oneomics)[-1])
colnames(oneomics)[1] = "locus_tag"

# adding BH adjusted pval
oneomics = oneomics %>% 
  mutate(padj_TP2 = p.adjust(pval_TP2, method = "BH"),
         padj_TP3 = p.adjust(pval_TP3, method = "BH"),
         padj_TP4 = p.adjust(pval_TP4, method = "BH"))

# filtering based on
# c (confidence >= 0.75 )
# mlr >= 0.2
cthr = 0.75
mlrthr = 0.2
oneomics = oneomics %>% 
  filter((c_TP2 >= cthr & mlr_TP2 >= mlrthr) |
           c_TP3 >= cthr & mlr_TP3 >= mlrthr |
           c_TP4 >= cthr & mlr_TP4 >= mlrthr)

# pivoting object
oneomicsLong = oneomics %>% 
  pivot_longer(cols = contains("TP"),
               names_to = c("type", "timepoint"),
               names_sep = "_")

# filtering out attributes that are not
# going to be used anymore
oneomicsWide = oneomics %>% 
  dplyr::select(locus_tag, contains("lfc"), contains("padj"))

colnames(oneomicsWide)[-1] = sub("lfc_", "mean_lfc_protein_lysate_", colnames(oneomicsWide)[-1])
colnames(oneomicsWide)[-1] = sub("padj_", "mean_padj_protein_lysate_", colnames(oneomicsWide)[-1])

# adjusting locus_tags
# according to dictProd
# "VNG5199H"  "VNG0606G"  "VNG0779C" "VNG0780H" 
# don't have a matching sequence in pfeiLocusTag
# those are going to be removed
oneomicsWide = left_join(oneomicsWide,
                         nrtxsep %>% dplyr::select(-product),
                         by = "locus_tag")

remove = oneomicsWide$locus_tag[is.na(oneomicsWide$representative) %>% which()]
oneomicsWide = oneomicsWide[!(oneomicsWide$locus_tag %in% remove),] %>% 
  mutate(locus_tag = representative) %>% 
  dplyr::select(-representative)
  
# for retrocompatibility reasons
# I will recreate the object names tpivstp1
# it is going to be similar
# but not the same
tpivstp1 = list()
for(i in paste0("TP", 2:4)){
  tpivstp1[[i]] = oneomicsWide %>% 
    dplyr::select(locus_tag,
           matches(i)) %>% 
    rename_with(.fn = ~ str_replace(string = .x, pattern = "_TP.*$", replacement = ""),
                .cols = -locus_tag) %>% 
    mutate(sigProt = case_when(abs(mean_lfc_protein_lysate) >= log2fcthreshold &
                                 mean_padj_protein_lysate < pthr ~ "yes",
                               TRUE ~ "no"),
           borderLineZeroProt = case_when(abs(mean_lfc_protein_lysate) < borderlinezero ~ "yes",
                                          TRUE ~ "no"))
}

# exploratory charts ####
# plotting densities
oneomicsLong %>%
  filter(type == "lfc") %>% 
  ggplot(aes(x=value, color = timepoint)) +
  geom_density() +
  xlab("log2(Fold Change)")

# plotting heatmap
Heatmap(oneomics %>% dplyr::select(contains("lfc")) %>% as.matrix())
