# alorenzetti 20210617

# description ####
# this script will take the read counts
# of an RNA-Seq RNase null mutant
# experiment in order to integrate the results
# into the TLR compendium.
# it is going to process the counts
# and obtain simple fold changes as well

# starting ####
# setting thresholds
padjthreshold = 0.05

# reading files
estCounts1503 = read_tsv(file = "data/1503_nonTrim13_tableEstCounts.tsv") %>% 
  select(-c(length, eff_length))
tpm1503 = read_tsv(file = "data/1503_nonTrim13_tableTpm.tsv") %>% 
  select(-c(length, eff_length))

# better headers
colnames(estCounts1503)[2:5] = c("d1503_A",
                                 "d1503_B",
                                 "dura3_A",
                                 "dura3_B")

colnames(tpm1503)[2:5] = c("d1503_A",
                           "d1503_B",
                           "dura3_A",
                           "dura3_B")

# # adding a pseudocount to tpms
# computing log2foldchange
# lfc1503 = tpm1503 %>%
#   mutate(across(.cols = where(is.numeric),
#                 .fns = ~ .x + 1)) %>%
#   rowwise() %>%
#   mutate(lfc1503 = log2(mean(d1503_A, d1503_B) / mean(dura3_A, dura3_B))) %>%
#   select(target_id, lfc1503)
# 
# # checking lfc distribution
# lfc1503 %>%
#   ggplot(aes(x = lfc1503)) +
#   geom_density()

# processing raw counts with DESeq2 ####
# defining sample properties
samples = colnames(estCounts1503)[-1]
strain = samples %>% 
  str_replace(string = .,
              pattern = "_.*",
              replacement = "")
bioRep = samples %>% 
  str_replace(string = .,
              pattern = ".*_",
              replacement = "")

colData = data.frame(row.names = samples,
                     strain = strain,
                     replicate = bioRep)

assay = estCounts1503[,-1] %>% 
  as.matrix() %>% 
  round(digits = 0)

# creating se object based on kallisto est counts
SE1503 = SummarizedExperiment(assay = list(counts=assay),
                              rowData = estCounts1503[,1],
                              colData = colData)
rownames(SE1503) = estCounts1503$target_id

dds1503 = SE1503
dds1503 = DESeqDataSet(dds1503,
                       design = ~ strain)

# removing genes with 0 counts across all samples
dds1503 = dds1503[rowSums(counts(dds1503)) > 1, ]

# generating model and getting the appropriate contrast
dds1503 = DESeq(dds1503)
res1503 = results(dds1503,
                  contrast = c("strain", "d1503", "dura3"),
                  alpha = padjthreshold) %>% 
  as_tibble(rownames = "target_id")

# checking lfc distribution
res1503 %>%
  ggplot(aes(x = log2FoldChange)) +
  geom_density()

# declaring final object containing log2 fold changes
# for each one of the genes
lfc1503 = res1503 %>% 
  select(representative = target_id,
         lfc1503 = log2FoldChange)
