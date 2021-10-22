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
estCounts2647 = read_tsv(file = "data/2647_nonTrim17_tableEstCounts.tsv") %>% 
  select(-c(length, eff_length))
tpm2647 = read_tsv(file = "data/2647_nonTrim17_tableTpm.tsv") %>% 
  select(-c(length, eff_length))

# better headers
colnames(estCounts2647)[2:5] = c("d2647_A",
                                 "d2647_B",
                                 "dura3_A",
                                 "dura3_B")

colnames(tpm2647)[2:5] = c("d2647_A",
                           "d2647_B",
                           "dura3_A",
                           "dura3_B")

# # adding a pseudocount to tpms
# computing log2foldchange using
# just the mean of tpms causes an
# unreallistic distribution
# therefore I will use DESeq2
# lfc2647 = tpm2647 %>%
#   mutate(across(.cols = where(is.numeric),
#                 .fns = ~ .x + 1)) %>% 
#   rowwise() %>% 
#   mutate(lfc2647 = log2(mean(d2647_A, d2647_B) / mean(dura3_A, dura3_B))) %>% 
#   select(target_id, lfc2647)
# 
# # checking lfc distribution
# lfc2647 %>%
#   ggplot(aes(x = lfc2647)) +
#   geom_density()

# processing raw counts with DESeq2 ####
# defining sample properties
samples = colnames(estCounts2647)[-1]
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

assay = estCounts2647[,-1] %>% 
  as.matrix() %>% 
  round(digits = 0)

# creating se object based on kallisto est counts
SE2647 = SummarizedExperiment(assay = list(counts=assay),
                          rowData = estCounts2647[,1],
                          colData = colData)
rownames(SE2647) = estCounts2647$target_id

dds2647 = SE2647
dds2647 = DESeqDataSet(dds2647,
                       design = ~ strain)

# removing genes with 0 counts across all samples
dds2647 = dds2647[rowSums(counts(dds2647)) > 1, ]

# generating model and getting the appropriate contrast
dds2647 = DESeq(dds2647)
res2647 = results(dds2647,
                  contrast = c("strain", "d2647", "dura3"),
                  alpha = padjthreshold) %>% 
  as_tibble(rownames = "target_id")

# checking lfc distribution
res2647 %>%
  ggplot(aes(x = log2FoldChange)) +
  geom_density()

# declaring final object containing log2 fold changes
# for each one of the genes
lfc2647 = res2647 %>% 
  select(representative = target_id,
         lfc2647 = log2FoldChange)
