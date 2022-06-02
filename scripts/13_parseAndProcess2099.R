# alorenzetti 202109

# description ####
# this script will get
# and parse RNase 2099 knockout
# microarray data in order
# to perform differential
# expression analysis

# part of the code was adapted
# from GEO2R platform from GEO@NCBI
# https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE45988

# downloading and saving data from GEO
# load series and platform data from GEO
if(!file.exists("data/GSE45988.RData")){
  gset = getGEO("GSE45988",
                GSEMatrix =TRUE,
                AnnotGPL=FALSE)
  
  save(gset, file = "data/GSE45988.RData")

}else{
  load("data/GSE45988.RData")
}

# taking object from list
gset = gset[[1]]

# since this is a tiling array experiment
# first I've got to find the genes
# intersecting each one of the probes

# adjusting data from gset to transform into
# a gr object; also adjusting probes
# spanning the end and beginning of a replicon
gsetfeat = fData(gset) %>%
  as_tibble() %>% 
  mutate(strand = case_when(str_detect(string = ID, pattern = "_fwd_") ~ "+",
                            TRUE ~ "-"),
         RANGE_START = RANGE_START %>% as.numeric(),
         RANGE_END = RANGE_END %>% as.numeric(),
         RANGE_END = case_when(RANGE_START > RANGE_END ~ RANGE_START,
                               TRUE ~ RANGE_END))

gsetgr = GRanges(seqnames = paste0(gsetfeat$RANGE_GB, ".1"),
                 ranges = IRanges(start = gsetfeat$RANGE_START,
                                  end = gsetfeat$RANGE_END),
                 strand = gsetfeat$strand)

names(gsetgr) = featureNames(gset)
  
# finding overlaps
ovlpsRes = GenomicRanges::findOverlaps(gsetgr, pfeiCDSgr,
                                       ignore.strand = F,
                                       type = "within") %>% 
  as_tibble()

ovlps = tibble(locus_tag = names(pfeiCDSgr)[ovlpsRes$subjectHits],
               probes = names(gsetgr)[ovlpsRes$queryHits])

# use all tiling array time points or
# only some of them (mid and late exponential phases)
entire="no"

# taking matrix from gset object
# and removing log10 transformation
gsetM = 10^exprs(gset) %>%
  as_tibble()

gsetM$ID = fData(gset)$ID

# joining with locus_tags names
# joining with non redundant transcriptome
# summarising using mean of probes
# and performing log2 transformation
gsetM = left_join(gsetM, ovlps,
                  by = c("ID" = "probes")) %>% 
  drop_na() %>% 
  left_join(., nrtxsep,
            by = "locus_tag") %>% 
  drop_na() %>% 
  select(-c(ID, product, locus_tag)) %>% 
  group_by(representative) %>% 
  summarise(across(.cols = starts_with("GSM"),
                   .fns = ~ log2(mean(.x))))

# creating a new ExpressionSet obj
newgset = gsetM[,-1] %>% as.matrix()
rownames(newgset) = gsetM$representative

gset = ExpressionSet(assayData = newgset,
                     phenoData = phenoData(gset))

# writing conditions
# with all
gset$strain = factor(x = c(rep("ura3",4), rep("2099",4)), levels = c("ura3", "2099"))
gset$timepoint = factor(x = paste0("TP", c(1:4, 1:4)), levels = c("TP1", "TP2", "TP3", "TP4"))

if(entire != "yes"){
  # subsetting matrix to include only
  # mid and late exponential phases
  gset = gset[,c(1,2,5,6)]
  
  gset$strain = factor(x = c(rep("ura3",2), rep("2099",2)), levels = c("ura3", "2099"))
  gset$timepoint = factor(x = paste0("TP", c(1:2, 1:2)), levels = c("TP1", "TP2"))
}

# defining design
design = model.matrix(~ strain + timepoint + 0, gset)

# fitting model
fit = lmFit(gset, design)
cont.matrix = makeContrasts(contrasts="strain2099-strainura3", levels=design)

# recalculating model coefficients
fit2 = contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 = eBayes(fit2, 0.01)
tT = topTable(fit2, number = dim(gset)[1])
tT$representative = rownames(tT)

# arranging table
tT = subset(tT, select=c("representative","adj.P.Val","P.Value","t","B","logFC"))
res2099 = tT %>% as_tibble()

# creating final obj
lfc2099 = res2099 %>% 
  arrange(representative) %>% 
  select(representative,
         lfc2099 = logFC)
