# alorenzetti 20200520

# description ######

# this script will generate a dictionary of locus_tag
# using sequence similarity search (rbbh)
# the goal is to convert names from original halo genome annotation 
# to the nonredundant transcriptome (part of runKallisto scripts)
# derived from the third party annotation effort by pfeiffer et al 2019
# please check https://github.com/alanlorenzetti/runKallisto

# following will be done
# 1. get and parse ncbi genbank annotation for halo (also R1 annotation)
# 2. get and parse ncbi third party annotation done by pfeiffer2019 
# 3. extract sequences for both of them
# 4. run blast to generate the dictionary

# pfeifer2019 
# https://www.ncbi.nlm.nih.gov/nuccore/BK010829
# https://www.ncbi.nlm.nih.gov/nuccore/BK010830
# https://www.ncbi.nlm.nih.gov/nuccore/BK010831

# Loading libs ################################################ 
source("scripts/loadingLibs.R")

# Loading files ################################################ 

# our in-house non redundant transcriptome is being
# used here
pfeiFile = "data/pfeiGen.fa"

# path to protein database used by proteomics analysis
protFile = "data/Halobacterium-20080205_VNG_cRAP_TargDecoy_plusRT.fasta"

# setting download flag in case we need to download
# R1 annoation
download = "n"

# downloading genome and annotation
# for R1 strain
if(download == "y"){
  genoUrl = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/069/025/GCF_000069025.1_ASM6902v1/GCF_000069025.1_ASM6902v1_genomic.fna.gz"
  annotUrl = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/069/025/GCF_000069025.1_ASM6902v1/GCF_000069025.1_ASM6902v1_genomic.gff.gz"
  
  system2(command = "curl", args = paste("-u anonymous:pw", genoUrl, "| gunzip"),
          stdout = "data/HsalinarumR1.fa")
  system2(command = "curl", args = paste("-u anonymous:pw", annotUrl, "| gunzip"),
          stdout = "data/HsalinarumR1.gff")
}

# reading the protein database used in proteomics analysis
protSeqsAA = readAAStringSet(protFile)

# Processing starts here ################################################ 
# removing all decoys and prot seqs
# keeping only halo ORFs
protSeqsAA = protSeqsAA[1:2646]

# adjusting names and writing sequences
names(protSeqsAA) = sub(" .*$", "", names(protSeqsAA), perl = T)
writeXStringSet(protSeqsAA, "data/protGen.fa", format = "fasta")

# dealing with Pfeiffer 2019 annotation ################################################ 
# reading non redundant transcriptome file
pfeiSeqs = readDNAStringSet(filepath = pfeiFile, format = "fasta", use.names = T)

# removing ncRNAs (tRNAs and rRNAs)
pfeiSeqs = pfeiSeqs[!str_detect(string = names(pfeiSeqs), pattern = "VNG_t|VNG_r|VNG_s")]

# converting to aa based on 11 table and writing
arc = getGeneticCode(id_or_name2 = "11")
pfeiSeqsAA = translate(pfeiSeqs, genetic.code = arc)
writeXStringSet(pfeiSeqsAA, "data/pfeiGenAA.fa", format = "fasta")

# dealing with R1 strain annotation ################################################ 
# reading R1 annotation
r1annot = rtracklayer::import("data/HsalinarumR1.gff")
r1genes = r1annot[r1annot$type == "gene",]

# keeping only protein coding genes
r1genes = r1genes[r1genes$gene_biotype == "protein_coding",]

# droping r1genes without old_locus_tags
# they will not be useful to get the halflives from Hundt et al. 2007
r1genes = r1genes[!is.na(r1genes$old_locus_tag),]

# loading R1 genome
r1geno = readDNAStringSet("data/HsalinarumR1.fa")
names(r1geno) = sub(" .*$", "", names(r1geno))

# I tried to extract the sequences
# but it was throwing and error
# the problem is: gene annotation has 
# a flaw on gene OE_RS14455
# its limits are out of possible sequence
# range; it is going to be removed
r1genes = r1genes[!r1genes$locus_tag == "OE_RS14455",]

# extracting and parsing sequence names
r1Seqs = BSgenome::getSeq(r1geno, r1genes)
names(r1Seqs) = sub(",.*$", "", r1genes$old_locus_tag)

# converting to aas
r1SeqsAA = translate(r1Seqs, genetic.code = arc)
writeXStringSet(r1SeqsAA, "data/R1GenAA.fa", format = "fasta")

# running RBBH analysis ################################################ 
# RBBH between datasets
if(!file.exists("data/res.RData")){
  res = blast_best(query_file = "data/protGen.fa",
                   subject_file = "data/pfeiGenAA.fa",
                   seq_type = "protein")
  save(res, file = "data/res.RData")
} else {
  load("data/res.RData")
}

if(!file.exists("data/res2.RData")){
  res2 = blast_best(query_file = "data/R1GenAA.fa",
                    subject_file = "data/pfeiGenAA.fa",
                    seq_type = "protein")
  save(res2, file = "data/res2.RData")
} else {
  load("data/res2.RData")
}

# allowing correspondence only if
# identity >= .98 and subsetting
# locus_tag cols
resFil = res %>%
  filter(perc_identity >= 98) %>% 
  dplyr::select(query_id, subject_id) %>% 
  ungroup()

# resFil will be transformed further
# but res2Fil is good to go after the
# following step
# it is going to be used by parseHalfLives.R
res2Fil = res2 %>% filter(perc_identity >= 98) %>% 
  dplyr::select(query_id, subject_id) %>% 
  mutate(query_id = sub("OE_", "OE", query_id)) %>% 
  ungroup()

dictR1 = res2Fil %>% 
  mutate(subject_id = sub("\\|.*$","", subject_id))
colnames(dictR1) = c("R1locus_tag", "pfeiLocus_tag")

# finalizing ################################################ 
dict = resFil %>% 
  mutate(query_id = sub("\\|.*$", "", query_id),
         subject_id = sub("\\|.*$", "", subject_id)) %>% 
  mutate(query_id = sub("-.*$", "", query_id)) %>% 
  mutate(query_id = sub("_", "", query_id))

# assigning products to dict
# loading pfeiffer annotation file
# loading pfeiffer2019 annotation
pfeiFile = "data/Hsalinarum-gene-annotation-pfeiffer2019.gff3"
pfei = rtracklayer::import(pfeiFile) %>% as_tibble()

############PFEI
# parsing annotation
pfeiGenes = pfei %>% filter(type == "gene" | type == "pseudogene")
pfeiCDS = pfei %>%
  filter(type == "CDS") %>% 
  dplyr::select(Parent, product) %>% 
  mutate(Parent = as.character(Parent))
pfeiGenes = left_join(pfeiGenes, pfeiCDS, by = c("ID" = "Parent")) %>% 
  dplyr::select(locus_tag, product.y)

# unifying products with dict
dictProd = left_join(x=dict, y=pfeiGenes, by=c("subject_id" = "locus_tag"))

# finding transposases
# filt = grepl(pattern = "insertion element|ISH|transposase", x = dictProd$product.y, perl = T)
# dictTransp = dictProd[filt,]