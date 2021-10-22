# alorenzetti 20210727

# description ####
# this script will run the tool
# non-b DNA motif search
# in order to get information
# about inverted repeats in halo genome

# the script will try
# to download and compile the required
# software. if not possible
# manual installation is required

# getting started ####
# this first part is only
# gonna be performed if
# the IR file does not exist yet
if(!file.exists("data/nonBpred_IR.gff")){
  # downloading software from github ####
  
  if(dir.exists("non-B_gfa")){
    system(command = "rm -r non-B_gfa")
  }
  
  if(dir.exists("data/non-B_gfa")){
    system(command = "rm -r data/non-B_gfa")
  }
  
  system(command = "git clone https://github.com/abcsFrederick/non-B_gfa")
  system(command = "mv non-B_gfa data")
  
  # compiling
  system(command = "cd data/non-B_gfa ; make ; cd -")
  
  # running software to predict IRs ####
  args = c("-seq data/Hsalinarum.fa",
           "-out data/nonBpred",
           "-skipAPR",
           "-skipSTR",
           "-skipDR",
           "-skipMR",
           "-skipGQ",
           "-skipZ",
           "-skipWGET")
  
  system2(command = "data/non-B_gfa/gfa",
          args = args,
          stdout = "data/nonBpred.log",
          stderr = "data/nonBpred.err")
}

# parsing results ####
# importing gff
irgr = rtracklayer::import("data/nonBpred_IR.gff")

# finding overlaps 
ovlpsRes = GenomicRanges::findOverlaps(pfeiCDSgr, irgr,
                                       ignore.strand = T) %>% 
  as_tibble()

ovlps = tibble(locus_tag = names(pfeiCDSgr)[ovlpsRes$queryHits],
               tf = irgr$ID[ovlpsRes$subjectHits])

irgenomewide = ovlps %>% 
  group_by(locus_tag) %>% 
  summarise(count = n())

# adding nrtx info
irs = left_join(x = nrtxsep,
                y = irgenomewide,
                by = "locus_tag") %>% 
  drop_na() %>% 
  group_by(representative) %>% 
  summarise(IRs = count[1])
