# alorenzetti 202007

# description #####
# this script will compute codon usage
# it should help us to understand the 
# translational landscape of halo

# coRdon will help us to compute codon usage ####
# taking loaded object and keeping only CDS
nrtxseqscds = nrtxseqs[str_detect(string = names(nrtxseqs), 
                                  pattern = "_t|_r|_s",
                                  negate = T)]

CT = codonTable(nrtxseqscds)

# getting codon adaptation index
# using ribosomal proteins as the reference set
# riboProts = dictProd[str_detect(string = dictProd$product.y, pattern = "ribosomal"),] %>% 
#   dplyr::select(subject_id) %>% 
#   unlist() %>% 
#   unname()

# now that we have protein abundance measures
# I am using the 5% most abundant proteins
# across all timepoints; it was generated
# in one of the last scripts in this series
# so for pragmatic reasons it is going to be
# input right here
mostAbundProtsLysate = c("VNG_0161G", "VNG_0166G", "VNG_0234C", "VNG_0258H", "VNG_0259G", 
                         "VNG_0287a", "VNG_0321G", "VNG_0330G", "VNG_0394C", "VNG_0457G", 
                         "VNG_0491G", "VNG_0524G", "VNG_0683C", "VNG_0771G", "VNG_0787G", 
                         "VNG_0790G", "VNG_1089G", "VNG_1104G", "VNG_1105G", "VNG_1128G", 
                         "VNG_1133G", "VNG_1143G", "VNG_1157G", "VNG_1294G", "VNG_1412H", 
                         "VNG_1414G", "VNG_1524C", "VNG_1541G", "VNG_1542G", "VNG_1668G", 
                         "VNG_1689G", "VNG_1690G", "VNG_1691G", "VNG_1692G", "VNG_1697G", 
                         "VNG_1698G", "VNG_1703G", "VNG_1711G", "VNG_1713G", "VNG_1715G", 
                         "VNG_1768G", "VNG_1793C", "VNG_1873G", "VNG_2001G", "VNG_2010G", 
                         "VNG_2093G", "VNG_2096G", "VNG_2122G", "VNG_2138G", "VNG_2139G", 
                         "VNG_2143G", "VNG_2226G", "VNG_2273H", "VNG_2293G", "VNG_2349G", 
                         "VNG_2443G", "VNG_2467G", "VNG_2486G", "VNG_2513G", "VNG_2514G", 
                         "VNG_2574G", "VNG_2600G", "VNG_2604G", "VNG_2616G", "VNG_2648G", 
                         "VNG_2649G", "VNG_2654G", "VNG_2657G", "VNG_2679G", "VNG_6270G", 
                         "VNG_6294G", "VNG_6309G", "VNG_6315G", "VNG_6316G", "VNG_6317G"
)

riboProtsLog = getID(CT) %in% mostAbundProtsLysate

cais = CAI(CT, id_or_name2 = "11",
           subsets = list(ribo = CT[riboProtsLog,]),
           alt.init = T) %>% 
  as.numeric() %>% 
  unlist() %>%
  unname()

# creating CAI object
cai = tibble(locus_tag = getID(CT),
             cai = cais)
