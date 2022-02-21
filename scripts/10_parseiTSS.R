# alorenzetti 202112

# description ####
# this script will take data
# from ten Caten et al. 2018 (Table S2) 
# on internal Transcription Start Sites
# and parse it accordingly to generate
# gff files for the genome browser

# loading original file and parsing ####
# SuppTable_02_iTSS.xls
itssoritbl = read_tsv(file = "data/SuppTable_02_iTSS.tsv")

itssgff = tibble(acc = itssoritbl$Replicon,
                 src = "tenCaten_et_al_2018",
                 genebiotype = "misc_feature",
                 start = itssoritbl$Position,
                 end = itssoritbl$Position,
                 score = ".",
                 strand = itssoritbl$Strand,
                 misc8 = ".",
                 att = paste0("ID=", itssoritbl$ID, ";",
                              "Name=", itssoritbl$ID, ";",
                              "Note=pos_uncertainty:", itssoritbl$Positional.Uncertainty, ",",
                              "common_name:", itssoritbl$Common.Name, ",",
                              "pvalue:", itssoritbl$p.Value)) %>% 
  mutate(acc = case_when(acc == "chr" ~ "NC_002607.1",
                         acc == "plasmid_pNRC100" ~ "NC_001869.1",
                         acc == "plasmid_pNRC200" ~ "NC_002608.1"))

# writing positive strand file
write_tsv(x = itssgff %>% filter(strand == "+"),
          col_names = F,
          file = "data/iTSS_gff_fwd.gff3")

# writing negative strand file
write_tsv(x = itssgff %>% filter(strand == "-"),
          col_names = F,
          file = "data/iTSS_gff_rev.gff3")  
