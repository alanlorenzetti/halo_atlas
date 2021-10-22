# alorenzetti 202109

# description ####
# this script will build a function
# to find enrichment of features
# within categorical variables
# determined by the user
# the rationale was taken from
# https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/ and
# http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html

# getting started ####
# the function will take
# a dataframe/tibble object (df)
# a vector to subset the df (subvec)
# the name of the col (testcol)
# containing the categorical
# and the level of categorical (lev)
# var to be tested

# declaring function
enrichAnalysis = function(df, subvec, testcol, lev){
  dfsub = df %>%
    filter(locus_tag %in% subvec)
  
  drawn = subvec %>% length()
  bu = dim(df)[1]
  wb = dfsub[,testcol] %>% 
    unlist() %>%
    unname() %>% 
    grepl(pattern = lev) %>%
    sum()
  wu = df[,testcol] %>% 
    unlist() %>%
    unname() %>% 
    grepl(pattern = lev) %>%
    sum()
  bu = df[,testcol] %>% 
    unlist() %>%
    unname() %>% 
    grepl(pattern = lev) %>%
    `!` %>% 
    sum()
  
  pval = phyper(q = wb, m = wu, n = bu, k = drawn, lower.tail = F)
  
  return(pval)
}

# testing enrichment gvp cluster
gvpenrich = list()
testcols = c("lsmSense", "asRNA", "tps")

for(i in testcols){
  gvpenrich[[i]] = enrichAnalysis(hmaFuncat, gvp1a, i, "yes")
}

# testing enrichment mobilome cluster
mobenrich = list()
mobilome = hmaFuncat %>% 
  filter(cog_category == "Mobilome: prophages, transposons") %>% 
  pull(locus_tag)

for(i in testcols){
  mobenrich[[i]] = enrichAnalysis(hmaFuncat, mobilome, i, "yes")
}

# testing green cluster tp2 vs tp1
greenclusterenrich = list()

for(i in testcols){
  greenclusterenrich[[i]] = enrichAnalysis(hmaFuncat, ptgs$TP2_vs_TP1$Q4$locus_tag, i, "yes")
}
