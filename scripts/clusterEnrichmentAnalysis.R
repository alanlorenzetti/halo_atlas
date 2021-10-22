# alorenzetti 202007

# description ####
# this script will take the reg rule clusters
# and dictFunCat objects and compute
# the hypergeometric test for a few
# types of categorical variables
# the rationale was extracted from
# https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/

# loading libs ####
source("scripts/loadingLibs.R")

# hypergeometric p threshold
qthr = 0.01

###################### regRules as primary clusters
# hypergeometric enrichment test of
# the following vars inside
# the regRule clusters
vars = c("lsmSense", "arCOG", "asRNA")

# creating list to store results
enrich = list()
for(i in names(joinedTibble)){
  regGroups = joinedTibble[[i]] %>%
    dplyr::select(regRule) %>% 
    unlist(use.names = F) %>%
    unique()
  for(j in regGroups){
    curRegGroup = joinedTibble[[i]] %>%
      filter(regRule == j)
    for(k in vars){
      curRegGroupVec = curRegGroup %>% 
        dplyr::select(k) %>% 
        unlist(use.names = F)
      
      curRegGroupLvs = curRegGroupVec %>% 
        unique()
      
      for(l in curRegGroupLvs){
        wb = sum(curRegGroupVec == l)
        vecu = dictFunCat %>% 
          dplyr::select(k) %>% 
          unlist(use.names = F)
        wu = sum(vecu == l)
        bu = sum(vecu != l)
        drawn = curRegGroupVec %>% length()
        
        pval = phyper(q= wb, m= wu, n= bu, k = drawn, lower.tail = F)
        
        tib = tibble(regRule = j,
                     criteria = k,
                     level = l,
                     pval = pval)
        enrich[[i]] = bind_rows(enrich[[i]], tib)
      }
    }
  }
}

###################### arCOGs as primary clusters
# hypergeometric enrichment test of
# the following vars inside
# the regRule clusters

# note that here, the universe genes 
# for regRule has to be extracted from
# the regRule object, so total genes would
# be ~600 instead of ~2400
vars = c("lsmSense", "regRule", "asRNA")

for(i in names(joinedTibble)){
  regGroups = joinedTibble[[i]] %>%
    dplyr::select(arCOG) %>% 
    unlist(use.names = F) %>%
    unique()
  for(j in regGroups){
    curRegGroup = joinedTibble[[i]] %>%
      filter(arCOG == j)
    for(k in vars){
      curRegGroupVec = curRegGroup %>% 
        dplyr::select(k) %>% 
        unlist(use.names = F)
      
      curRegGroupLvs = curRegGroupVec %>% 
        unique()
      
      for(l in curRegGroupLvs){
        wb = sum(curRegGroupVec == l)
        
        if(k == "regRule"){
          vecu = joinedTibble[[i]] %>% 
            dplyr::select(k) %>% 
            unlist(use.names = F)
        }else{
          vecu = dictFunCat %>% 
            dplyr::select(k) %>% 
            unlist(use.names = F)
        }
        
        wu = sum(vecu == l)
        bu = sum(vecu != l)
        drawn = curRegGroupVec %>% length()
        
        pval = phyper(q= wb, m= wu, n= bu, k = drawn, lower.tail = F)
        
        tib = tibble(regRule = j,
                     criteria = k,
                     level = l,
                     pval = pval)
        enrich[[i]] = bind_rows(enrich[[i]], tib)
      }
    }
  }
}

# correcting pvalues using BH method
# filtering by pval
for(i in names(enrich)){
  enrich[[i]]$qval = p.adjust(enrich[[i]]$pval)
  enrich[[i]] = enrich[[i]][enrich[[i]]$pval < qthr,]
}
