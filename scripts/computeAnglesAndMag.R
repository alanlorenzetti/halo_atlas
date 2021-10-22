# alorenzetti 20201001

# description #####
# this script will take the vectors
# on the two dimensional space (mRNA vs Prot)
# and compute direction and magnitude
# instructions to use atan2
# function to compute polar coordinates
# http://blog.patrikstas.com/2015/11/05/what-is-difference-between-atan-and-atan2/
# https://gamedev.stackexchange.com/questions/14602/what-are-atan-and-atan2-used-for-in-games

# loading libs ####
source("scripts/loadingLibs.R")

# creating main function ####
computeAngMag = function(df, tp){
  ti = sub(".*_vs_(.*)", "\\1", tp)
  tip1 = sub("(.*)_vs_.*", "\\1", tp)
  
  dfFil = df %>% 
    dplyr::select(locus_tag,
           starts_with("mean") &
             (ends_with(ti) | ends_with(tip1)) &
             (contains("protein_lysate") | contains("rna_total")))
  
  mag = sqrt(((dfFil[,5] - dfFil[,4])^2) + ((dfFil[,3] - dfFil[,2])^2)) %>% 
    unlist(use.names = F)
  
  ang = atan2(x = (dfFil[,5] - dfFil[,4]) %>% unlist(use.names = F),
              y = (dfFil[,3] - dfFil[,2]) %>% unlist(use.names = F)) %>% 
    RadToDeg()
  
  tib = tibble(mag, ang)
  colnames(tib) = paste0(colnames(tib), "_slide", ti, tip1)
  
  return(tib)
}

# running function for meaningful combinations
nuTibs = list()
normNuTibs = list()
combs = c("TP2_vs_TP1",
          "TP3_vs_TP2",
          "TP4_vs_TP3")

for(i in combs){
  nuTibs[[i]] = computeAngMag(abund, i)
}

for(i in combs){
  normNuTibs[[i]] = computeAngMag(abundNorm, i)
}

# unifying datasets
abund2 = bind_cols(abund, nuTibs)
abundNorm2 = bind_cols(abundNorm, normNuTibs)

# adding functional categorization to
# mag and angle dataframes
abund2Funcat = left_join(abund2,
                         dictFunCat,
                         by = c("locus_tag" = "pfeiLocusTag"))

abundNorm2Funcat = left_join(abundNorm2,
                             dictFunCat,
                             by = c("locus_tag" = "pfeiLocusTag"))

# testing variables
abundNorm2Funcat %>% 
  filter(utrSize != "Unknown") %>% 
  ggplot(aes(y = log10(mag_slideTP1TP2),
             x = ang_slideTP1TP2)) +
  scale_x_continuous(breaks = seq(-180, 180, 45)) +
  geom_point(alpha = 0.25) +
  coord_polar(theta = "x",
              direction = -1,
              start = 90 %>% DegToRad()) +
  ggtitle("TP1 -> TP2") +
  xlab("Direction") +
  ylab("log10(Magnitude)")

# angles and categories
abundNorm2Funcat %>% 
  dplyr::select(!contains("abundance")) %>% 
  pivot_longer(cols = contains("slide"),
               names_pattern = "(.*)_slide(.*)",
               names_to = c("var", "slide"),
               values_to = "value") %>% 
  pivot_wider(values_from = "value",
              names_from = "var") %>% 
  ggplot(aes(x = lsmSense, y = ang)) +
  geom_boxplot() +
  geom_beeswarm() +
  facet_wrap(. ~ slide)

# angles and continuous
abundNorm2Funcat %>% 
  dplyr::select(!contains("abundance")) %>% 
  pivot_longer(cols = contains("slide"),
               names_pattern = "(.*)_slide(.*)",
               names_to = c("var", "slide"),
               values_to = "value") %>% 
  pivot_wider(values_from = "value",
              names_from = "var") %>% 
  ggplot(aes(x = ang, y = HL)) +
  geom_point() +
  facet_wrap(. ~ slide)

# creating a list to store plots
magAngPlots = list()

magAngPlots[["TP1_TP2"]] = abundNorm2Funcat %>% 
  ggplot(aes(y = log10(mag_slideTP1TP2),
             x = ang_slideTP1TP2)) +
  scale_x_continuous(breaks = seq(-180, 180, 45)) +
  geom_point(alpha = 0.25) +
  coord_polar(theta = "x",
              direction = -1,
              start = 90 %>% DegToRad()) +
  ggtitle("TP1 -> TP2") +
  xlab("Direction") +
  ylab("log10(Magnitude)")

magAngPlots[["TP2_TP3"]] = abundNorm2Funcat %>% 
  ggplot(aes(y = log10(mag_slideTP2TP3),
             x = ang_slideTP1TP2)) +
  scale_x_continuous(breaks = seq(-180, 180, 45)) +
  geom_point(alpha = 0.25) +
  coord_polar(theta = "x",
              direction = -1,
              start = 90 %>% DegToRad()) +
  ggtitle("TP2 -> TP3") +
  xlab("Direction") +
  ylab("log10(Magnitude)")

magAngPlots[["TP3_TP4"]] = abundNorm2Funcat %>%
  ggplot(aes(y = log10(mag_slideTP3TP4),
             x = ang_slideTP3TP4)) +
  scale_x_continuous(breaks = seq(-180, 180, 45)) +
  geom_point(alpha = 0.25) +
  coord_polar(theta = "x",
              direction = -1,
              start = 90 %>% DegToRad()) +
  ggtitle("TP3 -> TP4") +
  xlab("Direction") +
  ylab("log10(Magnitude)")

# plotting 
ggarrange(plotlist = magAngPlots,
          nrow = 1,
          common.legend = T,
          legend = "right")

nearZero = abundNorm2 %>% 
  filter(ang_slideTP1TP2 > 85 & ang_slideTP1TP2 < 95) %>% 
  dplyr::select(locus_tag) %>% 
  drop_na() %>% 
  unlist(use.names = F) %>% 
  sort()
