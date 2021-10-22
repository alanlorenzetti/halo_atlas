# alorenze 20200525

# description ####
# this script will plot the timecourse
# trajectory of proteins, mRNAs and RPFs
# it also plots sets of proteins

# loading libs ####
source("scripts/loadingLibs.R")

# merging dataframes and adding suffix
# to identify timepoint
timecourse = full_join(alldfs$TP2, alldfs$TP3, by = "locus_tag", suffix = c("_TP2", "_TP3"))

# copying tp4 df and adding tp4 suffix
tp4df = alldfs$TP4
colnames(tp4df) = paste0(colnames(tp4df), "_TP4")
colnames(tp4df)[1] = "locus_tag"

# merging tp2-tp3 to tp4
timecourse = full_join(timecourse, tp4df, by = "locus_tag")

# removing quadrant cols and sigProtOn
timecourse = timecourse %>%
  dplyr::select(-contains("quad")) %>% 
  dplyr::select(-contains("sig")) %>% 
  dplyr::select(-starts_with("border"))

# adding RO 
# beta and TLR cols
timecourse = timecourse %>% 
  mutate(mRNA_TP32 = mRNA_TP2,
         mRNA_TP43 = mRNA_TP3,
         RPF_TP32 = RPF_TP2,
         RPF_TP43 = RPF_TP3,
         protein_TP32 = protein_TP3,
         protein_TP43 = protein_TP4) %>% 
  mutate(lfcse_mRNA_TP32 = lfcse_mRNA_TP2,
         lfcse_mRNA_TP43 = lfcse_mRNA_TP3,
         lfcse_RPF_TP32 = lfcse_RPF_TP2,
         lfcse_RPF_TP43 = lfcse_RPF_TP3) %>% 
  mutate(RO_TP2 = RPF_TP2 - mRNA_TP2,
         RO_TP3 = RPF_TP3 - mRNA_TP3,
         RO_TP4 = RPF_TP4 - mRNA_TP4,
         RO_TP32 = RPF_TP2 - mRNA_TP2,
         RO_TP43 = RPF_TP3 - mRNA_TP3) %>% 
  mutate(beta_TP2 = protein_TP2 - mRNA_TP2,
         beta_TP3 = protein_TP3 - mRNA_TP3,
         beta_TP4 = protein_TP4 - mRNA_TP4,
         beta_TP32 = protein_TP3 - mRNA_TP2,
         beta_TP43 = protein_TP4 - mRNA_TP3) %>% 
  mutate(TLR_TP2 = beta_TP2 - RO_TP2,
         TLR_TP3 = beta_TP3 - RO_TP3,
         TLR_TP4 = beta_TP4 - RO_TP4,
         TLR_TP32 = beta_TP32 - RO_TP2,
         TLR_TP43 = beta_TP43 - RO_TP3)

# adding square root of square sum of distances
timecourse = timecourse %>%
  rowwise() %>% 
  mutate(costFunc_ProtmRNA_Slide = sum((protein_TP3-mRNA_TP2)^2, (protein_TP4-mRNA_TP3)^2) %>% sqrt(),
         costFunc_RPFmRNA_Slide = sum((RPF_TP2-mRNA_TP2)^2, (RPF_TP3-mRNA_TP3)^2) %>% sqrt(),
         costFunc_ProtRPF_Slide = sum((protein_TP3-RPF_TP2)^2, (protein_TP4-RPF_TP3)^2) %>% sqrt(),
         costFunc_ProtmRNA_RPFmRNA_Slide = (costFunc_ProtmRNA_Slide - costFunc_RPFmRNA_Slide)^2 %>% sqrt()) %>%
  
  mutate(costFunc_ProtmRNA_AllTP = sum((protein_TP2-mRNA_TP2)^2, (protein_TP3-mRNA_TP3)^2, (protein_TP4-mRNA_TP4)^2) %>% sqrt(),
         costFunc_RPFmRNA_AllTP = sum((RPF_TP2-mRNA_TP2)^2, (RPF_TP3-mRNA_TP3)^2, (RPF_TP4-mRNA_TP4)^2) %>% sqrt(),
         costFunc_ProtRPF_AllTP = sum((protein_TP2-RPF_TP2)^2, (protein_TP3-RPF_TP3)^2, (protein_TP4-RPF_TP4)^2) %>% sqrt(),
         costFunc_ProtmRNA_RPFmRNA_AllTP = (costFunc_ProtmRNA_AllTP - costFunc_RPFmRNA_AllTP)^2 %>% sqrt()) %>% 
  ungroup()

# finding the 33% lowest values
PMAllTPthr = timecourse$costFunc_ProtmRNA_AllTP %>% quantile(probs = 0.33) %>% unname()
RMAllTPthr = timecourse$costFunc_RPFmRNA_AllTP %>% quantile(probs = 0.33) %>% unname()

PMSlidethr = timecourse$costFunc_ProtmRNA_Slide %>% quantile(probs = 0.33) %>% unname()
RMSlidethr = timecourse$costFunc_RPFmRNA_Slide %>% quantile(probs = 0.33) %>% unname()

# ignoring PMRM for PM or RM not passing the threshold
timecourse = timecourse %>% 
  mutate(costFunc_ProtmRNA_RPFmRNA_AllTP = case_when(costFunc_ProtmRNA_AllTP > PMAllTPthr |
                                                       costFunc_RPFmRNA_AllTP > RMAllTPthr ~ NA_real_,
                                                     TRUE ~ as.numeric(costFunc_ProtmRNA_RPFmRNA_AllTP)),
         costFunc_ProtmRNA_RPFmRNA_Slide = case_when(costFunc_ProtmRNA_Slide > PMSlidethr |
                                                       costFunc_RPFmRNA_Slide > RMSlidethr ~ NA_real_,
                                                     TRUE ~ as.numeric(costFunc_ProtmRNA_RPFmRNA_Slide)))

# assigning the lowest 66% to categorical variable
PMRMAllTPthr = timecourse$costFunc_ProtmRNA_RPFmRNA_AllTP %>% quantile(probs=.33, na.rm = T)
PMRMSlidethr = timecourse$costFunc_ProtmRNA_RPFmRNA_Slide %>% quantile(probs=.33, na.rm = T)

timecourse = timecourse %>% 
  mutate(pmrmPassAllTP = case_when(is.na(costFunc_ProtmRNA_RPFmRNA_AllTP) ~ "no",
                                   costFunc_ProtmRNA_RPFmRNA_AllTP > PMRMAllTPthr ~ "no",
                                   TRUE ~ "yes"),
         pmrmPassSlide = case_when(is.na(costFunc_ProtmRNA_RPFmRNA_Slide) ~ "no",
                                   costFunc_ProtmRNA_RPFmRNA_Slide > PMRMSlidethr ~ "no",
                                   TRUE ~ "yes"))

# pivoting dataframe
# and adding a TP0 with lfc and lfcse = 0
# it is an artificial time point and for
# aesthetics reasons, it is gonna be called
# TP1
tc = timecourse %>%
  dplyr::select(-starts_with("timePoint"),
                -starts_with("cost"),
#                -starts_with("TLR_"),
#                -starts_with("beta_"),
#                -starts_with("RO_"),
                -starts_with("pmrm")) %>% 
  mutate(protein_TP1 = 0,
         RPF_TP1 = 0,
         mRNA_TP1 = 0,
         lfcse_RPF_TP1 = 0,
         lfcse_mRNA_TP1 = 0,
         RO_TP1 = 0,
         beta_TP1 = 0,
         TLR_TP1 = 0) %>% 
  rename_at(vars(matches("^mRNA|^RPF|^protein|^RO|^beta|^TLR")),
            list(~sub("^","lfc_",.))) %>% 
  pivot_longer(cols = contains("TP"),
               names_to = c("measure", "libType", "timepoint"),
               names_pattern = "^(.*)_(.*)_(.*)$",
               values_to = "lfc") %>% 
  pivot_wider(names_from = measure,
              values_from = lfc)
tc$libType = factor(tc$libType, levels=c("mRNA","RPF","protein","RO","beta","TLR"))
#tc$libType = factor(tc$libType, levels=c("mRNA","RPF","protein", ))

# scatters
# alltp
ggplot(data = timecourse,
       aes(x=costFunc_ProtmRNA_AllTP,y=costFunc_RPFmRNA_AllTP)) +
  geom_point(alpha = 0.1) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("PM") + ylab("RM") +
  xlim(c(0,12.5)) + ylim(c(0,12.5))

# slide
ggplot(data = timecourse,
       aes(x=costFunc_ProtmRNA_Slide,y=costFunc_RPFmRNA_Slide)) +
  geom_point(alpha = 0.1) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("PM") + ylab("RM") +
  xlim(c(0,10)) + ylim(c(0,10))

# histogram
# alltp
ggplot(data = timecourse, aes(x = costFunc_ProtmRNA_RPFmRNA_AllTP)) +
  geom_histogram() +
  xlab("PMRM")

# slide
ggplot(data = timecourse, aes(x = costFunc_ProtmRNA_RPFmRNA_Slide)) +
  geom_histogram() +
  xlab("PMRM")

#lower all
genesOfInt = timecourse %>% arrange(costFunc_ProtmRNA_RPFmRNA_AllTP) %>% .$locus_tag %>% head(12)
#lower slide
genesOfInt = timecourse %>% arrange(costFunc_ProtmRNA_RPFmRNA_Slide) %>% .$locus_tag %>% head(12)

# trajectories alltp
tc %>%
  filter(tc$locus_tag %in% genesOfInt) %>% 
  filter(timepoint != "TP32" & timepoint != "TP43") %>% 
  filter(libType != "RO" & libType != "beta" & libType != "TLR") %>% 
  ggplot(aes(x=timepoint, y=lfc, colour = libType, group = libType)) +
  geom_line(size = 1, alpha = 0.75) +
  geom_linerange(aes(ymin = lfc-lfcse, ymax = lfc+lfcse)) +
  facet_wrap(~locus_tag) +
  xlab("Time Point") + ylab("LFC") +
  scale_color_manual(values = c("mRNA" = "#E15759",
                                "RPF" = "#F28E2B",
                                "protein" = "#4E79A7"),
                     name = "Variable")

# assuming is protein is consequence of previous time point contrast
# trajectories slide
# ggplot(data = tc[tc$locus_tag %in% genesOfInt,] %>% filter(timepoint == "TP32" | timepoint == "TP43"),
#        aes(x=timepoint, y=lfc, colour = libType, group = libType)) +
#   geom_line(size = 1, alpha = 0.75) +
#   geom_linerange(aes(ymin = lfc-lfcse, ymax = lfc+lfcse)) +
#   facet_wrap(~locus_tag) +
#   xlab("Time Point") + ylab("LFC") +
#   scale_color_manual(values = c("mRNA" = "#E15759",
#                                 "RPF" = "#F28E2B",
#                                 "protein" = "#4E79A7"),
#                      name = "Variable")
