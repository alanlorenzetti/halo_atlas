# 20220325 alorenzetti

# description ####
# this script will prepare
# files for deployment of
# the halo atlas shinyapp

# getting started ####
# copying latest files
# creating dir to store copied objx and files
if(!dir.exists("data")){dir.create("data")}

file.copy(from = "../results/ht_with_names_for_shiny.RData",
          to = "data/",
          overwrite = T)

# copying legend img
file.copy(from = "../plots/heatmap_legends.png",
          to = "data",
          overwrite = T)

# copying static version
file.copy(from = "../plots/abundanceHeatmap_expanded_en.pdf",
          to = "data",
          overwrite = T)

# deploying
library(rsconnect)
# rsconnect::setAccountInfo(name='alorenzetti', token='7080306C3625803FDA453437CEC051E3', secret='')
deployApp()
