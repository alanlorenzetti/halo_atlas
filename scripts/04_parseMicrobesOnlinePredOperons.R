# alorenzetti 20200622

# description #####
# this script will take the predicted operons
# by microbes online, and parse them to my own format
# copy and save as csv the table available at
# http://www.microbesonline.org/operons/gnc64091.html

# starting processing ####
# reading csv table 
if(!file.exists("data/microbesOnlineOperons.RData")){
  df = read_delim("http://www.microbesonline.org/operons/gnc64091.named", delim="\t")
  save(df, file = "data/microbesOnlineOperons.RData")
} else {
  load("data/microbesOnlineOperons.RData")
}

# filtering in only needed cols
df = df %>%
  dplyr::select(SysName1,SysName2,bOp)

# function to parse the table
parseOperons = function(df){
  
  m =  dim(df)[1]
  n =  dim(df)[2]

  opCounter = 1
  
  output = tibble()
  
  for(i in 1:m){
    sys1cur = df[i,"SysName1"] %>% unlist() %>% unname()
    sys2cur = df[i,"SysName2"] %>% unlist() %>% unname()
    bopcur = df[i,"bOp"] %>% unlist() %>% unname()
    
    if(i == 1){
      sys1prev = NULL
      sys2prev = NULL
      bopprev = NULL
      
      if(bopcur == TRUE){
        line = tibble(locus_tag = c(sys1cur,sys2cur),
                      operon = rep(paste0("Operon", opCounter)))
      }
      if(bopcur == FALSE){
        line = tibble(locus_tag = c(sys1cur,sys2cur),
                      operon = c(paste0("Operon", opCounter),
                                 paste0("Operon", opCounter + 1)))
        opCounter = opCounter + 2
      }
    }
    
    if(i != 1){
      sys1prev = df[i-1,"SysName1"] %>% unlist() %>% unname()
      sys2prev = df[i-1,"SysName2"] %>% unlist() %>% unname()
      bopprev = df[i-1,"bOp"] %>% unlist() %>% unname()
      m = dim(output)[1]
      opCounterprev = output[m,"operon"] %>% sub("^Operon","",.) %>% as.numeric()
      
      if(sys1cur == sys2prev){
        if(bopcur == TRUE){
          opCounter=opCounterprev
          line = tibble(locus_tag = sys2cur,
                        operon = paste0("Operon", opCounter))
        }
        if(bopcur == FALSE){
          opCounter = opCounterprev + 1
          line = tibble(locus_tag = sys2cur,
                        operon = paste0("Operon", opCounter))
        }
      }
      
      if(sys1cur != sys2prev){
        if(bopcur == TRUE){
          opCounter = opCounterprev + 1
          line = tibble(locus_tag = c(sys1cur,sys2cur),
                        operon = rep(paste0("Operon", opCounter)))
        }
        if(bopcur == FALSE){
          opCounter = opCounterprev + 1
          line = tibble(locus_tag = c(sys1cur,sys2cur),
                        operon = c(paste0("Operon", opCounter),
                                   paste0("Operon", opCounter + 1)))
        }
      }
    }
    output = bind_rows(output,line)
  }
  
  output = output %>%
    group_by(operon) %>%
    add_tally() %>%
    filter(n > 1) %>% 
    mutate(genes = paste0("Operon:",paste0(locus_tag, collapse = ","))) %>% 
    ungroup() %>% 
    dplyr::select(locus_tag = locus_tag,
           operon = genes)
  
  return(output)
}

# exec function
operons = parseOperons(df)
