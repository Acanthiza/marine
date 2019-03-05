
# Housekeeping

  # clear R workspace
  rm(list=ls())
  
  # load required packages
  library("tidyverse")
  library("readxl")
  library("RODBC")
  
  
# Data
  con <- odbcConnectAccess2007("I:/DATABASES/BRUV/BDBSA_DataOutput/BRUV_outsideBDBSAcompile.accdb")
  
  dat <- tibble(tables = list(sqlTables(con))) %>%
    tidyr::unnest() %>%
    dplyr::filter(TABLE_TYPE == "TABLE") %>%
    dplyr::mutate(data = map(TABLE_NAME, sqlFetch, channel = con)
                  , newFile = paste0("data/copied_",TABLE_NAME,".csv")
                  )
  
  close(con)
  rm(con)
  
  # copying data to data facilitates future reproducability
  purrr::walk2(dat$data, dat$newFile, write_csv)
  
  
  # remake with copied data
  opcode <- read_csv("data/copied_BRUVS_MaxNdata.csv")
  analysisInfo <- read_csv("data/copied_LU_AnalysisInfo.csv")
  luFish <- read_csv("data/copied_LU_FishedSpecies2.csv") %>%
    dplyr::mutate(Groper = ifelse(SPECIES == "Achoerodus gouldii", 1, 0))
  
  opcodes <- analysisInfo %>%
    dplyr::select(RES_NUM,FieldSeason,SZFullname,InOutSZ,OpCode) %>%
    unique()
  
  dat <- opcode %>%
    dplyr::left_join(analysisInfo) %>%
    dplyr::mutate(spp = paste(Genus, Species, sep = " ")) %>%
    dplyr::left_join(luFish, by = c("spp" = "SPECIES")) %>%
    dplyr::group_by(RES_NUM,FieldSeason,SZFullname,InOutSZ,OpCode) %>%
    dplyr::filter(Groper == 1) %>%
    dplyr::summarise(sum = sum(MaxN)) %>%
    dplyr::full_join(opcodes) %>%
    dplyr::mutate(sum = ifelse(is.na(sum),0,sum)) %>%
    dplyr::group_by(RES_NUM,FieldSeason,SZFullname,InOutSZ) %>%
    dplyr::summarise(mean = mean(sum)
                     , sd = sd(sum)
                     , n = n()
                     , se = sd/sqrt(n)
                     , ymax = mean + se
                     , ymin = mean - se
                     ) %>%
    write_csv("tbl/groper.csv")
  
  
  plots <- dat %>%
    group_by(RES_NUM,SZFullname) %>%
    tidyr::nest() %>%
    dplyr::mutate(plot = map2(data, SZFullname, function(.x,.y) .x %>%
                                ggplot(aes(InOutSZ,mean, fill = InOutSZ)) +
                                
                                geom_col() +
                                geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
                                facet_grid(~FieldSeason) +
                                ggtitle(.y)
                              )
                  , fileName = paste0("fig/",SZFullname,".png")
                  )
  
  purrr::walk2(plots$plot,plots$fileName, ~ ggsave(filename = .y, plot = .x))
  
  

  