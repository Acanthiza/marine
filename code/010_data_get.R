#------Setup-----

  library("tidyverse")
  library("readxl")
  

#------Import data------

  datAnimals <- read_excel(list.files("data/",pattern = ".xlsx",full.names = TRUE)) %>%
    dplyr::rename(Sanctuary = InsideOutside) %>%
    dplyr::mutate(Subregion = as.character(Subregion)) %>%
    tidyr::separate(TAXONOMIC_GROUP,"phylum",sep=",",remove=FALSE) %>%
    dplyr::mutate(tax = if_else(phylum == "Chordata","Fish","Invert"))
    
  datAlgae <- read_excel(list.files("data/",pattern = ".xlsx",full.names = TRUE),sheet=3) %>%
    dplyr::rename(Sanctuary = InsideOutside) %>%
    dplyr::mutate(Subregion = as.character(Subregion)
                  , points = SumOfAbundance
                  , Macroalgae = if_else(grepl("macro",tolower(Species)),points,0)
                  , SiteName = `Site Name`
                  , Year = SYear
                  )
  