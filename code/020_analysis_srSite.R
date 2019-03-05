
#-----Setup-------

  library("tidyverse")
  library("DataExplorer")
  library("gridExtra")
  library("GGally")
  library("rstan")
  library("rstanarm")
  library("lubridate")
  library("ggridges")
  source("code/010_data_get.R")

  # make stan run on multiple cores
  options(mc.cores = parallel::detectCores())

  # reset ggplot default theme (rstanarm changes it for some reason)
  theme_set(theme_grey())
  
  
#------Data prep------
  
  datSR <- datAnimals %>%
    dplyr::filter(Method1 != 0
                  , Method2 != 0
                  , LocalBlock != 2
                  ) %>%
    dplyr::count(GSiteID,SurveyID,TransectID,LocalBlock,SurveyDate,SiteName,Subregion,tax,SpeciesID) %>%
    dplyr::group_by(GSiteID,SurveyID,TransectID,LocalBlock,SurveyDate,SiteName,Subregion,tax) %>%
    dplyr::summarise(ind = sum(n)
                     , srSite = n() ########## SET Y HERE ##########
                     ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Year = year(SurveyDate)
                  , medYear = median(Year)
                  , time = (Year - medYear)
                  )

  
#-----Set values------
  
  Y <- "srSite" # dependent variable and used in file names for output
  
  dir.create(paste0("out/",Y)) # make sure folder exists for outputs
  
  # variables to explore
  varExp <- c(get("Y")
              , colnames(datSR)
              ) %>%
    unique() %>%
    grep(pattern = "medTime",invert = TRUE,value=TRUE)
  
  datExp <- datSR %>%
    dplyr::select(varExp) %>%
    na.omit() %>%
    mutate_at(vars(1), as.numeric) %>%
    mutate_if(is.difftime,as.numeric)%>%
    write_csv(paste0("out/",Y,"/",Y,"_explore.csv"))
  
  
#-----Variables explore-----
  
  #  Missing values
  
  plot_missing(datExp)
  ggsave(paste0("out/",Y,"/",Y,"_explore_missing.png"))
  
  
  # Discrete variables
  
  plot_bar(datExp)
  ggsave(paste0("out/",Y,"/",Y,"_explore_discreteCount.png"))
  
  
  # Continuous variables
  
  plot_histogram(datExp)
  ggsave(paste0("out/",Y,"/",Y,"_explore_continuousCount.png"))
  
  
  
  # Everything vs. Y
  
  plot_boxplot(datExp, by = Y)
  ggsave(paste0("out/",Y,"/",Y,"_explore_CatvsY.png"))
  
  plot_scatterplot(datExp, by = Y)
  ggsave(paste0("out/",Y,"/",Y,"_explore_ContvsY.png"))
  
  
  # Outliers
  
  ggplot(datExp %>%
           dplyr::select(varExp) %>%
           dplyr::select_if(is.numeric) %>%
           tidyr::gather(variable,value,2:ncol(.)) %>%
           dplyr::mutate(id = row_number())
         , aes(value,id)
         ) +
    geom_point(alpha = 0.5) +
    facet_wrap(~variable, scales = "free")
  
  ggsave(paste0("out/",Y,"/",Y,"_explore_outliersNumeric.png"))
  
  ggplot(datExp %>%
           dplyr::select(varExp) %>%
           dplyr::select_if(is.character) %>%
           tidyr::gather(variable,value,2:ncol(.)) %>%
           dplyr::mutate(id = row_number())
         , aes(value,id)
         ) +
    geom_point(alpha = 0.5) +
    facet_wrap(~variable, scales = "free")
  
  ggsave(paste0("out/",Y,"/",Y,"_explore_outliersCharacter.png"))
  
  
  # Homogeneity of variance
    # check residuals - plot residuals vs. fitted values
  
  # Normally distributed
    # check distribution of residuals - they should be normally distributed
  
  
  # Collinearity among the covariates
  
  collinearity <- ggpairs(datExp %>%
                            dplyr::select(varExp) %>%
                            dplyr::select_if(function(x) length(levels(as.factor(x))) < 14|is.numeric(x))
                          )
  
  png(paste0("out/",Y,"/",Y,"_explore_collinearity.png")
      , width=2000
      , height=2000
      , res=300
      , bg = "transparent"
      )
  
  collinearity
  
  dev.off()
  
  
#----Variables model------
  
  # variables to model
  varMod <- c(get("Y")
              , "time"
              , "SiteName"
              , "tax"
              )
  
  datMod <- datExp %>%
    dplyr::select(varMod)

#---- Model------

  # If the model has not already been run, run it and save the result to the R directory
  # It it exists, load it
  if(file.exists(paste0("out/",Y,"/",Y,"_mod.rds")) == FALSE){
    
    mod <-
      
      stan_glm(as.formula(paste(get("Y"),"~ time*SiteName + time*tax"))
               , data = datMod
               , family = neg_binomial_2
               , chains = 6
               #, adapt_delta = 0.99
               , iter = 2000
               )
    
    saveRDS(mod,paste0("out/",Y,"/",Y,"_mod.rds"))
    
    } else {
      
      mod <- readRDS(paste0("out/",Y,"/",Y,"_mod.rds"))
      
      }


#------Diagnostics------
  
  # Trace plot
  if(file.exists(paste0("out/",Y,"/",Y,"_diagnostic_trace.png")) == FALSE){
    ggsave(paste0("out/",Y,"/",Y,"_diagnostic_trace.png"),stan_trace(mod))
  }

  # Rhat plot
  if(file.exists(paste0("out/",Y,"/",Y,"_diagnostic_Rhat.png")) == FALSE){
    ggsave(paste0("out/",Y,"/",Y,"_diagnostic_Rhat.png"),plot(mod, "rhat_hist"))
  }
  
  # prop0 plot
  if(file.exists(paste0("out/",Y,"/",Y,"_diagnostic_prop0.png")) == FALSE){
    prop_zero <- function(x) mean(x == 0)
    ggsave(paste0("out/",Y,"/",Y,"_diagnostic_prop0.png")
           , pp_check(mod
                      , plotfun = "stat"
                      , stat = "prop_zero"
                      )
           )
    }
  
  # Model fit
  if(file.exists(paste0("out/",Y,"/",Y,"_diagnostic_fit.png")) == FALSE){
    a <- pp_check(mod) + coord_cartesian(xlim = c(0,2)) + labs(title = "0-2")
    b <- pp_check(mod) + coord_cartesian(xlim = c(0,quantile(mod$fitted.values,probs=0.5))) + labs(title = "0-median pred")
    c <- pp_check(mod) + coord_cartesian(xlim = c(0,max(mod$fitted.values))) + labs(title = "0-max pred")
    d <- grid.arrange(a,b,c)
    ggsave(paste0("out/",Y,"/",Y,"_diagnostic_fit.png"),d)
  }
  
  # Scatterplot of two test statistics
  if(file.exists(paste0("out/",Y,"/",Y,"_diagnostic_scatter2d.png")) == FALSE){
    ggsave(paste0("out/",Y,"/",Y,"_diagnostic_scatter2d.png")
           , pp_check(mod, plotfun = "stat_2d", stat = c("mean", "sd"))
    )
  }
  
  
#---- Predict----
  
  # Use the model to predict results over variables of interest
  modPred <- datMod %>%
    dplyr::select(get("varMod")[-1]) %>%
    unique() %>%
    dplyr::mutate(col = row.names(.)) %>%
    dplyr::left_join(as_tibble(posterior_predict(mod
                                                 , newdata = .
                                                 #, re.form = NA
                                                 )
                               ) %>%
                       tibble::rownames_to_column(var = "row") %>%
                       tidyr::gather(col,value,2:ncol(.))
                     ) %>%
    (function(x) dplyr::bind_cols(x %>% dplyr::select(-value),value = as.numeric(x$value))) %>%
    dplyr::mutate(Year = unique(datSR$medYear) + time) %>%
    write_csv(paste0("out/",Y,"/",Y,"_res_pred.csv"))
  
  # summarise the results
  modRes <- as_tibble(modPred) %>%
    dplyr::group_by_(.dots=c(varMod[-1],"Year")) %>%
    dplyr::summarise(n = n()
                     , nCheck = nrow(as_tibble(mod))
                     , modMedian = quantile(value,0.5)
                     , modMean = mean(value)
                     , modci90lo = quantile(value, 0.05)
                     , modci90up = quantile(value, 0.95)
                     ) %>%
    write_csv(paste0("out/",Y,"/",Y,"_res_res.csv"))
  
  tabRes <- modRes %>%
    dplyr::group_by(SiteName,tax) %>%
    dplyr::filter(Year == max(Year)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(text = paste0(round(modMedian,1)," (",round(modci90lo,1)," to ",round(modci90up,1),")")) %>%
    write_csv(paste0("out/",Y,"/",Y,"_res_last.csv"))
  
  
#------Residuals--------
  
  residuals <- tibble(residual = mod$residuals) %>%
    dplyr::bind_cols(datMod) %>%
    dplyr::group_by(time) %>%
    dplyr::mutate(sdResid = sd(residual)
                  , standResid = residual/sdResid
                  ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-sdResid)
  
  # Everything vs. residuals
  
  ggplot(residuals, aes(time,standResid)) +
    geom_jitter(size = 1
                , alpha=0.5
                ) +
    geom_hline(aes(yintercept = 0),linetype = 2)
  
  ggsave(paste0("out/",Y,"/",Y,"_diagnostic_TimevsResiduals.png"))
  
  plot_scatterplot(residuals
                   , by = "standResid"
                   , size = 0.5
                   , alpha = 0.5
                   )
  
  ggsave(paste0("out/",Y,"/",Y,"_diagnostic_CatvsResiduals.png"))
  
 
#------Slope-----
  
  modSlope <- modPred %>%
    dplyr::group_by(SiteName,tax,row) %>%
    dplyr::mutate(maxYear = max(Year)
                  , minYear = min(Year)
                  , diffYear = maxYear-minYear
                  ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(diffYear >= 10) %>%
    dplyr::mutate(yearFilter = if_else(Year == maxYear
                                       , "max"
                                       , if_else(Year == minYear
                                                 , "min"
                                                 , "NA"
                                                 )
                                       )
                  ) %>%
    dplyr::filter(yearFilter != "NA") %>%
    dplyr::select(SiteName,tax,row,yearFilter,value) %>%
    tidyr::spread(yearFilter,value) %>%
    dplyr::mutate(diff = max-min)
  
  ggplot(modSlope
         , aes(diff,SiteName,fill=factor(..quantile..))
         ) +
    stat_density_ridges(jittered_points = TRUE
                        , point_size = 0.05
                        , point_alpha = 0.5
                        , position = "raincloud"
                        , geom = "density_ridges_gradient"
                        , calc_ecdf = TRUE
                        , quantiles = c(0.025,0.975)
                        ) +
    scale_fill_manual(name = "Probability"
                      , values = c("red","grey","blue")
                      , labels = c("<0.025","0.025-0.975",">0.975")
                      ) +
    geom_vline(xintercept = 0, linetype = 2) +
    facet_wrap(~tax)
  
  ggsave(paste0("out/",Y,"/",Y,"_res_slope.png")
         , width = 19
         , height = 28
         , units = "cm"
         )
  
  
#-----Model vs original plot-----
  
  ggplot() +
    geom_point(data = datExp
               , aes(Year,srSite,colour=tax)
               , alpha = 0.5
               ) +
    geom_line(data = modRes
              , aes(Year,modMedian,colour=tax)
              ) +
    geom_ribbon(data = modRes
                , aes(Year,ymin=modci90lo,ymax=modci90up,fill=tax)
                , alpha = 0.3
                ) +
    facet_wrap(~SiteName) +
    theme(strip.text = element_text(size = 5, angle = 0)
          , axis.text.x = element_text(angle=90, hjust=0)
          )
  
  ggsave(paste0("out/",Y,"/",Y,"_res_plot.png")
         , width = 19
         , height = 28
         , units = "cm"
         )
  