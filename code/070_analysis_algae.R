
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
  
  
  
#-----Set values------
  
  dir.create("out")
  
  Y <- "propMAlgae" # dependent variable and used in file names for output
  
  if(!file.exists(paste0("out/",Y,"/",Y,"_mod.rds"))) {
    
    dir.create(paste0("out/",Y)) # make sure folder exists for outputs
    
    
#------Data prep------

  dat <- datAlgae %>%
    dplyr::filter() %>%
    dplyr::group_by(GSiteID,SurveyID,TransectID,SurveyDate,SiteName,Composition) %>%
    dplyr::summarise(points = sum(points)
                     , Macroalgae = sum(Macroalgae) 
                     ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Year = year(SurveyDate)
                  , medYear = ((max(Year)-min(Year))/2)+min(Year)
                  , time = (Year - medYear)
                  , UQ(rlang::sym(Y)) :=  Macroalgae/points ########## SET Y HERE ##########
                  )
    

#-----Variables explore-----
    
    # variables to explore
    varExp <- c(get("Y")
                , colnames(dat)
                ) %>%
      unique() %>%
      grep(pattern = "medYear",invert = TRUE, value=TRUE)
    
    datExp <- dat %>%
      dplyr::select(varExp) %>%
      #na.omit() %>%
      mutate_at(vars(1), as.numeric) %>%
      write_csv(paste0("out/",Y,"/",Y,"_explore.csv"))
  
  
#-----Explore data------
  
    #  Missing values
    plot_missing(datExp)
    ggsave(paste0("out/",Y,"/",Y,"_explore_missing.png"))
    
    
    # Count discrete variables
    ggplot(datExp %>%
             dplyr::mutate_if(is.factor,as.character) %>%
             dplyr::select_if(is.character) %>%
             tidyr::gather(variable,value,1:ncol(.))
           ) +
      geom_histogram(aes(value),stat="count") +
      facet_wrap(~variable, scales = "free") +
      theme(axis.text.x=element_text(angle=90, vjust=0.5))
      
    ggsave(paste0("out/",Y,"/",Y,"_explore_CountDiscrete.png"))
    
    
    # Count continuous variables
    ggplot(datExp %>%
             dplyr::select_if(is.numeric) %>%
             tidyr::gather(variable,value,1:ncol(.))
           , aes(value)
           ) +
      geom_histogram() +
      facet_wrap(~variable, scales = "free") +
      theme(axis.text.x=element_text(angle=90, vjust=0.5))
    
    ggsave(paste0("out/",Y,"/",Y,"_explore_CountContinuous.png"))
  
  
    # Y vs. Discrete
    ggplot(datExp %>%
             dplyr::mutate(UQ(rlang::sym(Y)) := factor(!!ensym(Y))) %>%
             dplyr::mutate_if(is.factor,as.character) %>%
             dplyr::select_if(is.character) %>%
             dplyr::mutate(UQ(rlang::sym(Y)) := as.numeric(!!ensym(Y))) %>%
             tidyr::gather(variable,value,2:ncol(.))
           ) +
      geom_boxplot(aes(value,!!ensym(Y))) +
      facet_wrap(~variable, scales = "free") +
      theme(axis.text.x=element_text(angle=90, vjust=0.5))
    
    ggsave(paste0("out/",Y,"/",Y,"_explore_YvsDiscrete.png"))
    
    
    # Y vs. Numeric
    ggplot(datExp %>%
             dplyr::select(varExp) %>%
             dplyr::select_if(is.numeric) %>%
             tidyr::gather(variable,value,2:ncol(.)) %>%
             dplyr::arrange(!!ensym(Y))
           , aes(!!ensym(Y),value)
           ) +
      geom_jitter(alpha = 0.5
                  , width = 0.1
                  , height = 0.1
                  ) +
      facet_wrap(~variable, scales = "free") +
      theme(axis.text.x=element_text(angle=90, vjust=0.5))
    
    ggsave(paste0("out/",Y,"/",Y,"_explore_YvsNumeric.png"))
    
    
    # Homogeneity of variance
    # check residuals - plot residuals vs. fitted values
    
    # Normally distributed
    # check distribution of residuals - they should be normally distributed
    
    
    # Collinearity among the covariates
    p <- ggpairs(datExp %>%
              dplyr::select(varExp) %>%
              dplyr::select_if(function(x) length(levels(as.factor(x))) < 14|is.numeric(x))
            )
    
    ggsave(plot = p
           , paste0("out/",Y,"/",Y,"_explore_collinearity.png")
           , width = 20
           , height = 28
           , units = "cm"
           )
  
  
#----Variables model------
  
    # variables to model
    varMod <- c(get("Y")
                , "time"
                , "SiteName"
                )
    
    datMod <- datExp %>%
      dplyr::select(varMod) %>%
      dplyr::mutate(UQ(rlang::sym(Y)) :=  if_else(!!ensym(Y)==0,0.001,!!ensym(Y))
                    , UQ(rlang::sym(Y)) :=  if_else(!!ensym(Y)>0.999,0.999,!!ensym(Y))
                    )

#-----Model------

  # If the model has not already been run, run it and save the result to the R directory
      mod <-
      
      stan_glm(
                  as.formula(paste(get("Y")," ~ time*SiteName")) # gaussian or beta
                  , data = datMod
                  , family = mgcv::betar
                  )
    
    saveRDS(mod,paste0("out/",Y,"/",Y,"_mod.rds"))
    

#------Diagnostics------
  
    # Model fit
    a <- pp_check(mod) + coord_cartesian(xlim = c(0,quantile(mod$fitted.values,probs=0.05))) + labs(title = "0-quantile(0.1)")
    b <- pp_check(mod) + coord_cartesian(xlim = c(0,quantile(mod$fitted.values,probs=0.5))) + labs(title = "0-quantile(0.5)")
    c <- pp_check(mod) + coord_cartesian(xlim = c(quantile(mod$fitted.values,probs=0.5),max(mod$fitted.values))) + labs(title = "quantile(0.5)-max")
    d <- pp_check(mod) + labs(title = "0-max pred")
    e <- grid.arrange(a,b,c,d)
    ggsave(paste0("out/",Y,"/",Y,"_diagnostic_fit.png")
           , width = 20
           , height = 28
           , units = "cm"
           , e
           )
  
    # Trace plot
    ggsave(paste0("out/",Y,"/",Y,"_diagnostic_trace.png"),stan_trace(mod))
  

    # Rhat plot
    ggsave(paste0("out/",Y,"/",Y,"_diagnostic_Rhat.png"),plot(mod, "rhat_hist"))
  
  
    # prop0 plot
    prop_zero <- function(x) mean(x == 0)
    ggsave(paste0("out/",Y,"/",Y,"_diagnostic_prop0.png")
           , pp_check(mod
                      , plotfun = "stat"
                      , stat = "prop_zero"
                      )
           )
  
    # Scatterplot of two test statistics
    ggsave(paste0("out/",Y,"/",Y,"_diagnostic_scatter2d.png")
           , pp_check(mod, plotfun = "stat_2d", stat = c("mean", "sd"))
           )
  
  
#---- Predict----
  
  # Use the model to predict results over variables of interest
  modPred <- datMod %>%
    dplyr::select(time,SiteName) %>%
    unique() %>%
    dplyr::arrange(time,SiteName) %>%
    dplyr::mutate(col = row.names(.)) %>%
    dplyr::left_join(as_tibble(posterior_predict(mod
                                                 , newdata = .
                                                 , re.form = NA
                                                 )
                               ) %>%
                       tibble::rownames_to_column(var = "row") %>%
                       tidyr::gather(col,value,2:ncol(.))
                     ) %>%
    (function(x) dplyr::bind_cols(x %>% dplyr::select(-value),value = as.numeric(x$value))) %>%
    dplyr::mutate(Year = unique(dat$medYear) + time) %>%
    write_csv(paste0("out/",Y,"/",Y,"_res_pred.csv"))
  
  # summarise the results
   modRes <- as_tibble(modPred) %>%
    dplyr::group_by(time,Year,SiteName) %>%
    dplyr::summarise(n = n()
                     , nCheck = nrow(as_tibble(mod))
                     , modMedian = quantile(value,0.5)
                     , modMean = mean(value)
                     , modci90lo = quantile(value, 0.05)
                     , modci90up = quantile(value, 0.95)
                     , ci = modci90up-modci90lo
                     , text = paste0(round(ci,2)," (",round(modci90lo,2)," to ",round(modci90up,2),")")
                     ) %>%
     dplyr::mutate_if(is.numeric,round,2) %>%
     dplyr::ungroup() %>%
     write_csv(paste0("out/",Y,"/",Y,"_res_res.csv"))
  
  tabRes <- modRes %>%
    dplyr::group_by(time,Year,SiteName) %>%
    dplyr::filter(Year == max(Year)) %>%
    dplyr::ungroup() %>%
    write_csv(paste0("out/",Y,"/",Y,"_res_last.csv"))
  
  
#------Residuals--------
  
  modResid <- tibble(residual = mod$residuals) %>%
    dplyr::bind_cols(datMod) %>%
    dplyr::group_by_at(vars(names(.)[-c(1:2)])) %>%
    dplyr::mutate(sdResid = sd(residual)
                  , standResid = residual/sdResid
                  ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-sdResid) %>%
    write_csv(paste0("out/",Y,"/",Y,"_res_modResid.csv"))
  
  # Y vs. residuals
  ggplot(modResid, aes(!!ensym(Y),residual)) +
    geom_point(size = 2
              #, alpha=0.5
              , position=position_jitter(width=0.3,height=0.2)
              ) +
    geom_hline(aes(yintercept = 0),linetype = 2) +
    facet_wrap(~SiteName) +
    scale_colour_viridis_c()
    #scale_colour_gradient2(mid = "yellow",midpoint=(max(datExp$Transect)-min(datExp$Transect))/2)
  
  ggsave(paste0("out/",Y,"/",Y,"_diagnostic_",Y,"vsResiduals.png"))
  
  
  
  # Time vs. residuals
  ggplot(modResid, aes(time,residual)) +
    geom_point(size = 2
              #, alpha=0.5
              , position=position_jitter(width=0.3,height=0.2)
              ) +
    geom_hline(aes(yintercept = 0),linetype = 2) +
    facet_wrap(~SiteName) +
    scale_colour_viridis_c()
    #scale_colour_gradient2(mid = "yellow",midpoint=(max(datExp$Transect)-min(datExp$Transect))/2)
  
  ggsave(paste0("out/",Y,"/",Y,"_diagnostic_TimevsResiduals.png"))
  
  
  # Everything vs. residuals
  plot_scatterplot(residuals
                   , by = "residual"
                   , size = 0.5
                   , alpha = 0.5
                   )
  
  ggsave(paste0("out/",Y,"/",Y,"_diagnostic_CatvsResiduals.png"))

  

#-----Model vs original plot-----
  
  ggplot() +
    #geom_jitter(alpha = 0.5) +
    geom_point(data = datExp
              , aes(Year
                    ,!!ensym(Y)
                    #, colour = Transect
                    #, label = Transect
                    )
              , alpha = 0.5
              , size = 2
              ) +
    geom_line(data = modRes
              , aes(Year,modMedian)
              ) +
    geom_ribbon(data = modRes
                , aes(Year,ymin=modci90lo,ymax=modci90up)
                , alpha = 0.3
                ) +
    facet_wrap(~SiteName) +
    theme(axis.text.x = element_text(angle=90, hjust=0)) +
    scale_colour_viridis_c() +
    #scale_colour_gradient2(mid = "yellow",midpoint=(max(datExp$Transect)-min(datExp$Transect))/2) +
    treeColFill
  
  ggsave(paste0("out/",Y,"/",Y,"_res_plot.png")
         , width = 19
         , height = 28
         , units = "cm"
         )
  
#------Load for Report------
  
  } else {
    
    datExp <- read_csv(paste0("out/",Y,"/",Y,"_explore.csv"))
    
    mod <- readRDS(paste0("out/",Y,"/",Y,"_mod.rds"))
    
    modResid <- read_csv(paste0("out/",Y,"/",Y,"_res_residuals.csv"))
    
    modPred <- read_csv(paste0("out/",Y,"/",Y,"_res_pred.csv"))
    
    modRes <- read_csv(paste0("out/",Y,"/",Y,"_res_res.csv"))
    
  }