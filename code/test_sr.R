

  library(tidyverse)
  library(rstan)  
  library(rstanarm)
  library(ggridges)
  
  nSites <- 3
  nObservers <- 3
  nReplicates <- 20
  
  options(mc.cores = round(parallel::detectCores()*3/4,0))
  
  dat <- tibble(site = 1:nSites) %>%
    dplyr::mutate(siteAdj = map_dbl(site,~round(rnorm(1,0,3)))) %>%
    dplyr::group_by(site,siteAdj) %>%
    tidyr::expand(obs = 1:nObservers) %>%
    dplyr::group_by(site,siteAdj,obs) %>%
    tidyr::expand(rep = 1:nReplicates) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(row = row_number()
                  , value = map_dbl(site,~(rpois(1,15+siteAdj)))
                  , obs = factor(obs)
                  , site = factor(site)
                  )
  
  datRes <- dat %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(n = n()
                     , nCheck = nrow(as_tibble(mod))
                     , modMedian = quantile(value,0.5)
                     , modMean = mean(value)
                     , modci90lo = quantile(value, 0.05)
                     , modci90up = quantile(value, 0.95)
                     ) %>%
    dplyr::ungroup()
  
  ggplot(dat,aes(value,site)) + geom_density_ridges()
  
  mod <- stan_glmer(value ~ site + (1|obs)
                    , data = dat
                    , family = neg_binomial_2()
                    , adapt_delta = 0.99
                    )
  
  datPred <- dat %>%
    dplyr::select(site) %>%
    unique() %>%
    dplyr::mutate(col = row.names(.)) %>%
    dplyr::left_join(as_tibble(posterior_predict(mod
                                                 , newdata = .
                                                 , re.form = NA
                                                 )
                               ) %>%
                       tibble::rownames_to_column(var = "row") %>%
                       tidyr::gather(col,value,2:ncol(.))
                     ) %>%
    (function(x) dplyr::bind_cols(x %>% dplyr::select(-value),value = as.numeric(x$value)))
  
  
  # summarise the results
  modRes <- datPred %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(n = n()
                     , nCheck = nrow(as_tibble(mod))
                     , modMedian = quantile(value,0.5)
                     , modMean = mean(value)
                     , modci90lo = quantile(value, 0.05)
                     , modci90up = quantile(value, 0.95)
                     ) %>%
    dplyr::ungroup()
  
  # prob of new values
  newValues <- tribble(
    ~obs, ~site, ~newValue,
    "DB", 1, 12,
    "DB", 2, 14,
    "AB", 1, 7,
    "AB", 3, 4,
    "CD", 2, 20
    ) %>%
    dplyr::mutate(site = factor(site, levels = levels(dat$site)))
  
  prob <- datPred %>%
    dplyr::left_join(newValues) %>%
    dplyr::group_by(site,obs,newValue) %>%
    dplyr::summarise(n = n()
                     , med = median(value)
                     , lo95 = quantile(value,probs=0.025)
                     , up95 = quantile(value,probs=0.975)
                     , aboveLike = 100*sum(newValue>=value)/n
                     , belowLike = 100*sum(newValue<=value)/n
                     )
  
  ggplot(datPred, aes(value,fill=site)) +
    geom_density() +
    #geom_density_ridges(data = dat, aes(x=count,y=site),alpha=0.5) +
    geom_vline(data = newValues
               , aes(xintercept = newValue, colour = obs)
               ,size=2
               ) +
    facet_wrap(~site) +
    scale_colour_viridis_d(option = "A") +
    scale_fill_viridis_d()
  
  