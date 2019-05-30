

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
                  , count = map_dbl(site,~(rpois(1,15+siteAdj)))
                  , obs = factor(obs)
                  , site = factor(site)
                  )
  
  ggplot(dat,aes(count,site)) + geom_density_ridges()
  
  mod <- stan_glmer(count ~ site + (1|obs)
                    , data = dat
                    , family = poisson
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
  
  # prob of new value
  newValue <- 17
  
  prob <- datPred %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(n = n()
                     , med = median(value)
                     , up95 = quantile(value,probs=0.975)
                     , lo95 = quantile(value,probs=0.025)
                     , aboveLike = 100*sum(value>newValue)/n
                     , belowLike = 100*sum(value<newValue)/n
                     , equalLike = 100*sum(value==newValue)/n
                     )
  
  ggplot(datPred, aes(x=value,y=site)) +
    geom_density_ridges() +
    #geom_density_ridges(data = dat, aes(x=count,y=site),alpha=0.5) +
    geom_vline(aes(xintercept = newValue),linetype=2)
  
  