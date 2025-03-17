## Functions to calculate probability of quasi-extinction (PQE) through 2060, 2080

## Code originally developed by Eric Buhle
## Code adapted/modified by Michele Buonanduci

# Quasi-extinction risk is calculated as the probability that the 4-year (approximately one generation) 
# running mean spawner abundance falls below the quasi-extinction threshold (QET = 50) 
# at least once during the forecast period

# PQE calculated from single-population forecast models
PQE_calc <- function(mod, fish_data_fore, QET = 50) 
{
  # Draws for S_t
  draws_S <- as_draws_rvars(as.matrix(mod, "S"))
  
  # Data through 2060
  dat60 <- fish_data_fore %>% 
    mutate(end_year = "2060", S = draws_S$S) %>%
    filter(forecast, year <= 2060)
  
  # Data through 2080
  dat80 <- fish_data_fore %>% 
    mutate(end_year = "2080", S = draws_S$S) %>%
    filter(forecast, year <= 2080)
  
  # Calculate PQE
  PQE <- rbind(dat60, dat80) %>% 
    arrange(end_year, year) %>%
    dplyr::select(end_year, year, S) %>% 
    group_by(end_year, year) %>% 
    dplyr::summarize(S = rvar_sum(S)) %>% 
    mutate(mean4_S = c(rep(as_rvar(NA), 3), 
                        rvar_apply(array(4:n()), .margin = 1, 
                                   .f = function(i) rvar_mean(S[(i-3):i]))),
           .after = S) %>% 
    dplyr::summarize(quasi_extinct = rvar_any(mean4_S < QET, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(PQE = binconf(sum(quasi_extinct), ndraws(quasi_extinct), alpha = 0.1)) %>% 
    do.call(data.frame, .) # unpack col with nested data frame
  
  return(PQE)
}


# PQE calculated from multi-population forecast models
PQE_pp_calc <- function(mod, fish_data_fore, pop_ID, QET = 50) 
{
  # Draws for S_t
  draws_S <- as_draws_rvars(as.matrix(mod, "S"))
  
  # Data through 2060
  dat60 <- fish_data_fore %>% 
    filter(pop == pop_ID) %>%
    mutate(end_year = "2060", 
           S = draws_S$S[which(fish_data_fore$pop == pop_ID)]) %>%
    filter(forecast, year <= 2060)
  
  # Data through 2080
  dat80 <- fish_data_fore %>% 
    filter(pop == pop_ID) %>%
    mutate(end_year = "2080", 
           S = draws_S$S[which(fish_data_fore$pop == pop_ID)]) %>%
    filter(forecast, year <= 2080)
  
  # Calculate PQE
  PQE <- rbind(dat60, dat80) %>% 
    arrange(end_year, year) %>%
    dplyr::select(end_year, year, S) %>% 
    group_by(end_year, year) %>% 
    dplyr::summarize(S = rvar_sum(S)) %>% 
    mutate(mean4_S = c(rep(as_rvar(NA), 3), 
                        rvar_apply(array(4:n()), .margin = 1, 
                                   .f = function(i) rvar_mean(S[(i-3):i]))),
           .after = S) %>% 
    dplyr::summarize(quasi_extinct = rvar_any(mean4_S < QET, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(PQE = binconf(sum(quasi_extinct), ndraws(quasi_extinct), alpha = 0.1)) %>% 
    do.call(data.frame, .) # unpack col with nested data frame
  
  return(PQE)
}

