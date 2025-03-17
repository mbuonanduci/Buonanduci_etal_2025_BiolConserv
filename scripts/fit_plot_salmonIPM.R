## Fit multi-population IPM to three Willapa Bay index reaches:
## Ellsworth Cr, Canon River, Lower Salmon Cr

## Code originally developed by Eric Buhle
## Code adapted/modified by Michele Buonanduci

library(salmonIPM)       # available at https://github.com/ebuhle/salmonIPM
library(tidyverse)       # data wrangling
library(posterior)       # working with posterior samples
library(here)            # file system paths
library(scales)          # plotting - commas
library(cowplot)         # plotting - multiple panels
library(magicaxis)       # plotting - pretty labels
library(Hmisc)           # PQE - binomial confidence intervals

# Functions to calculate probability of quasi-extinction
source(here("Buonanduci_etal_2025", "scripts","func_PQE.R")) 

# Use multiple cores
options(mc.cores = parallel::detectCores(logical = FALSE))

# Toggle to fit new models (TRUE) or load previous model fits (FALSE)
fit_new <- FALSE


# Load chum data & covariates for Willapa index reaches (1984 - 2022) -----------
all_dat <- here("Buonanduci_etal_2025", "data", "Willapa_chum_input.csv") %>% read_csv()


# Load projected covariates (2023-2080) -------------

# WA coast sea surface temp, spring following spawning
WA_SST_proj <- read_csv(here("Buonanduci_etal_2025", "data","WA_SST_projections.csv"))

WA_SST_proj_canesm2 <- WA_SST_proj %>%
  filter(model == 'canesm2') %>% dplyr::select(-model)
WA_SST_proj_ccsm4 <- WA_SST_proj %>%
  filter(model == 'ccsm4') %>%  dplyr::select(-model)
WA_SST_proj_giss_e2_h <- WA_SST_proj %>%
  filter(model == 'giss_e2_h') %>%  dplyr::select(-model)
WA_SST_proj_noresm1_m <- WA_SST_proj %>%
  filter(model == 'noresm1_m') %>%  dplyr::select(-model)

# WA coast upwelling (CUTI), spring following spawning 
WA_upwell_proj <- read_csv(here("Buonanduci_etal_2025", "data", "WA_upwell_projections.csv"))

# Old-growth structural index (OGSI) across Ellsworth Creek watershed, year of spawning
# Note: These are simple approximations of (a) a future intensive harvest scenario and
#       (b) an ecological forest management (EFM) scenario
OGSI_scenarios <- read_csv(here("Buonanduci_etal_2025", "data", "Ellsworth_OGSI_scenarios.csv"))

OGSI_harvest <- OGSI_scenarios %>%
  filter(scenario == 'harvest') %>% dplyr::select(-scenario)
OGSI_EFM <- OGSI_scenarios %>%
  filter(scenario == 'EFM') %>% dplyr::select(-scenario)


# Create new data frames that include rows for Ellsworth projections ---------
# Center and scale covariates using mean/sd of OBSERVED data

ells_dat_proj <- tibble(pop = 1, A = 1.2, watershed = "Ellsworth", year = c(2023:2080), S_obs = NA, 
                        n_age3_obs = 0, n_age4_obs = 0, n_age5_obs = 0,
                        n_H_obs = 0, n_W_obs = 0, fit_p_HOS = FALSE, B_take_obs = 0, 
                        F_rate = 0.08, # assume harvest fraction of 8% 
                        Spartina_acre = 0.0183, # assume Spartina acreage holds at near-zero
                        hatch_release = 2.5e6, # assume hatchery releases stay around 2.5mil goal
                        forecast = TRUE)

# Function to center and scale covariates using mean/sd of OBSERVED data
ells_rescale <- function(dat){
  dat %>%
    mutate(SST_spr_mean = (SST_spr_mean - mean(all_dat$SST_spr_mean))/sd(all_dat$SST_spr_mean)) %>%
    mutate(upwell_spr = (upwell_spr - mean(all_dat$upwell_spr))/sd(all_dat$upwell_spr)) %>%
    mutate(Spartina_acre = (Spartina_acre - mean(all_dat$Spartina_acre))/sd(all_dat$Spartina_acre)) %>%
    mutate(hatch_release = (hatch_release - mean(all_dat$hatch_release))/sd(all_dat$hatch_release)) %>%
    mutate(OGSI_mean = (OGSI_mean - mean(all_dat$OGSI_mean))/sd(all_dat$OGSI_mean)) %>%
    dplyr::select(-OGSI_mean_sbuff50, -OGSI_mean_sbuff100, -OGSI_mean_sbuff250)
}


## SST from CanESM2 climate model -------

# Harvest scenario
ells_dat_proj_canesm2_harv <- ells_dat_proj %>%
  left_join(WA_SST_proj_canesm2) %>% 
  left_join(WA_upwell_proj) %>% 
  left_join(OGSI_harvest)
ells_dat_proj_canesm2_harv <- mutate(all_dat, forecast = FALSE) %>%
  bind_rows(ells_dat_proj_canesm2_harv)
ells_dat_proj_canesm2_harv <- ells_dat_proj_canesm2_harv %>% ells_rescale()

# EFM scenario
ells_dat_proj_canesm2_EFM <- ells_dat_proj %>%
  left_join(WA_SST_proj_canesm2) %>%
  left_join(WA_upwell_proj) %>%
  left_join(OGSI_EFM)
ells_dat_proj_canesm2_EFM <- mutate(all_dat, forecast = FALSE) %>%
  bind_rows(ells_dat_proj_canesm2_EFM)
ells_dat_proj_canesm2_EFM <- ells_dat_proj_canesm2_EFM %>% ells_rescale()

## SST from CCSM4 climate model -------

# Harvest scenario
ells_dat_proj_ccsm4_harv <- ells_dat_proj %>%
  left_join(WA_SST_proj_ccsm4) %>%
  left_join(WA_upwell_proj) %>%
  left_join(OGSI_harvest)
ells_dat_proj_ccsm4_harv <- mutate(all_dat, forecast = FALSE) %>%
  bind_rows(ells_dat_proj_ccsm4_harv)
ells_dat_proj_ccsm4_harv <- ells_dat_proj_ccsm4_harv %>% ells_rescale()

# EFM scenario
ells_dat_proj_ccsm4_EFM <- ells_dat_proj %>%
  left_join(WA_SST_proj_ccsm4) %>%
  left_join(WA_upwell_proj) %>%
  left_join(OGSI_EFM)
ells_dat_proj_ccsm4_EFM <- mutate(all_dat, forecast = FALSE) %>%
  bind_rows(ells_dat_proj_ccsm4_EFM)
ells_dat_proj_ccsm4_EFM <- ells_dat_proj_ccsm4_EFM %>% ells_rescale()

## SST from GISS-E2-H climate model ---------

# Harvest scenario
ells_dat_proj_giss_e2_h_harv <- ells_dat_proj %>%
  left_join(WA_SST_proj_giss_e2_h) %>%
  left_join(WA_upwell_proj) %>%
  left_join(OGSI_harvest)
ells_dat_proj_giss_e2_h_harv <- mutate(all_dat, forecast = FALSE) %>%
  bind_rows(ells_dat_proj_giss_e2_h_harv)
ells_dat_proj_giss_e2_h_harv <- ells_dat_proj_giss_e2_h_harv %>% ells_rescale()

# EFM scenario
ells_dat_proj_giss_e2_h_EFM <- ells_dat_proj %>%
  left_join(WA_SST_proj_giss_e2_h) %>%
  left_join(WA_upwell_proj) %>%
  left_join(OGSI_EFM)
ells_dat_proj_giss_e2_h_EFM <- mutate(all_dat, forecast = FALSE) %>%
  bind_rows(ells_dat_proj_giss_e2_h_EFM)
ells_dat_proj_giss_e2_h_EFM <- ells_dat_proj_giss_e2_h_EFM %>% ells_rescale()

## SST from NorESM1-M climate model --------

# Harvest scenario
ells_dat_proj_noresm1_m_harv <- ells_dat_proj %>%
  left_join(WA_SST_proj_noresm1_m) %>%
  left_join(WA_upwell_proj) %>%
  left_join(OGSI_harvest)
ells_dat_proj_noresm1_m_harv <- mutate(all_dat, forecast = FALSE) %>%
  bind_rows(ells_dat_proj_noresm1_m_harv)
ells_dat_proj_noresm1_m_harv <- ells_dat_proj_noresm1_m_harv %>% ells_rescale()

# EFM scenario
ells_dat_proj_noresm1_m_EFM <- ells_dat_proj %>%
  left_join(WA_SST_proj_noresm1_m) %>%
  left_join(WA_upwell_proj) %>%
  left_join(OGSI_EFM)
ells_dat_proj_noresm1_m_EFM <- mutate(all_dat, forecast = FALSE) %>%
  bind_rows(ells_dat_proj_noresm1_m_EFM)
ells_dat_proj_noresm1_m_EFM <- ells_dat_proj_noresm1_m_EFM %>% ells_rescale()



# Fit IPMs to observed data --------------

# Approach for covariates!
# Freshwater/estuarine covariates that affect habitat quality:
#   Affect the intrinsic productivity parameter, alpha
# Freshwater/estuarine covariates that affect habitat quantity:
#   Affect the maximum recruitment parameter, Rmax
# Oceanographic and hatchery covariates:
#   Affect recruit abundance (post-density dependence)

if(fit_new){
  
  fit_wshed <- salmonIPM(life_cycle = "SS", pool_pops = TRUE, SR_fun = "Ricker", 
                         par_models = list(R ~ SST_spr_mean + upwell_spr + hatch_release,
                                           Rmax ~ Spartina_acre,
                                           alpha ~ OGSI_mean), 
                         prior = list(mu_alpha ~ normal(1, 2)), # default is N(2,5)
                         center = TRUE, scale = TRUE, 
                         fish_data = all_dat, 
                         chains = 4, iter = 2000, warmup = 1000, seed = 123,
                         control = list(adapt_delta = 0.999, max_treedepth = 12)) 
  
  fit_sbuff250 <- salmonIPM(life_cycle = "SS", pool_pops = TRUE, SR_fun = "Ricker", 
                            par_models = list(R ~ SST_spr_mean + upwell_spr + hatch_release,
                                              Rmax ~ Spartina_acre,
                                              alpha ~ OGSI_mean_sbuff250), 
                            prior = list(mu_alpha ~ normal(1, 2)), # default is N(2,5)
                            center = TRUE, scale = TRUE, 
                            fish_data = all_dat, 
                            chains = 4, iter = 2000, warmup = 1000, seed = 123,
                            control = list(adapt_delta = 0.999, max_treedepth = 12)) 
  
  fit_sbuff100 <- salmonIPM(life_cycle = "SS", pool_pops = TRUE, SR_fun = "Ricker", 
                            par_models = list(R ~ SST_spr_mean + upwell_spr + hatch_release,
                                              Rmax ~ Spartina_acre,
                                              alpha ~ OGSI_mean_sbuff100), 
                            prior = list(mu_alpha ~ normal(1, 2)), # default is N(2,5)
                            center = TRUE, scale = TRUE, 
                            fish_data = all_dat, 
                            chains = 4, iter = 2000, warmup = 1000, seed = 123,
                            control = list(adapt_delta = 0.999, max_treedepth = 12)) 
  
  fit_sbuff50 <- salmonIPM(life_cycle = "SS", pool_pops = TRUE, SR_fun = "Ricker", 
                           par_models = list(R ~ SST_spr_mean + upwell_spr + hatch_release,
                                             Rmax ~ Spartina_acre,
                                             alpha ~ OGSI_mean_sbuff50), 
                           prior = list(mu_alpha ~ normal(1, 2)), # default is N(2,5)
                           center = TRUE, scale = TRUE, 
                           fish_data = all_dat, 
                           chains = 4, iter = 2000, warmup = 1000, seed = 123,
                           control = list(adapt_delta = 0.999, max_treedepth = 12)) 
}

## Save/load model fits -------

# Save new model fits
if(fit_new){
  saveRDS(fit_wshed,    file = here("Buonanduci_etal_2025", "R_objects", "fit_wshed.rds"))
  saveRDS(fit_sbuff250, file = here("Buonanduci_etal_2025", "R_objects", "fit_sbuff250.rds"))
  saveRDS(fit_sbuff100, file = here("Buonanduci_etal_2025", "R_objects", "fit_sbuff100.rds"))
  saveRDS(fit_sbuff50,  file = here("Buonanduci_etal_2025", "R_objects", "fit_sbuff50.rds"))
}

# Load previously fit model objects
if(!fit_new){
  fit_wshed <- here("Buonanduci_etal_2025", "R_objects", "fit_wshed.rds") %>% readRDS()
  fit_sbuff250 <- here("Buonanduci_etal_2025", "R_objects", "fit_sbuff250.rds") %>% readRDS()
  fit_sbuff100 <- here("Buonanduci_etal_2025", "R_objects", "fit_sbuff100.rds") %>% readRDS()
  fit_sbuff50 <- here("Buonanduci_etal_2025", "R_objects", "fit_sbuff50.rds") %>% readRDS()
}

# Check posterior summaries and Rhat values
print(fit_wshed, pars = stan_pars("IPM_SS_pp", "hyper"), prob = c(0.025, 0.5, 0.975))
print(fit_sbuff250, pars = stan_pars("IPM_SS_pp", "hyper"), prob = c(0.025, 0.5, 0.975))
print(fit_sbuff100, pars = stan_pars("IPM_SS_pp", "hyper"), prob = c(0.025, 0.5, 0.975))
print(fit_sbuff50, pars = stan_pars("IPM_SS_pp", "hyper"), prob = c(0.025, 0.5, 0.975))

# Get number of divergent transitions
get_divergent_iterations(fit_wshed) %>% sum()
get_divergent_iterations(fit_sbuff250) %>% sum()
get_divergent_iterations(fit_sbuff100) %>% sum()
get_divergent_iterations(fit_sbuff50) %>% sum()



# Fit IPMs to observed + forecast data --------------

if(fit_new){
  
  ## SST from CanESM2 climate model projections ------
  
  # Harvest scenario
  fit_proj_canesm2_harv <- salmonIPM(life_cycle = "SS", pool_pops = TRUE, SR_fun = "Ricker", 
                                     par_models = list(R ~ SST_spr_mean + upwell_spr + hatch_release,
                                                       Rmax ~ Spartina_acre,
                                                       alpha ~ OGSI_mean),    
                                     prior = list(mu_alpha ~ normal(1, 2)), # default is N(2,5)
                                     center = FALSE, scale = FALSE,  # covariates manually scaled based on observed data
                                     fish_data = ells_dat_proj_canesm2_harv, 
                                     chains = 4, iter = 2000, warmup = 1000, seed = 123,
                                     control = list(adapt_delta = 0.999, max_treedepth = 12)) 
  
  # EFM scenario
  fit_proj_canesm2_EFM <- salmonIPM(life_cycle = "SS", pool_pops = TRUE, SR_fun = "Ricker", 
                                    par_models = list(R ~ SST_spr_mean + upwell_spr + hatch_release,
                                                      Rmax ~ Spartina_acre,
                                                      alpha ~ OGSI_mean),
                                    prior = list(mu_alpha ~ normal(1, 2)), # default is N(2,5)
                                    center = FALSE, scale = FALSE,  # covariates manually scaled based on observed data
                                    fish_data = ells_dat_proj_canesm2_EFM, 
                                    chains = 4, iter = 2000, warmup = 1000, seed = 123,
                                    control = list(adapt_delta = 0.999, max_treedepth = 12)) 
  
  ## SST from noresm1_m climate model projections ----
  
  # Harvest scenario
  fit_proj_noresm1_m_harv <- salmonIPM(life_cycle = "SS", pool_pops = TRUE, SR_fun = "Ricker", 
                                       par_models = list(R ~ SST_spr_mean + upwell_spr + hatch_release,
                                                         Rmax ~ Spartina_acre,
                                                         alpha ~ OGSI_mean),
                                       prior = list(mu_alpha ~ normal(1, 2)), # default is N(2,5)
                                       center = FALSE, scale = FALSE,  # covariates manually scaled based on observed data
                                       fish_data = ells_dat_proj_noresm1_m_harv, 
                                       chains = 4, iter = 2000, warmup = 1000, seed = 123,
                                       control = list(adapt_delta = 0.999, max_treedepth = 12)) 
  
  # EFM scenario
  fit_proj_noresm1_m_EFM <- salmonIPM(life_cycle = "SS", pool_pops = TRUE, SR_fun = "Ricker", 
                                      par_models = list(R ~ SST_spr_mean + upwell_spr + hatch_release,
                                                        Rmax ~ Spartina_acre,
                                                        alpha ~ OGSI_mean),
                                      prior = list(mu_alpha ~ normal(1, 2)), # default is N(2,5)
                                      center = FALSE, scale = FALSE,  # covariates manually scaled based on observed data
                                      fish_data = ells_dat_proj_noresm1_m_EFM, 
                                      chains = 4, iter = 2000, warmup = 1000, seed = 1234,
                                      control = list(adapt_delta = 0.999, max_treedepth = 12)) 
  
  ## SST from CCSM4 climate model projections ----
  
  # Harvest scenario
  fit_proj_ccsm4_harv <- salmonIPM(life_cycle = "SS", pool_pops = TRUE, SR_fun = "Ricker", 
                                   par_models = list(R ~ SST_spr_mean + upwell_spr + hatch_release,
                                                     Rmax ~ Spartina_acre,
                                                     alpha ~ OGSI_mean),
                                   prior = list(mu_alpha ~ normal(1, 2)), # default is N(2,5)
                                   center = FALSE, scale = FALSE,  # covariates manually scaled based on observed data
                                   fish_data = ells_dat_proj_ccsm4_harv, 
                                   chains = 4, iter = 2000, warmup = 1000, seed = 1234,
                                   control = list(adapt_delta = 0.999, max_treedepth = 12)) 
  
  # EFM scenario
  fit_proj_ccsm4_EFM <- salmonIPM(life_cycle = "SS", pool_pops = TRUE, SR_fun = "Ricker",
                                  par_models = list(R ~ SST_spr_mean + upwell_spr,
                                                    Rmax ~ Spartina_acre,
                                                    alpha ~ OGSI_mean),
                                  prior = list(mu_alpha ~ normal(1, 2)), # default is N(2,5)
                                  center = FALSE, scale = FALSE,  # covariates manually scaled based on observed data
                                  fish_data = ells_dat_proj_ccsm4_EFM,
                                  chains = 4, iter = 2000, warmup = 1000, seed = 123,
                                  control = list(adapt_delta = 0.999, max_treedepth = 12))
  
  ## SST from GISS-E2-H climate model projections ----
  
  # Harvest scenario
  fit_proj_giss_e2_h_harv <- salmonIPM(life_cycle = "SS", pool_pops = TRUE, SR_fun = "Ricker", 
                                       par_models = list(R ~ SST_spr_mean + upwell_spr + hatch_release,
                                                         Rmax ~ Spartina_acre,
                                                         alpha ~ OGSI_mean),
                                       prior = list(mu_alpha ~ normal(1, 2)), # default is N(2,5)
                                       center = FALSE, scale = FALSE,  # covariates manually scaled based on observed data
                                       fish_data = ells_dat_proj_giss_e2_h_harv, 
                                       chains = 4, iter = 2000, warmup = 1000, seed = 123,
                                       control = list(adapt_delta = 0.999, max_treedepth = 12)) 
  
  # EFM scenario
  fit_proj_giss_e2_h_EFM <- salmonIPM(life_cycle = "SS", pool_pops = TRUE, SR_fun = "Ricker", 
                                      par_models = list(R ~ SST_spr_mean + upwell_spr + hatch_release,
                                                        Rmax ~ Spartina_acre,
                                                        alpha ~ OGSI_mean),
                                      prior = list(mu_alpha ~ normal(1, 2)), # default is N(2,5)
                                      center = FALSE, scale = FALSE,  # covariates manually scaled based on observed data
                                      fish_data = ells_dat_proj_giss_e2_h_EFM, 
                                      chains = 4, iter = 2000, warmup = 1000, seed = 123,
                                      control = list(adapt_delta = 0.999, max_treedepth = 12)) 
  
}


## Save/load model fits -------

# Save new model fits
if(fit_new){
  saveRDS(fit_proj_canesm2_harv, file = here("Buonanduci_etal_2025", "R_objects", "fit_proj_canesm2_harv.rds"))
  saveRDS(fit_proj_canesm2_EFM, file = here("Buonanduci_etal_2025", "R_objects", "fit_proj_canesm2_EFM.rds"))
  saveRDS(fit_proj_ccsm4_harv, file = here("Buonanduci_etal_2025", "R_objects", "fit_proj_ccsm4_harv.rds"))
  saveRDS(fit_proj_ccsm4_EFM, file = here("Buonanduci_etal_2025", "R_objects", "fit_proj_ccsm4_EFM.rds"))
  saveRDS(fit_proj_noresm1_m_harv, file = here("Buonanduci_etal_2025", "R_objects", "fit_proj_noresm1_m_harv.rds"))
  saveRDS(fit_proj_noresm1_m_EFM, file = here("Buonanduci_etal_2025", "R_objects", "fit_proj_noresm1_m_EFM.rds"))
  saveRDS(fit_proj_giss_e2_h_harv, file = here("Buonanduci_etal_2025", "R_objects", "fit_proj_giss_e2_h_harv.rds"))
  saveRDS(fit_proj_giss_e2_h_EFM, file = here("Buonanduci_etal_2025", "R_objects", "fit_proj_giss_e2_h_EFM.rds"))
}

# Load previously fit model objects
if(!fit_new){
  fit_proj_canesm2_harv <- here("Buonanduci_etal_2025", "R_objects", "fit_proj_canesm2_harv.rds") %>% readRDS()
  fit_proj_canesm2_EFM <- here("Buonanduci_etal_2025", "R_objects", "fit_proj_canesm2_EFM.rds") %>% readRDS()
  fit_proj_ccsm4_harv <- here("Buonanduci_etal_2025", "R_objects", "fit_proj_ccsm4_harv.rds") %>% readRDS()
  fit_proj_ccsm4_EFM <- here("Buonanduci_etal_2025", "R_objects", "fit_proj_ccsm4_EFM.rds") %>% readRDS()
  fit_proj_noresm1_m_harv <- here("Buonanduci_etal_2025", "R_objects", "fit_proj_noresm1_m_harv.rds") %>% readRDS()
  fit_proj_noresm1_m_EFM <- here("Buonanduci_etal_2025", "R_objects", "fit_proj_noresm1_m_EFM.rds") %>% readRDS()
  fit_proj_giss_e2_h_harv <- here("Buonanduci_etal_2025", "R_objects", "fit_proj_giss_e2_h_harv.rds") %>% readRDS()
  fit_proj_giss_e2_h_EFM <- here("Buonanduci_etal_2025", "R_objects", "fit_proj_giss_e2_h_EFM.rds") %>% readRDS()
}

# Get number of divergent transitions
get_divergent_iterations(fit_proj_canesm2_harv) %>% sum()
get_divergent_iterations(fit_proj_canesm2_EFM) %>% sum()
get_divergent_iterations(fit_proj_ccsm4_harv) %>% sum()
get_divergent_iterations(fit_proj_ccsm4_EFM) %>% sum()
get_divergent_iterations(fit_proj_noresm1_m_harv) %>% sum()
get_divergent_iterations(fit_proj_noresm1_m_EFM) %>% sum()
get_divergent_iterations(fit_proj_giss_e2_h_harv) %>% sum()
get_divergent_iterations(fit_proj_giss_e2_h_EFM) %>% sum()


# Plotting  --------------------

## Extract posteriors --------

# Extract posteriors: watershed-scale model
post_wshed <- as_draws_rvars(fit_wshed) %>% 
  mutate_variables(`log(alpha[1])` = log(alpha[1]), `log(alpha[2])` = log(alpha[2]), `log(alpha[3])` = log(alpha[3]), 
                   `log(Rmax[1])` = log(Rmax[1]), `log(Rmax[2])` = log(Rmax[2]), `log(Rmax[3])` = log(Rmax[3]),
                   beta_R = as.vector(beta_R), beta_Rmax = as.vector(beta_Rmax),
                   beta_alpha = as.vector(beta_alpha)) %>% 
  as_draws_df() %>% dplyr::select(-`.chain`, -`.iteration`, -`.draw`) %>%
  mutate(scale = "watershed")

# Extract posteriors: 50-m stream buffer model
post_sbuff50 <- as_draws_rvars(fit_sbuff50) %>% 
  mutate_variables(`log(alpha[1])` = log(alpha[1]), `log(alpha[2])` = log(alpha[2]), `log(alpha[3])` = log(alpha[3]), 
                   `log(Rmax[1])` = log(Rmax[1]), `log(Rmax[2])` = log(Rmax[2]), `log(Rmax[3])` = log(Rmax[3]),
                   beta_R = as.vector(beta_R), beta_Rmax = as.vector(beta_Rmax),
                   beta_alpha = as.vector(beta_alpha)) %>% 
  as_draws_df() %>% dplyr::select(-`.chain`, -`.iteration`, -`.draw`) %>%
  mutate(scale = "50-m stream buffer")

# Extract posteriors: 100-m stream buffer model
post_sbuff100 <- as_draws_rvars(fit_sbuff100) %>% 
  mutate_variables(`log(alpha[1])` = log(alpha[1]), `log(alpha[2])` = log(alpha[2]), `log(alpha[3])` = log(alpha[3]), 
                   `log(Rmax[1])` = log(Rmax[1]), `log(Rmax[2])` = log(Rmax[2]), `log(Rmax[3])` = log(Rmax[3]),
                   beta_R = as.vector(beta_R), beta_Rmax = as.vector(beta_Rmax),
                   beta_alpha = as.vector(beta_alpha)) %>% 
  as_draws_df() %>% dplyr::select(-`.chain`, -`.iteration`, -`.draw`) %>%
  mutate(scale = "100-m stream buffer")

# Extract posteriors: 250-m stream buffer model
post_sbuff250 <- as_draws_rvars(fit_sbuff250) %>% 
  mutate_variables(`log(alpha[1])` = log(alpha[1]), `log(alpha[2])` = log(alpha[2]), `log(alpha[3])` = log(alpha[3]), 
                   `log(Rmax[1])` = log(Rmax[1]), `log(Rmax[2])` = log(Rmax[2]), `log(Rmax[3])` = log(Rmax[3]),
                   beta_R = as.vector(beta_R), beta_Rmax = as.vector(beta_Rmax),
                   beta_alpha = as.vector(beta_alpha)) %>% 
  as_draws_df() %>% dplyr::select(-`.chain`, -`.iteration`, -`.draw`) %>%
  mutate(scale = "250-m stream buffer")


post <- bind_rows(post_wshed,
                  post_sbuff50,
                  post_sbuff100,
                  post_sbuff250) %>%
  mutate(scale = factor(scale, levels = c("50-m stream buffer", "100-m stream buffer",
                                          "250-m stream buffer", "watershed")))



## Plot estimated true spawner abundance ------

indices_ells <- which(all_dat$watershed == "Ellsworth")
indices_can <- which(all_dat$watershed == "Canon")
indices_ls <- which(all_dat$watershed == "LowerSalmon")

draws_all <- as_draws_rvars(as.matrix(fit_wshed, c("S", "tau")))

S_ells <- tibble(location = "Ellsworth",
                 year = all_dat$year[indices_ells],
                 S_obs = all_dat$S_obs[indices_ells]) %>%
  mutate(S = draws_all$S[indices_ells],
         tau = draws_all$tau, 
         S_ppd = rvar_rng(rlnorm, n(), log(S), tau))

S_can <- tibble(location = "Canon",
                year = all_dat$year[indices_can],
                S_obs = all_dat$S_obs[indices_can]) %>%
  mutate(S = draws_all$S[indices_can],
         tau = draws_all$tau, 
         S_ppd = rvar_rng(rlnorm, n(), log(S), tau))

S_ls <- tibble(location = "Lower Salmon",
               year = all_dat$year[indices_ls],
               S_obs = all_dat$S_obs[indices_ls]) %>%
  mutate(S = draws_all$S[indices_ls],
         tau = draws_all$tau, 
         S_ppd = rvar_rng(rlnorm, n(), log(S), tau))

S_all <- bind_rows(S_ells, S_can, S_ls)


plt_S <- ggplot(data = S_all, aes(x = year)) +
  geom_line(aes(y = median(S), color = location), lwd = 1) +
  geom_ribbon(aes(ymin = t(quantile(S, 0.025)), ymax = t(quantile(S, 0.975)), fill = location), alpha = 0.5) +
  geom_ribbon(aes(ymin = t(quantile(S_ppd, 0.025)), ymax = t(quantile(S_ppd, 0.975)), fill = location), alpha = 0.3) +
  geom_point(aes(y = S_obs, color = location)) +
  facet_wrap(~location, nrow = 3)+
  labs(x = "Year", y = "Spawners") + 
  scale_color_manual(values = c("cadetblue3", "black", "tan")) +
  scale_fill_manual(values = c("cadetblue3", "black", "tan")) +
  scale_x_continuous(minor_breaks = unique(all_dat$year), expand = expansion(0.01)) +
  scale_y_continuous(labels = comma) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank())

# percent of observations falling within 95% credible interval for observation model?
S_all$S_ppd_hi <- t(quantile(S_all$S_ppd, 0.975))
S_all$S_ppd_lo <- t(quantile(S_all$S_ppd, 0.025))
S_all <- S_all %>%
  mutate(S_obs_ppd = ifelse(S_obs <= S_ppd_hi & S_obs >= S_ppd_lo, 1, 0))
sum(S_all$S_obs_ppd) / length(S_all$S_obs_ppd)


## Combined spawner, alpha, Rmax plots --------

plt_alpha <- ggplot(data = post_wshed) +
  geom_density(mapping = aes(x = `mu_alpha`, y = after_stat(density), color = "Basin-wide"), lwd = 2, alpha = 0.05) +
  geom_density(mapping = aes(x = `log(alpha[1])`, y = after_stat(density), color = "Ellsworth")) +
  geom_density(mapping = aes(x = `log(alpha[2])`, y = after_stat(density), color = "Canon")) +
  geom_density(mapping = aes(x = `log(alpha[3])`, y = after_stat(density), color = "Lower Salmon")) +
  scale_color_manual(name = "Population",
                     values = c("Basin-wide" = "#7080904D", "Canon" = "cadetblue3", "Ellsworth" = "black", "Lower Salmon" = "tan")) +
  labs(y = "Posterior density", x = expression(log(alpha))) +
  xlim(c(-0.3, 2)) +
  theme_bw() +
  theme(plot.margin = margin(15, 5.5, 5.5, 5.5, "pt"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = "none", 
        panel.grid = element_blank())


plt_Rmax <- ggplot(data = post_wshed) +
  geom_density(mapping = aes(x = `mu_Rmax`, y = after_stat(density), color = "Basin-wide"), lwd = 2, alpha = 0.05) +
  geom_density(mapping = aes(x = `log(Rmax[1])`, y = after_stat(density), color = "Ellsworth")) +
  geom_density(mapping = aes(x = `log(Rmax[2])`, y = after_stat(density), color = "Canon")) +
  geom_density(mapping = aes(x = `log(Rmax[3])`, y = after_stat(density), color = "Lower Salmon")) +
  scale_color_manual(name = NULL,
                     values = c("Basin-wide" = "#7080904D", "Canon" = "cadetblue3", "Ellsworth" = "black", "Lower Salmon" = "tan")) +
  labs(y = "Posterior density", x = expression(log(R[max]))) +
  xlim(c(7,14)) +
  theme_bw() +
  guides(color = guide_legend(ncol = 2)) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position = "bottom", 
        legend.direction = "vertical",
        legend.key.size = unit(10, "pt"),
        panel.grid = element_blank())

plot_grid(plt_S, plot_grid(plt_alpha, plt_Rmax, align = "v",
                           nrow = 2, rel_heights = c(0.8, 1)) , 
          nrow = 1, rel_widths =  c(1, 0.6), 
          labels = c("(a)", "(b)"), label_size = 10)

ggsave(here("Buonanduci_etal_2025", "figures", "fig_spawners.png"),
       width = 6.5, height = 4, units = "in", bg = "white", dpi = 350)



## Combined covariate & covariate effects plots -----

# Covariates
p_a <- ggplot(data = all_dat, aes(year, OGSI_mean, color = watershed)) +
  geom_point() + geom_line() +
  scale_color_manual(name = NULL, 
                     labels = c("Canon", "Ellsworth", "Lower Salmon"),
                     values = c("cadetblue3", "black", "tan")) +
  labs(y = "Forest structural\ncomplexity index", x = "Brood year") +
  ylim(c(10, 26)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = c(0.2, 0.2),
        legend.text = element_text(size = 8),
        legend.background = element_blank(),
        legend.key.size = unit(9, "pt"),
        panel.grid = element_blank())

p_b <- ggplot(data = all_dat, aes(year, Spartina_acre)) +
  geom_point(color = "slategray") + geom_line(color = "slategray") +
  scale_y_continuous(labels = comma) +
  labs(y = expression(paste(italic("Spartina"), " acreage")), x = "Brood year") +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid = element_blank())

p_c <- ggplot(data = all_dat, aes(year, hatch_release)) +
  geom_point(color = "slategray") + geom_line(color = "slategray") +
  scale_y_continuous(breaks = c(0, 2.5e6, 5e6, 7.5e6),
                     labels = c("0", "2.5", "5.0", "7.5")) +
  labs(y = "Hatchery releases\n(millions)", x = "Brood year") +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid = element_blank())

p_d <- ggplot(data = all_dat, aes(year, SST_spr_mean)) +
  geom_point(color = "slategray") + geom_line(color = "slategray") +
  labs(y = "Sea surface\ntemperature (\u00b0C)", x = "Brood year") +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid = element_blank())

p_e <- ggplot(data = all_dat, aes(year, upwell_spr)) +
  geom_point(color = "slategray") + geom_line(color = "slategray") +
  labs(y = "Coastal upwelling\n(m/day)", x = "Brood year") +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid = element_blank())

# Posterior probabilities of effects being > or < 0
prob_beta_alpha <- round(sum(post_wshed$beta_alpha > 0) / nrow(post_wshed), 2) # forest structure
prob_beta_Rmax <- round(sum(post_wshed$beta_Rmax < 0) / nrow(post_wshed), 3) # Spartina
prob_beta_R3 <- round(sum(post_wshed$`beta_R[3]` > 0) / nrow(post_wshed), 2) # hatchery releases
prob_beta_R1 <- round(sum(post_wshed$`beta_R[1]` < 0) / nrow(post_wshed), 2) # SST
prob_beta_R2 <- round(sum(post_wshed$`beta_R[2]` > 0) / nrow(post_wshed), 2) # upwelling

p_f <- ggplot(data = post_wshed) +
  geom_density(mapping = aes(x = `beta_alpha`, y = after_stat(density)), color = "slategray", lwd = 1) +
  labs(y = "Posterior density", x = expression(paste("Effect size (", alpha, ")"))) +
  geom_vline(xintercept = 0, color = "slategray", alpha = 0.2) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid = element_blank())
yl <- layer_scales(p_f)$y$range$range; xl <- layer_scales(p_f)$x$range$range
p_f <- p_f +
  annotate("text", size = 2.8,
           x = xl[1] + 0.83*diff(xl), y = yl[1] + 0.8*diff(yl), 
           label = expression(paste("P(", beta, " > 0) = 0.97")))

p_g <- ggplot(data = post_wshed) +
  geom_density(mapping = aes(x = `beta_Rmax`, y = after_stat(density)), color = "slategray", lwd = 1) +
  labs(y = "Posterior density", x = expression(paste("Effect size (", R[max], ")"))) +
  geom_vline(xintercept = 0, color = "slategray", alpha = 0.2) + 
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid = element_blank())
yl <- layer_scales(p_g)$y$range$range; xl <- layer_scales(p_g)$x$range$range
p_g <- p_g +
  annotate("text", size = 2.8,
           x = xl[1] + 0.83*diff(xl), y = yl[1] + 0.8*diff(yl), 
           label = expression(paste("P(", beta, " < 0) = 0.99")))

# range for Rt effects
Rt_range <- c(post_wshed$`beta_R[1]`, post_wshed$`beta_R[2]`, post_wshed$`beta_R[3]`) %>% range()

p_h <- ggplot(data = post_wshed) +
  geom_density(mapping = aes(x = `beta_R[3]`, y = after_stat(density)), color = "slategray", lwd = 1) +
  labs(y = "Posterior density", x = expression(paste("Effect size (", R[t], ")"))) +
  xlim(Rt_range) +
  geom_vline(xintercept = 0, color = "slategray", alpha = 0.2) + 
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid = element_blank())
yl <- layer_scales(p_h)$y$range$range; xl <- layer_scales(p_h)$x$range$range
p_h <- p_h +
  annotate("text", size = 2.8,
           x = xl[1] + 0.83*diff(xl), y = yl[1] + 0.8*diff(yl), 
           label = expression(paste("P(", beta, " > 0) = 0.61")))

p_i <- ggplot(data = post_wshed) +
  geom_density(mapping = aes(x = `beta_R[1]`, y = after_stat(density)), color = "slategray", lwd = 1) +
  labs(y = "Posterior density", x = expression(paste("Effect size (", R[t], ")"))) +
  xlim(Rt_range) +
  geom_vline(xintercept = 0, color = "slategray", alpha = 0.2) + 
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid = element_blank())
yl <- layer_scales(p_i)$y$range$range; xl <- layer_scales(p_i)$x$range$range
p_i <- p_i +
  annotate("text", size = 2.8,
           x = xl[1] + 0.83*diff(xl), y = yl[1] + 0.8*diff(yl), 
           label = expression(paste("P(", beta, " < 0) = 0.79")))

p_j <- ggplot(data = post_wshed) +
  geom_density(mapping = aes(x = `beta_R[2]`, y = after_stat(density)), color = "slategray", lwd = 1) +
  labs(y = "Posterior density", x = expression(paste("Effect size (", R[t], ")"))) +
  xlim(Rt_range) +
  geom_vline(xintercept = 0, color = "slategray", alpha = 0.2) + 
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid = element_blank())
yl <- layer_scales(p_j)$y$range$range; xl <- layer_scales(p_j)$x$range$range
p_j <- p_j +
  annotate("text", size = 2.8,
           x = xl[1] + 0.83*diff(xl), y = yl[1] + 0.8*diff(yl), 
           label = expression(paste("P(", beta, " > 0) = 0.81")))

# Combine

plot_grid( plot_grid(NULL, p_a, p_b, p_c, p_d, p_e, rel_heights = c(0.1, 1, 1, 1, 1, 1),
                     nrow = 6, align = "v", label_size = 10, label_y = 1.08,
                     labels = c("", "(a)", "(b)", "(c)", "(d)", "(e)")),
           plot_grid(NULL, p_f, p_g, p_h, p_i, p_j, rel_heights = c(0.1, 1, 1, 1, 1, 1),
                     nrow = 6, align = "v", label_size = 10, label_y = 1.08,
                     labels = c("", "(f)", "(g)", "(h)", "(i)", "(j)")),
           ncol = 2, rel_widths = c(1, 0.7))


ggsave(here("Buonanduci_etal_2025", "figures", "fig_covariates.png"),
       width = 6.5, height = 7.5, units = "in", bg = "white", dpi = 350)




## Effect of forest structure spatial scale -----

all_OGSI <- all_dat %>% 
  dplyr::select(year, watershed, OGSI_mean, OGSI_mean_sbuff250, OGSI_mean_sbuff100, OGSI_mean_sbuff50) %>% 
  rename("OGSI_mean_wshed" = OGSI_mean) %>%
  pivot_longer(cols = starts_with("OGSI"),
               names_to = "scale",
               names_prefix = "OGSI_mean_",
               values_to = "OGSI_mean") %>%
  mutate(scale = factor(scale, levels = c("sbuff50", "sbuff100", "sbuff250", "wshed")))

# Plot OGSI
p1 <- ggplot(data = all_OGSI, aes(year, OGSI_mean, color = watershed, linetype = scale)) +
  geom_line() +
  scale_color_manual(name = NULL, values = c("cadetblue3", "black", "tan")) +
  scale_linetype_manual(name = NULL, 
                        values = c(1, 2, 3, 6),
                        labels = c("50-m stream buffer", "100-m stream buffer",
                                   "250-m stream buffer", "Watershed")) +
  labs(y = "Forest structural complexity index", x = "Brood year") +
  ylim(c(12,31)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = c(0.4, 0.85),
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.key.size = unit(14, "pt"),
        panel.grid = element_blank())

# Plot effect sizes
p2 <- ggplot(data = post, aes(linetype = scale)) +
  geom_density(mapping = aes(x = `beta_alpha`, y = after_stat(density)), color = "slategray") +
  scale_linetype_manual(name = NULL, 
                        values = c(1, 2, 3, 6),
                        labels = c("50-m stream buffer", "100-m stream buffer",
                                   "250-m stream buffer", "Watershed")) +
  labs(y = "Posterior density", x = expression(paste("Effect size (", alpha, ")"))) +
  geom_vline(xintercept = 0, color = "slategray", alpha = 0.2) + 
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(5.5, 5.5, 5.5, 9), units = "pt"))


# Combine

plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 0.6),
          align = "h", label_size = 10,
          labels = c("(a)", "(b)"))

ggsave(here("Buonanduci_etal_2025", "figures", "supp_fig_OGSI_scale.png"),
       width = 6.5, height = 3.2, units = "in", bg = "white", dpi = 350)


## Plot projected covariates ------------

p_sst <- ggplot() +
  geom_line(data = WA_SST_proj_canesm2, aes(year, SST_spr_mean, color = 'CanESM2')) +
  geom_line(data = WA_SST_proj_ccsm4, aes(year, SST_spr_mean, color = 'CCSM4')) +
  geom_line(data = WA_SST_proj_noresm1_m, aes(year, SST_spr_mean, color = 'NorESM1-M')) +
  geom_line(data = WA_SST_proj_giss_e2_h, aes(year, SST_spr_mean, color = 'GISS-E2-H')) +
  geom_line(data = all_dat[indices_ells,], aes(year, SST_spr_mean, color = 'Observed')) +
  ylab("Sea surface\ntemperature (\u00b0C)") + xlab("Brood year") +
  scale_color_manual(name = NULL,
                     breaks = c('Observed', 'CanESM2', 'CCSM4', 'NorESM1-M', 'GISS-E2-H'),
                     values = c('Observed' = 'black', 'CanESM2' = '#99DDFF', 'CCSM4' = '#DDDDDD', 'NorESM1-M' = '#44BB99', 'GISS-E2-H' = '#EE8866')) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = c(0.2, 0.78),
        legend.text = element_text(size = 8),
        legend.background = element_blank(),
        legend.key.size = unit(8, "pt"))

p_upw <- ggplot() +
  geom_line(data = WA_upwell_proj, aes(year, upwell_spr, color = 'Future')) +
  geom_line(data = all_dat[indices_ells,], aes(year, upwell_spr, color = 'Observed')) +
  xlab("Brood year") + ylab("Coastal upwelling\n(m/day)") +
  scale_color_manual(name = NULL,
                     breaks = c('Observed', 'Future'),
                     values = c('Observed' = 'black', 'Future' = 'slategray3')) +
  ylim(c(-0.055, 0.47)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = c(0.16, 0.89),
        legend.text = element_text(size = 8),
        legend.background = element_blank(),
        legend.key.size = unit(8, "pt"))

p_ogsi <- ggplot() +
  geom_line(data = all_dat[indices_ells,], aes(year, OGSI_mean, color = 'Observed')) +
  geom_line(data = OGSI_harvest, aes(year, OGSI_mean, color = 'Industrial')) +
  geom_line(data = OGSI_EFM, aes(year, OGSI_mean, color = 'Ecological')) +
  ylab("Forest structural\ncomplexity index") + xlab("Brood year") +
  scale_color_manual(name = NULL,
                     breaks = c('Observed', 'Ecological', 'Industrial'),
                     values = c('Observed' = 'black', 'Industrial' = 'gray50', 'Ecological' = '#DADA0E')) +
  ylim(c(17, 26)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = c(0.17, 0.85),
        legend.text = element_text(size = 8),
        legend.background = element_blank(),
        legend.key.size = unit(8, "pt"))

# Combine

plot_grid( p_ogsi, p_sst, p_upw, ncol = 1,
           align = "hv", label_size = 10, 
           labels = c("(a)", "(b)", "(c)"))

ggsave(here("Buonanduci_etal_2025", "figures", "fig_covariates_future.png"),
       width = 3.25, height = 5.5, units = "in", bg = "white", dpi = 350)



## Plot spawner abundance forecasts --------

# Indices for rows corresponding to Ellsworth in forecasting dataframe
indices_ells_proj <- which(ells_dat_proj_canesm2_EFM$watershed == "Ellsworth")

# Function to extract S with uncertainty intervals
func_get_S <- function(fit_proj, dat_proj, indices){
  
  draws <- as_draws_rvars(as.matrix(fit_proj, c("S", "tau")))
  
  S <- tibble(year = dat_proj$year[indices],
              S_obs = dat_proj$S_obs[indices]) %>%
    mutate(S = draws$S[indices],
           tau = draws$tau, 
           S_ppd = rvar_rng(rlnorm, n(), log(S), tau))
  
  return(S)
  
}

S_canesm2_harv <- func_get_S(fit_proj_canesm2_harv, ells_dat_proj_canesm2_harv, indices_ells_proj) %>%
  mutate(climate = "CanESM2", management = "harvest")
S_canesm2_EFM <- func_get_S(fit_proj_canesm2_EFM, ells_dat_proj_canesm2_EFM, indices_ells_proj) %>%
  mutate(climate = "CanESM2", management = "no harvest")

S_ccsm4_harv <- func_get_S(fit_proj_ccsm4_harv, ells_dat_proj_ccsm4_harv, indices_ells_proj) %>%
  mutate(climate = "CCSM4", management = "harvest")
S_ccsm4_EFM <- func_get_S(fit_proj_ccsm4_EFM, ells_dat_proj_ccsm4_EFM, indices_ells_proj) %>%
  mutate(climate = "CCSM4", management = "no harvest")

S_noresm1_m_harv <- func_get_S(fit_proj_noresm1_m_harv, ells_dat_proj_noresm1_m_harv, indices_ells_proj) %>%
  mutate(climate = "NorESM1-M", management = "harvest")
S_noresm1_m_EFM <- func_get_S(fit_proj_noresm1_m_EFM, ells_dat_proj_noresm1_m_EFM, indices_ells_proj) %>%
  mutate(climate = "NorESM1-M", management = "no harvest")

S_giss_e2_h_harv <- func_get_S(fit_proj_giss_e2_h_harv, ells_dat_proj_giss_e2_h_harv, indices_ells_proj) %>%
  mutate(climate = "GISS-E2-H", management = "harvest")
S_giss_e2_h_EFM <- func_get_S(fit_proj_giss_e2_h_EFM, ells_dat_proj_giss_e2_h_EFM, indices_ells_proj) %>%
  mutate(climate = "GISS-E2-H", management = "no harvest")

S_all <- bind_rows(S_canesm2_harv, S_canesm2_EFM,
                   S_ccsm4_harv, S_ccsm4_EFM,
                   S_noresm1_m_harv, S_noresm1_m_EFM,
                   S_giss_e2_h_harv, S_giss_e2_h_EFM)

# Plotting

plt_forecasts <- ggplot(data = filter(S_all, year <= 2080), aes(x = year)) +
  geom_ribbon(aes(ymin = t(quantile(S, 0.05)), ymax = t(quantile(S, 0.95)), fill = management), alpha = 0.3) +
  geom_line(aes(y = median(S), color = management), lwd = 0.7) +
  facet_wrap(~climate) +
  scale_color_manual(name = "Management\nscenario", values = c("gray50", "#AAAA00"), labels = c("Industrial", "Ecological")) +
  scale_fill_manual(name = "Management\nscenario", values = c("gray50", "#AAAA00"), labels = c("Industrial", "Ecological")) +
  scale_x_continuous(name = "Year") +
  scale_y_log10(name = "Spawner abundance",
                breaks = function(l) maglab(na.omit(l), log = TRUE)$tickat,
                labels = label_log()) +
  coord_cartesian(ylim =c(50, 1.5e5), xlim = c(1983, 2081) ) +
  theme_bw() + theme(axis.text = element_text(size = 8),
                     axis.title = element_text(size = 8),
                     strip.text = element_text(size = 8),
                     legend.position = "none", 
                     strip.background = element_blank(),
                     panel.grid = element_blank())


## Calculate and plot PQE ----

# Quasi-extinction threshold (QET) = 50
PQE_canesm2_harv <- PQE_pp_calc(fit_proj_canesm2_harv, ells_dat_proj_canesm2_harv, pop_ID = 1) %>%
  mutate(SST_model = 'CanESM2', model = 'harvest')
PQE_canesm2_EFM <- PQE_pp_calc(fit_proj_canesm2_EFM, ells_dat_proj_canesm2_EFM, pop_ID = 1) %>%
  mutate(SST_model = 'CanESM2', model = 'no harvest')

PQE_ccsm4_harv <- PQE_pp_calc(fit_proj_ccsm4_harv, ells_dat_proj_ccsm4_harv, pop_ID = 1) %>%
  mutate(SST_model = 'CCSM4', model = 'harvest')
PQE_ccsm4_EFM <- PQE_pp_calc(fit_proj_ccsm4_EFM, ells_dat_proj_ccsm4_EFM, pop_ID = 1) %>%
  mutate(SST_model = 'CCSM4', model = 'no harvest')

PQE_giss_e2_h_harv <- PQE_pp_calc(fit_proj_giss_e2_h_harv, ells_dat_proj_giss_e2_h_harv, pop_ID = 1) %>%
  mutate(SST_model = 'GISS-E2-H', model = 'harvest')
PQE_giss_e2_h_EFM <- PQE_pp_calc(fit_proj_giss_e2_h_EFM, ells_dat_proj_giss_e2_h_EFM, pop_ID = 1) %>%
  mutate(SST_model = 'GISS-E2-H', model = 'no harvest')

PQE_noresm1_m_harv <- PQE_pp_calc(fit_proj_noresm1_m_harv, ells_dat_proj_noresm1_m_harv, pop_ID = 1) %>%
  mutate(SST_model = 'NorESM1-M', model = 'harvest')
PQE_noresm1_m_EFM <- PQE_pp_calc(fit_proj_noresm1_m_EFM, ells_dat_proj_noresm1_m_EFM, pop_ID = 1) %>%
  mutate(SST_model = 'NorESM1-M', model = 'no harvest')

PQE_all <-  bind_rows(PQE_canesm2_harv, PQE_canesm2_EFM,
                      PQE_ccsm4_harv, PQE_ccsm4_EFM,
                      PQE_giss_e2_h_harv, PQE_giss_e2_h_EFM,
                      PQE_noresm1_m_harv, PQE_noresm1_m_EFM)

plt_pqe <- PQE_all %>% 
  filter(end_year == "2080") %>%
  ggplot(aes(x = SST_model, y = PQE.PointEst, color = model, fill = model)) +
  geom_col(position = "dodge", color = NA, alpha = 0.5) +
  scale_color_manual(name = "Management\nscenario", values = c("gray50", "#AAAA00"), labels = c("Industrial\nharvest", "Ecological\nforest mgmt.")) + 
  scale_fill_manual(name = "Management\nscenario", values = c("gray50", "#AAAA00"), labels = c("Industrial\nharvest", "Ecological\nforest mgmt.")) + 
  geom_pointrange(aes(ymin = PQE.Lower, ymax = PQE.Upper), size = 0.1, linewidth = 0.7, 
                  position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0.05, color = "gray50", lty = 3) +
  labs(x = "Climate model", y = "Quasi-extinction\nprobability") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position = c(0.6, 0.8), 
        legend.direction = "horizontal")


# Combine
plot_grid( plt_forecasts, plt_pqe, rel_heights = c(1, 0.6),
           align = "hv", label_size = 10, 
           labels = c("(a)", "(b)"),
           ncol = 1)

ggsave(here("Buonanduci_etal_2025", "figures", "fig_forecast_PQE.png"),
       width = 6.5, height = 5, units = "in", bg = "white", dpi = 350)


# Probability of population dropping below 1000 (min observed in Ellsworth = 1160)
P1000_canesm2_harv <- PQE_pp_calc(fit_proj_canesm2_harv, ells_dat_proj_canesm2_harv, pop_ID = 1, QET = 1000) %>%
  mutate(SST_model = 'CanESM2', model = 'harvest')
P1000_canesm2_EFM <- PQE_pp_calc(fit_proj_canesm2_EFM, ells_dat_proj_canesm2_EFM, pop_ID = 1, QET = 1000) %>%
  mutate(SST_model = 'CanESM2', model = 'no harvest')

P1000_ccsm4_harv <- PQE_pp_calc(fit_proj_ccsm4_harv, ells_dat_proj_ccsm4_harv, pop_ID = 1, QET = 1000) %>%
  mutate(SST_model = 'CCSM4', model = 'harvest')
P1000_ccsm4_EFM <- PQE_pp_calc(fit_proj_ccsm4_EFM, ells_dat_proj_ccsm4_EFM, pop_ID = 1, QET = 1000) %>%
  mutate(SST_model = 'CCSM4', model = 'no harvest')

P1000_giss_e2_h_harv <- PQE_pp_calc(fit_proj_giss_e2_h_harv, ells_dat_proj_giss_e2_h_harv, pop_ID = 1, QET = 1000) %>%
  mutate(SST_model = 'GISS-E2-H', model = 'harvest')
P1000_giss_e2_h_EFM <- PQE_pp_calc(fit_proj_giss_e2_h_EFM, ells_dat_proj_giss_e2_h_EFM, pop_ID = 1, QET = 1000) %>%
  mutate(SST_model = 'GISS-E2-H', model = 'no harvest')

P1000_noresm1_m_harv <- PQE_pp_calc(fit_proj_noresm1_m_harv, ells_dat_proj_noresm1_m_harv, pop_ID = 1, QET = 1000) %>%
  mutate(SST_model = 'NorESM1-M', model = 'harvest')
P1000_noresm1_m_EFM <- PQE_pp_calc(fit_proj_noresm1_m_EFM, ells_dat_proj_noresm1_m_EFM, pop_ID = 1, QET = 1000) %>%
  mutate(SST_model = 'NorESM1-M', model = 'no harvest')

P1000_all <-  bind_rows(P1000_canesm2_harv, P1000_canesm2_EFM,
                        P1000_ccsm4_harv, P1000_ccsm4_EFM,
                        P1000_giss_e2_h_harv, P1000_giss_e2_h_EFM,
                        P1000_noresm1_m_harv, P1000_noresm1_m_EFM)

P1000_all %>% 
  filter(end_year == "2080") %>%
  ggplot(aes(x = SST_model, y = PQE.PointEst, color = model, fill = model)) +
  geom_col(position = "dodge", color = NA, alpha = 0.5) +
  scale_color_manual(name = "Management\nscenario", values = c("gray50", "#AAAA00"), labels = c("Industrial\nharvest", "Ecological\nforest mgmt.")) + 
  scale_fill_manual(name = "Management\nscenario", values = c("gray50", "#AAAA00"), labels = c("Industrial\nharvest", "Ecological\nforest mgmt.")) + 
  geom_pointrange(aes(ymin = PQE.Lower, ymax = PQE.Upper), size = 0.1, linewidth = 0.7, 
                  position = position_dodge(width = 0.9)) +
  labs(x = "Climate model", y = "Probability of 4-year running mean spawner\nabundance dropping below 1,000 individuals") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position = c(0.7, 0.85), 
        legend.direction = "horizontal")

ggsave(here("Buonanduci_etal_2025", "figures", "supp_fig_forecast_PQE_1000.png"),
       width = 6.5, height = 3, units = "in", bg = "white", dpi = 350)
