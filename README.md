# Buonanduci_etal_2025_BiolConserv
Data and code accompanying the manuscript 'Forest restoration can bolster salmon population persistence under climate change' by Buonanduci, Buhle, Case, Howe, Robertson, VanBuskirk, and Ettinger, published in Biological Conservation.


[![CC BY 4.0][cc-by-shield]][cc-by]

This information is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by]. Any user of these data ("User" hereafter) is required to cite it appropriately in any publication that results from its use. These data may be actively used by others for ongoing research, so coordination may be necessary to prevent duplicate publication. The User is urged to contact the authors of these data for questions about methodology or results.  The User is encouraged to consider collaboration or co-authorship with authors where appropriate. Misinterpretation of data may occur if used out of context of the original study. Substantial efforts are made to ensure accuracy of the data and documentation, however complete accuracy of data sets cannot be guaranteed. All data are made available as is. Data may be updated periodically and it is the responsibility of the User to check for new versions of the data. The authors and the repository where these data were obtained shall not be liable for damages resulting from any use or misinterpretation of the data.

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg



For reproducibility, the following files are made available:

## Data files

#### Willapa_chum_input.csv
This file contains the data used to fit a hierarchical integrated population model (IPM) to three index populations of chum salmon (*Oncorhynchus keta*) in Willapa Bay. The following columns are included:

- **pop**: unique population identifier
- **watershed**: name of the index reach watershed
- **A**: surveyed length of each index reach (river miles)
- **year**: calendar year corresponding to each spawner abundance observation
- **S_obs**: observed total number of spawners
- **n_age3_obs**: observed frequencies of age-3 spawners
- **n_age4_obs**: observed frequencies of age-4 spawners
- **n_age5_obs**: observed frequencies of age-5 spawners
- **n_H_obs**: observed frequencies of hatchery-origin spawners (all 0 because hatchery-reared chum are not marked)
- **n_W_obs**: observed frequencies of natural-origin spawners (all 0 because hatchery-reared chum are not marked)
- **fit_p_HOS**: logical (T/F) indicating whether the model should estimate proportion hatchery-origin spawners (all FALSE because hatchery-reared chum are not marked)
- **B_take_obs**: number of adults taken for hatchery broodstock (assumed 0)
- **F_rate**: total commercial fisheries harvest rate (proportion) 
- **SST_spr_mean**: average spring (Apr-Jun) sea surface temperature along the Washington coast in the year following spawning
- **upwell_spr**: average spring (Apr-Jun) CUTI upwelling index along the Washington coast in the year following spawning
- **hatch_release**: number of hatchery-reared juvenile chum released in Willapa Bay in the year following spawning
- **Spartina_acre**: acreage of *Spartina alterniflora* in Willapa Bay in the year following spawning
- **OGSI_mean**: average watershed-scale old-growth structural index in the year of spawning
- **OGSI_mean_sbuff250**: average old-growth structural index within 250-m stream buffers in the year of spawning
- **OGSI_mean_sbuff100**: average old-growth structural index within 100-m stream buffers in the year of spawning
- **OGSI_mean_sbuff50**: average old-growth structural index within 50-m stream buffers in the year of spawning

#### WA_SST_projections.csv
This file contains sea surface temperature projections for the Washington coast. The following columns are included:

- **year**: calendar year of spawning
- **model**: global climate model from which sea surface temperature projections are derived
- **SST_spr_mean**: average projected spring (Apr-Jun) sea surface temperature along the Washington coast in the year following spawning

#### WA_upwell_projections.csv
This file contains coastal upwelling projections for the Washington coast. The following columns are included:

- **year**: calendar year of spawning
- **upwell_spr**: average projected spring (Apr-Jun) CUTI upwelling index along the Washington coast in the year following spawning

#### Ellsworth_OGSI_scenarios.csv
This file contains data representing two hypothetical future forest management scenarios for the Ellsworth Creek watershed. The following columns are included:

- **year**: calendar year of spawning
- **scenario**: hypothethical future forest management scenario; either "EFM" (ecological forest management) or "harvest" (industrial timber harvest)
- **OGSI_mean**: hypothetical future average watershed-scale old-growth structural index in the year of spawning



## R script files

#### fit_plot_salmonIPM.R
Code for fitting integrated population models using [salmonIPM](https://github.com/ebuhle/salmonIPM) and plotting the results.

#### func_PQE.R
Functions for calculating probabilities of quasi-extinction (PQE) from future population projections.
