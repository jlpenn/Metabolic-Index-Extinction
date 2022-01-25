# Metabolic-Index-Extinction

This repository contains the plotting and analysis scripts and model code for J. Penn, C. Deutsch "Avoiding Ocean Mass Extinction from Climate warming" (2022).

Please contact me at jpenn@princeton.edu if you have any questions :)

Attribution: Please cite as

Penn, J., Deutsch, C. (2022). Analysis scripts and model code for "Avoiding Ocean Mass Extinction from Climate warming".

Penn, J., Deutsch, C. (2022). "Avoiding Ocean Mass Extinction from Climate warming" (under revision).

Description: All code is written in MATLAB (v2018-2019a) and Python 3. Manuscript figures can be reproduced by running the script 'get_6me_Fig.m'. Model output and observational data used for plotting are available in the accompanying .MAT file (6me_output.mat) and can be unpacked using the script 'get_6me_output.m'.

Model code projects global extinction and extirpation risk of marine animal species due to climate warming and ocean oxygen loss by combining species ecophysiological traits and Earth System Model (ESM) simulations of future climate states under different greenhouse gas emissions scenarios. Use 'getCMIP6_clim.m' and 'getCMIP5_clim.m' to create multi-decadal climate anomalies from CMIP input fields. CMIP inputs are annual average climate fields from the greenhouse gas forcing experiments minus pre-industrial time-averaged fields (1850-1900). ‘SixExtinct_spatial.m’ sets up the extinction model by selecting the greenhouse gas emissions scenario, CMIP model, and extinction model parameters. It calls the function 'get6extinct_func.m' which computes each species three-dimensional Metabolic Index and global habitat volume over time, and uses these to project the spatiotemporal patterns of extinction and extirpation. To do this, 'get6extinct_func.m' calls 'getPDF.m' and 'getTmin2.m', which generate the frequency distributions of Metabolic Index traits and species cold limits, respectively. The Metabolic Index is computed using the function 'metab_index.m'. 'SixExtinct_global.m' computes the extinction proportion for the entire global community using the global habitat volumes computed for each species in 'SixExtinct_spatial.m'. Model richness patterns versus latitude are computed using the script 'Model_richness.m', which sets up the calculation (e.g., grid resolution, summation depth) for 'getRichness_func.m'. This latter script uses the observational climatology and the trait distributions to sum the number of model species in each latitude band based on available habitat (using the Metabolic Index). Input data file paths for model code will need to be changed to those found on the local machine.

Three scripts analyze the IUCN redlist data: 'get_iucn_marine.m' computes the fraction of marine species at risk of extinction globally, 'get_iucn_marine_regional.m' does this for different marine regions, and 'get_iucn_threats.m' tallies the number of marine species threatened by different stressors identified by the IUCN.

Data: The data analyzed in this paper are publically available (see Data and materials availability). In particular, input data for the extinction model requires downloading ESM output from the CMIP archive as well as observational climatologies from the World Ocean Atlas. See the README_CMIPdata for a description of CMIP downloading and processing steps and the section 'Materials and Methods' for details, including model, variable, and scenario names.
