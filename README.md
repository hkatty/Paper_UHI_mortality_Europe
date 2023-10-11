# Urban heat islands' impact on human mortality risk in European cities

This repository contains core code associated with the paper:

Huang et al. 2023, Economic valuation of temperature-related mortality attributed to urban heat islands in European cities.

## Data required

Please note that to run the analysis scripts, additional data need to be downloaded from their respective sources as listed below:

- UrbClim temperature and mask data from Copernicus Climate Change Service (https://doi.org/10.24381/cds.c6459d3a)
- Temperature-mortality relationships from Masselot et al. 2023 (https://doi.org/10.5281/zenodo.7672108)
- Elevation map from MERIT DEM (http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM)
- Population density from NASA SEDAC (https://doi.org/10.7927/H49C6VHW)
- Land imperviousness data from Copernicus Land Monitoring Service (https://land.copernicus.eu/pan-european/high-resolution-layers/imperviousness)
- Köppen-Geiger climate classification by Beck et al. 2018 (https://doi.org/10.1038/sdata.2018.214) from GloH2O (https://www.gloh2o.org/koppen)
- Eurostat data from https://ec.europa.eu/eurostat, including:
  - NUTS 3 level:
    - Annual mortality by age and sex (data code: demo_r_magec3)
    - Population by age and sex (data code: demo_r_pjangrp3)
    - Annual population structures by age (data code: demo_r_pjanind3)
  - NUTS 2 level:
    - Life expectancy by age and sex (data code: demo_r_mlifexp)
  - City level:
    - Average area of living accommodation in m<sup>2</sup> per person (data code: urb_clivcon SA1022V)
    - Average annual rent for housing per m<sup>2</sup> (data code: urb_clivcon SA1049V)
    - Cost of monthly public transport ticket (data code: urb_ctran TT1080V)
- Data on mortality associated with PM2.5 from the supplementary materials associated with Khomenko et al. 2021 (https://doi.org/10.1016/S2542-5196(20)30272-2)
- Data on mortality associated with ozone from the European Environment Agency (https://www.eea.europa.eu/data-and-maps/data/air-quality-health-risk-assessments, permalink: https://www.eea.europa.eu/ds\_resolveuid/86ef37b3bf844299b978867c86cf99e7)

Data generated from the study can be found on Zenodo at https://doi.org/10.5281/zenodo.7986841.

## Code organisation

Code in this repository is grouped as described below and are currently designed to be modified and executed individually as needed. Scripts in subsequent groups may be dependent on outputs from a script in the previous group, though the data outputs on Zenodo as specified above (https://doi.org/10.5281/zenodo.7986841) can also be used in place of certain processing steps. Unless otherwise noted, the scripts within each group are not dependent on each other.
- Files starting with "`00`" are pre-processing scripts. Please note that the bash script `000_prep_data_urbclim.sh` would need to run first to prepare the UrbClim temperature data.
- Files starting with "`01`" are attribution scripts applying the exposure-response relationships to the temperature data. Two scripts are provided:
  - `01_attribute.R` is for the best estimate exposure-response relationships and outputs timeseries of spatial maps
  - `01_attribute_simulations.R` is for the Monte Carlo simulated ensemble and outputs only the timeseries of difference between urban and rural averages
  - Attribution outputs are in attributable fraction.
- Files starting with "`02`" are core analysis scripts, translating attributable fractions into attributable deaths, aggregating across age groups, and performing spatial-temporal analyses.
  - The "base analysis" scripts (`02_base_analysis.py` and `02_base_analysis_simulated.py`) perform aggregations mostly on the temporal dimension, leaving the spatial dimension in the output (if present in the input).
  - The "timeseries" scripts (`02_urbanrural_avg_timeseries.py` and `02-2_timeseries_temporal_analysis.py`) are spatially aggregated into urban and rural averages, as well as the difference between the two. Note that `02-2_timeseries_temporal_analysis.py` is dependent on the outputs of `02_urbanrural_avg_timeseries.py`.
- Files starting with "`03`" are additional analysis scripts, providing outputs quoted in the paper. Supplementary tables S4-S11, listing the risks for each city, are outputs of `03_confidence_interval.py`.
- Files starting with "`Fig`" are used to produce the figures in the main text of the paper, along with any additional associated analyses.

## References

Beck, H. E. et al. (2018) Present and future Köppen-Geiger climate classification maps at 1-km resolution. *Scientific Data*, **5**, 1–12. https://doi.org/10.1038/sdata.2018.214.

Khomenko, S. et al. (2021) Premature mortality due to air pollution in European cities: a health impact assessment. *The Lancet Planetary Health*, **5**, e121–e134. https://doi.org/10.1016/S2542-5196(20)30272-2.

Masselot, P. et al. (2023) Excess mortality attributed to heat and cold: a health impact assessment study in 854 cities in Europe. *The Lancet Planetary Health*, **7**, e271–e281. https://doi.org/10.1016/S2542-5196(23)00023-2.

