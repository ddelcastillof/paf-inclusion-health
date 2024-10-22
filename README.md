# Estimating population-attributable fractions due to social exclusion groups in HICs
The aim of this repo is to estimate global population attributable fractions (PAFs) for inclusion health groups.
 
Data sources at present are not fully referenced. In its current version, the data is obtained from a Google Sheet file which contains:
- Inclusion health group's population sizes
- SMRs for inclusion health groups based on a previous systematic review: [Health of people experiencing co-occurring homelessness, imprisonment, substance use, sex work and/or severe mental illness in high-income countries: a systematic review and meta-analysis](https://doi.org/10.1016/S0140-6736(17)31869-X)
- Other data required for the analyses, e.g. GDB death estimates and population sizes, were obtained from public sources.

## Authors
Originally created by @rwaldridge and modified by me.
Forked from @UCL-Public-Health-Data-Science

## To-do list
- [ ] Fix missing homeless data for Argentina and Belgium (not a priority, probably requires primary data sources).
- [ ] Make the age groups for the SMRs more comparable.
- [ ] Make calendar years more comparable (separate years and decide what to do with the data, check all cause versus all-cause mortality data)
- [x] Add better data on male and female numbers for each country. At the moment, the assumption is based on UK data and 80% M and 20% F
- [ ] Add in population numbers for each country by age and sex
- [x] Model the inclusion health groups overlap using MCM simulation
- [ ] Model the inclusion health distribution by age and sex for each country

## Packages used
- googlesheets4
- metafor
- tidyverse
- skimr
- mc2d
- ggplot2
- patchwork
- rjags
- coda
