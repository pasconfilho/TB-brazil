# Pulmonary tuberculosis in Brazil before, during and after the COVID-19 pandemic, 2015-2025

Code and aggregated data for the article *Pulmonary tuberculosis notifications, deaths and treatment outcomes in Brazil before, during and after the COVID-19 pandemic, 2015-2025: a nationwide time-series and spatial analysis* (manuscript under review).

A single commented R script, `R/tb_brazil_analysis.R`, reproduces every analysis in the paper from public data:

| Step | What it produces | Manuscript |
|---|---|---|
| Case selection | flow of records, checks against official tabulations | Figure S1 |
| Annual series | crude and age-standardised rates, Prais-Winsten annual percent change | Table S2 |
| Monthly series | segmented quasi-Poisson regression (Newey-West SE); counterfactual with simulated 95% prediction intervals, for Brazil, regions and states | Table 3, Figures 1-2, Table S3, Figure S7 |
| Sensitivity | flat and shorter baselines, negative binomial, window to 2024, states excluded, alternative case definitions, 12-month moving-block bootstrap | Table S1 |
| Rate ratios | 2022-2025 and 2020-2021 vs 2015-2019 by sex, age, race/skin colour (census-interpolated denominators, extended to 2025 as a sensitivity analysis, with rate differences), region and state; by age with laboratory-confirmed and with all new pulmonary cases | Figure 3, Table S5, Table S7 |
| Profile | characteristics by period, laboratory indicators, HIV, entry types (new, relapse, return after loss to follow-up), completeness of the PDL field | Table 1, Table 2, Figure S2 |
| Outcomes | cure, loss to follow-up and death by year, period and subgroup, with closed or all cases as denominator; modified Poisson model for loss to follow-up | Table 4, Figure 5, Table S4 |
| Spatial | empirical Bayes rates, Moran's I for rates, LISA with FDR control, immediate geographic regions, funnel plot with limits adjusted for overdispersion | Figure 4, Figure S4, Table S6 |
| Mortality | deaths with TB as underlying cause (SIM), observed/expected, age-standardised mortality and an ageing-adjusted counterfactual | Figure S5, Table S1 |

The figures in the manuscript were drawn with Python (matplotlib); those scripts are in `python/` (see `python/README.md`). The R script draws simpler versions of the main figures.

## Data

**Individual tuberculosis records (Sinan).** Not included in the repository because of their size (about 25 MB compressed). Set `DOWNLOAD_SINAN <- TRUE` at the top of the script to download the national files for 2014-2025 from DATASUS with [microdatasus](https://github.com/rfsaldanha/microdatasus) (development version, which supports Sinan tuberculosis), or place a CSV with the same variables at `data-raw/tb_sinan_BR_2014_2025.csv.gz`. Files for 2020-2025 are preliminary and DATASUS updates them, so counts may differ slightly from the paper (files downloaded on 23 November 2026). The script prints two validation counts against official tabulations.

**Aggregated files in `data/`** (all public):

| File | Content | Source |
|---|---|---|
| `pop_uf_ibge_projections_2024.csv` | Population by state, age group and sex, 2015-2025 | IBGE, population projections, 2024 revision |
| `race_uf_census2010.csv`, `race_uf_census2022.csv` | Population by state and race/skin colour | IBGE censuses 2010 (SIDRA table 136) and 2022 (SIDRA table 9605) |
| `municipal_population_census2022.csv` | Population by municipality and immediate geographic region | IBGE, 2022 census; regional division 2017 |
| `immediate_regions.csv` | Names of the 510 immediate geographic regions | IBGE |
| `sim_tb_deaths_monthly_uf.csv` | Deaths with underlying cause A15-A19, by month and state of residence | Mortality Information System (SIM), DATASUS TabNet |
| `sim_tb_deaths_by_icd_year.csv` | Same deaths by three-character ICD-10 code and year | SIM, TabNet |
| `sim_tb_deaths_by_age_year.csv` | Same deaths by age group and year | SIM, TabNet |

Municipal boundaries are downloaded at run time from [geodata-br](https://github.com/tbrugz/geodata-br) (IBGE meshes). See `data-raw/README.md` for the exact queries.

## How to run

```r
# R >= 4.2
install.packages(c("data.table", "MASS", "sandwich", "lmtest", "sf", "spdep", "ggplot2", "scales", "R.utils"))
# optional, only to download Sinan files:
# remotes::install_github("rfsaldanha/microdatasus", ref = "dev")
```

Then, from the repository root:

```bash
Rscript R/tb_brazil_analysis.R
```

The run takes about two minutes on a laptop once the data are in place (the Sinan download takes much longer). Results go to `output/tables/*.csv`, `output/figures/*.png` and `output/results.rds`; `output/sessionInfo.txt` records package versions. Seeds are fixed, so simulation-based intervals are reproducible up to Monte Carlo error.

## Definitions used

- **Case**: pulmonary form only (`FORMA == 1`), new case (`TRATAMENTO == 1`), laboratory-confirmed by positive smear, positive culture or *M. tuberculosis* detected on the rapid molecular test; residence in Brazil; date of diagnosis 2015-2025.
- **Periods**: pre-pandemic 2015-2019 (monthly baseline January 2015-February 2020), pandemic March 2020-December 2021, post-pandemic 2022-2025.
- **PDL**: persons deprived of liberty (`POP_LIBER == 1`).
- **Loss to follow-up**: outcome codes 2 and 10 (including primary loss) among closed records (codes 1-10).

## Ethics

Public, anonymised secondary data (Brazilian National Health Council Resolution 510/2016); no ethics approval required.

## Licence and citation

Code: MIT licence (see `LICENSE`). Data in `data/` come from IBGE and DATASUS and remain subject to their terms of use. If you use this code, please cite the article (citation to be added on publication).
