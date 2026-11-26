# Raw data and how the files in `data/` were obtained

This folder holds files that are downloaded, not versioned:

- `tb_sinan_BR_2014_2025.csv.gz`: Sinan tuberculosis records, created by `R/tb_brazil_analysis.R` when `DOWNLOAD_SINAN <- TRUE` (with `provenance_sinan.csv`, the list of DATASUS files used).
- `geojs-100-mun.json`: municipal boundaries, downloaded by the script from
  `https://raw.githubusercontent.com/tbrugz/geodata-br/master/geojson/geojs-100-mun.json`.

## Sources of the aggregated files in `data/`

**Population projections** (`pop_uf_ibge_projections_2024.csv`): IBGE, *Projeções da População do Brasil e Unidades da Federação por sexo e idade: 2000-2070*, 2024 revision. Single-year ages were summed into 0-14, 15-24, 25-34, 35-44, 45-54, 55-64 and 65+; rows "Homens" and "Mulheres" are totals by sex; `BR` is Brazil.

**Race/skin colour** (`race_uf_census2010.csv`, `race_uf_census2022.csv`): IBGE SIDRA, table 136 (2010 census) and table 9605 (2022 census), resident population by colour or race (branca, preta, amarela, parda, indígena) and state. In the script, the state composition is interpolated linearly between the two censuses, kept at the 2022 composition from 2023, and applied to the annual state projections.

**Municipal population** (`municipal_population_census2022.csv`): IBGE, 2022 census, resident population by municipality, with the code of the immediate geographic region (IBGE regional division of 2017). `immediate_regions.csv` gives the region names.

**Deaths (SIM)** (`sim_tb_deaths_monthly_uf.csv`, `sim_tb_deaths_by_icd_year.csv`, `sim_tb_deaths_by_age_year.csv`): DATASUS TabNet, *Óbitos por residência*, `http://tabnet.datasus.gov.br/cgi/deftohtm.exe?sim/cnv/obt10uf.def`:

- Row: *Ano/mês do Óbito*; Column: *Unidade da Federação*; Content: *Óbitos por residência*
- Period: 2015-2025 (final up to 2024, 2025 preliminary at the time of access, 24 November 2026)
- Selection: *Categoria CID-10* = A15 to A19
- For the second file, row *Categoria CID-10* and column *Ano do Óbito*.
- For `sim_tb_deaths_by_age_year.csv`, row *Ano do Óbito*, column *Faixa Etária det* and selection *Grupo CID-10* = Tuberculose (A15-A19), accessed 25 November 2026; detailed age bands were summed into 0-14, 15-24, 25-34, 35-44, 45-54, 55-64 and 65+ years (deaths of unknown age kept in a separate column).

These are deaths with tuberculosis as the underlying cause, all forms. Deaths with tuberculosis only as an associated cause are not included.
