# Python scripts used for the manuscript figures

The R script in `R/` is the reference implementation: it reproduces every number in the paper from the public data. The figures printed in the manuscript were drawn with matplotlib, and the scripts are kept here for transparency.

- `pipeline/`: the original Python analysis (`analysis_br.py`, `sens_br.py`, `spatial_br.py`, `revisions_br.py`, `sim_analysis.py`, `quality.py`). It produces the intermediate objects the figure scripts read (`res_br.pkl`, `mun_spatial.gpkg`, `sim_results.json`). `sync_from_R.py` replaces the annual percent changes and spatial results with those of the R script, which uses the exact Prais-Winsten estimator and unsimplified municipal boundaries.
- `figures/`: one script per figure (Figure 1: `fig1_v2b.py`; Figure 2: `fig2_v2.py`; Figures 3-5 and S2: `figures_br.py`; Figure S1: `figures_supp_br.py`; Figure S3: `fig_qualidade.py`; Figure S4 and Table S6: `fig_mun_br.py`; Figure S5: `fig_sim_v2.py`; Figure S6: `fig1_v2.py`).

These are working scripts: paths point to the authors' working directory and some comments are in Portuguese. To run them, set the paths at the top of each file to your copies of the Sinan file and of the files in `data/` (the race files in `data/` use English column names; the original scripts expect the IBGE Portuguese names Branca, Preta, Amarela, Parda, Indígena). Requirements: Python 3.10+, pandas, numpy, statsmodels, scipy, geopandas, libpysal, esda, matplotlib.
