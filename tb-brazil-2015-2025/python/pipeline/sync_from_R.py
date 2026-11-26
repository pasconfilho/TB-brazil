"""Replace APC and spatial results with the values from the R script (source of truth for v2.0)."""
import pickle, pandas as pd, geopandas as gpd
T = "/home/claude/work/gh/tb-brazil-2015-2025/output/tables/"
R = pickle.load(open("res_br.pkl", "rb"))
ap = pd.read_csv(T + "table_s2_apc.csv").rename(columns={"window": "period"})
R["apc"] = ap[["series", "period", "APC", "lo", "hi", "p", "trend"]]
oa = pd.read_csv(T + "outcomes_apc.csv").set_index("indicator")
R["cure_apc"] = oa.loc["Cure"].to_dict(); R["ltfu_apc"] = oa.loc["Loss to follow-up"].to_dict()
M = pd.read_csv(T + "municipal_results.csv", dtype={"cod7": str}).set_index("cod7")
lab = {"High-High": "High-high", "Low-Low": "Low-low", "High-Low": "High-low", "Low-High": "Low-high", "Not significant": "Not significant"}
for f in ["mun_spatial.gpkg", "mun_spatial.parquet"]:
    g = gpd.read_file(f) if f.endswith("gpkg") else pd.read_parquet(f)
    m = M.reindex(g.cod7)
    assert m.lisa_np.notna().all(), "unmatched municipalities"
    g["lisa_np"] = m.lisa_np.map(lab).values; g["seb_np"] = m.seb_np.values; g["eb_np"] = m.eb_np.values; g["eb_all"] = m.eb_all.values
    g.to_file(f, driver="GPKG") if f.endswith("gpkg") else g.to_parquet(f)
S = pd.read_csv(T + "spatial_summary.csv").value.tolist()
R["spatial"]["np"]["I"] = S[2]; R["spatial"]["all"]["I"] = S[4]
pickle.dump(R, open("res_br.pkl", "wb"))
print(R["apc"].round(2).head(3)); print(pd.Series(g.lisa_np).value_counts())
