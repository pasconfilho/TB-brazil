import pickle, warnings
import numpy as np, pandas as pd, geopandas as gpd
import libpysal
from esda.moran import Moran_Rate, Moran_Local_Rate
from esda.smoothing import Spatial_Empirical_Bayes, Empirical_Bayes

warnings.filterwarnings("ignore")
R = pickle.load(open("/home/claude/work/br/res_br.pkl", "rb"))
mun = R["mun"].copy()
g = gpd.read_file("/home/claude/work/br/ext/mun_br_simpl.gpkg")
g["mun"] = g.cod7.str[:6]
g = g.merge(mun.reset_index().drop(columns=["cod7"]), on="mun", how="left")
g = g[g["pop2022"].notna()].reset_index(drop=True)
print("municipios no mapa", len(g))
w = libpysal.weights.Queen.from_dataframe(g, use_index=False)
print("islands", len(w.islands))
# reforça ilhas com vizinho mais próximo
if w.islands:
    knn = libpysal.weights.KNN.from_dataframe(g.set_geometry(g.geometry.centroid), k=1)
    w = libpysal.weights.attach_islands(w, knn)
w.transform = "r"
py = (g["pop2022"] * 11).values.astype(float)
out = {}
for nm, col in [("all", "cases"), ("np", "cases_np")]:
    e = g[col].values.astype(float)
    g[f"eb_{nm}"] = Empirical_Bayes(e, py).r.flatten() * 1e5
    g[f"seb_{nm}"] = Spatial_Empirical_Bayes(e, py, w).r.flatten() * 1e5
    mi = Moran_Rate(e, py, w, permutations=999)
    lm = Moran_Local_Rate(e, py, w, permutations=999, seed=2026)
    q = lm.q; sig = lm.p_sim < 0.05
    lab = np.where(~sig, "Not significant", np.select([q == 1, q == 2, q == 3, q == 4], ["High-high", "Low-high", "Low-low", "High-low"], default="Not significant"))
    g[f"lisa_{nm}"] = lab
    out[nm] = dict(I=mi.I, p=mi.p_sim, counts=pd.Series(lab).value_counts().to_dict())
    hh = g[g[f"lisa_{nm}"] == "High-high"]
    out[nm]["hh_pop"] = hh["pop2022"].sum(); out[nm]["hh_cases"] = hh[col].sum(); out[nm]["hh_n"] = len(hh)
    out[nm]["hh_uf"] = hh.cod7.str[:2].value_counts().to_dict()
print(out)
R["spatial"] = out
g.drop(columns="geometry").to_parquet("/home/claude/work/br/mun_spatial.parquet")
g[["mun", "cod7", "name", "geometry", "eb_all", "eb_np", "seb_all", "seb_np", "lisa_all", "lisa_np", "pct_ppl", "pop2022", "cases", "cases_np", "rgi"]].to_file("/home/claude/work/br/mun_spatial.gpkg", driver="GPKG")
pickle.dump(R, open("/home/claude/work/br/res_br.pkl", "wb"))
