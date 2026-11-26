"""Análises adicionais solicitadas na revisão (v2.0)."""
import pickle, warnings, numpy as np, pandas as pd, statsmodels.api as sm, statsmodels.formula.api as smf
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore"); rng = np.random.default_rng(2026)
R = pickle.load(open("res_br.pkl", "rb")); cn = pd.read_parquet("cn_br.parquet"); mb = R["monthly"].copy(); mb["t"] = np.arange(132)
months = pd.period_range("2015-01", "2025-12", freq="M"); t_p, t_q = 62, 84; YEARS = list(range(2015, 2026))
OUT = {}
# ---------- 1. denominadores por raça interpolados 2010-2022 por UF ----------
UF = {'11': 'RO', '12': 'AC', '13': 'AM', '14': 'RR', '15': 'PA', '16': 'AP', '17': 'TO', '21': 'MA', '22': 'PI', '23': 'CE', '24': 'RN', '25': 'PB', '26': 'PE',
      '27': 'AL', '28': 'SE', '29': 'BA', '31': 'MG', '32': 'ES', '33': 'RJ', '35': 'SP', '41': 'PR', '42': 'SC', '43': 'RS', '50': 'MS', '51': 'MT', '52': 'GO', '53': 'DF'}
RC = ["Branca", "Preta", "Amarela", "Parda", "Indígena"]; EN = dict(zip(RC, ["White", "Black", "Asian", "Brown (pardo)", "Indigenous"]))
r10 = pd.read_csv("ext/race_uf_2010.csv", dtype={"cuf": str}).set_index("cuf")[RC]; r22 = pd.read_csv("ext/race_uf_2022.csv", dtype={"cuf": str}).set_index("cuf")[RC]
s10 = r10.div(r10.sum(1), axis=0); s22 = r22.div(r22.sum(1), axis=0); popuf = R["pop_uf"]
racepop = {EN[c]: pd.Series(0.0, index=YEARS) for c in RC}
for cuf, u in UF.items():
    for y in YEARS:
        w = min((y - 2010) / 12, 1.0); sh = s10.loc[cuf] * (1 - w) + s22.loc[cuf] * w   # interpola 2010->2022; mantém 2022 depois
        for c in RC: racepop[EN[c]][y] += sh[c] * popuf.loc[y, u]
racepop = pd.DataFrame(racepop); OUT["racepop"] = racepop
def rt(a, pa, b, pb):
    r = (a / pa) / (b / pb); se = np.sqrt(1 / a + 1 / b); return r, r * np.exp(-1.96 * se), r * np.exp(1.96 * se)
rows = []
for rc in racepop.columns:
    n = cn[cn.raca == rc].groupby("ano").size().reindex(YEARS, fill_value=0); p = racepop[rc]
    a, pa = n.loc[2022:2025].sum(), p.loc[2022:2025].sum(); b, pb = n.loc[2015:2019].sum(), p.loc[2015:2019].sum()
    c_, pc = n.loc[2020:2021].sum(), p.loc[2020:2021].sum()
    r_, l_, h_ = rt(a, pa, b, pb); rp, lp, hp = rt(c_, pc, b, pb)
    rows.append(dict(race=rc, n=n.sum(), rate=n.sum() / p.sum() * 1e5, rate_pre=b / pb * 1e5, rate_post=a / pa * 1e5, RR=r_, lo=l_, hi=h_, RRpand=rp, RRpand_lo=lp, RRpand_hi=hp))
race2 = pd.DataFrame(rows).set_index("race"); wr = race2.loc["White", "rate"]
race2["RRw"] = race2.rate / wr
for rc in race2.index:
    k = race2.loc[rc, "n"]; kb = race2.loc["White", "n"]; se = np.sqrt(1 / k + 1 / kb)
    race2.loc[rc, "RRw_lo"] = race2.loc[rc, "RRw"] * np.exp(-1.96 * se); race2.loc[rc, "RRw_hi"] = race2.loc[rc, "RRw"] * np.exp(1.96 * se)
OUT["race2"] = race2; print(race2.round(3))
# ---------- 2. contrafactual: definições alternativas ----------
d = pd.read_parquet("tb_br.parquet", columns=["DT_DIAG", "TRATAMENTO", "FORMA", "ID_MN_RESI", "BACILOSC_E", "CULTURA_ES", "TEST_MOLEC"])
d["mes"] = pd.to_datetime(d.DT_DIAG, errors="coerce").dt.to_period("M"); d = d[d.ID_MN_RESI.fillna("").str[:2].isin(UF.keys()) & (d.FORMA == "1")]
d["lab"] = (d.BACILOSC_E == "1") | (d.CULTURA_ES == "1") | d.TEST_MOLEC.isin(["1", "2"])
def cf(y, trend=True, start=0, nsim=10000, end=132):
    m = mb.copy(); m["y"] = y; pre = m[(m.t >= start) & (m.t < t_p)]; cols = (["t"] if trend else []) + ["s1", "c1", "s2", "c2"]; XC = ["Intercept"] + cols
    mp = smf.glm("y ~ " + " + ".join(cols), data=pre, family=sm.families.Poisson(), offset=pre.off).fit(scale="X2")
    b = mp.params[XC].values; V = mp.cov_params().loc[XC, XC].values; ph = max(mp.scale, 1.0001); X = sm.add_constant(m[cols], has_constant="add").values
    MU = np.exp(rng.multivariate_normal(b, V, nsim) @ X.T + m.off.values); YS = rng.poisson(rng.gamma(MU / (ph - 1), ph - 1)); e = np.exp(X @ b + m.off.values)
    o = {}
    for lab, sel in [("Pandemic", (m.t >= t_p) & (m.t < t_q)), ("Post", (m.t >= t_q) & (m.t < end)), ("Cumulative", (m.t >= t_p) & (m.t < end))]:
        ii = np.where(sel.values)[0]; ob = m.y.values[ii].sum(); s_ = YS[:, ii].sum(axis=1)
        o[lab] = dict(OE=ob / e[ii].sum(), lo=ob / np.percentile(s_, 97.5), hi=ob / np.percentile(s_, 2.5))
    return o, mp, e
alt = {"All new pulmonary cases (any confirmation)": d[d.TRATAMENTO == "1"], "All laboratory-confirmed pulmonary cases (new and retreatment)": d[d.lab]}
for k, sub in alt.items():
    y = sub.groupby("mes").size().reindex(months, fill_value=0).values; o, _, _ = cf(y); OUT[k] = o
    print(k, {p: (round(v["OE"], 3), round(v["lo"], 3), round(v["hi"], 3)) for p, v in o.items()})
# ---------- 3. autocorrelação dos resíduos pré-pandemia e bootstrap em blocos ----------
deaths = pd.read_csv("sim/sim_tb_mensal_uf.csv").Total.values
SER = {"All": mb.y.values, "Excluding PDL": mb.y_nonppl.values, "PDL": mb.y_ppl.values, "Deaths": deaths}
ac = {}
for k, y in SER.items():
    o, mp, e = cf(y, nsim=2000); r = mp.resid_pearson
    lb = acorr_ljungbox(r, lags=[12], return_df=True); a1 = np.corrcoef(r[:-1], r[1:])[0, 1]
    # bootstrap em blocos (12 meses) dos resíduos de log(y/e) no pré-período, aplicado ao esperado projetado
    m = mb.copy(); m["y"] = y; pre = m.t < t_p; lr = np.log(np.maximum(y[pre], 0.5) / e[pre]); lr = lr - lr.mean()
    B = 4000; L = 12; nb = int(np.ceil(132 / L)); sims = np.empty((B, 132))
    Vsim = rng.multivariate_normal(mp.params.values, mp.cov_params().values, B); X = sm.add_constant(m[["t", "s1", "c1", "s2", "c2"]]).values
    for bb in range(B):
        st = rng.integers(0, len(lr) - L + 1, nb); eps = np.concatenate([lr[s:s + L] for s in st])[:132]
        sims[bb] = np.exp(X @ Vsim[bb] + m.off.values + eps)
    post = (m.t >= t_q).values; ob = y[post].sum(); sp = sims[:, post].sum(1)
    ac[k] = dict(rho1=a1, LB12_p=float(lb.lb_pvalue.iloc[0]), OE_post=o["Post"]["OE"], PI_param=(o["Post"]["lo"], o["Post"]["hi"]),
                 PI_block=(ob / np.percentile(sp, 97.5), ob / np.percentile(sp, 2.5)))
    print(k, {kk: (np.round(v, 3) if not isinstance(v, tuple) else tuple(np.round(v, 3))) for kk, v in ac[k].items()})
OUT["autocorr"] = ac
# ---------- 4. LISA com correção FDR ----------
from esda.moran import Moran_Local_Rate
import geopandas as gpd, libpysal
g = gpd.read_file("mun_spatial.gpkg"); w = libpysal.weights.Queen.from_dataframe(g, use_index=False)
if w.islands: w = libpysal.weights.attach_islands(w, libpysal.weights.KNN.from_dataframe(g.set_geometry(g.geometry.centroid), k=1))
w.transform = "r"
lm = Moran_Local_Rate(g.cases_np.values.astype(float), (g.pop2022 * 11).values.astype(float), w, permutations=9999, seed=2026)
rej = multipletests(lm.p_sim, alpha=0.05, method="fdr_bh")[0]
hh05 = (lm.q == 1) & (lm.p_sim < 0.05); hhf = (lm.q == 1) & rej
tot_np = g.cases_np.sum(); P = g.pop2022.sum()
OUT["lisa"] = dict(hh05_n=int(hh05.sum()), hh05_pop=float(g.pop2022[hh05].sum() / P), hh05_cases=float(g.cases_np[hh05].sum() / tot_np),
                   hhfdr_n=int(hhf.sum()), hhfdr_pop=float(g.pop2022[hhf].sum() / P), hhfdr_cases=float(g.cases_np[hhf].sum() / tot_np))
print(OUT["lisa"])
pickle.dump(OUT, open("revisions_br.pkl", "wb"))
