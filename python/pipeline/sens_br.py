"""Análises de sensibilidade do contrafactual nacional."""
import pickle, warnings
import numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf

warnings.filterwarnings("ignore")
rng = np.random.default_rng(2026)
R = pickle.load(open("/home/claude/work/br/res_br.pkl", "rb"))
cn = pd.read_parquet("/home/claude/work/br/cn_br.parquet")
mb = R["monthly"]
months = pd.period_range("2015-01", "2025-12", freq="M")
t_p = list(months).index(pd.Period("2020-03", "M")); t_q = list(months).index(pd.Period("2022-01", "M"))
t_25 = list(months).index(pd.Period("2025-01", "M"))
POP = mb[["pop", "t", "dias", "s1", "c1", "s2", "c2", "off"]].copy()

def frame(sub):
    m = POP.copy(); m["y"] = sub.groupby("mes").size().reindex(months, fill_value=0).values; return m

def cf(m, trend=True, start=0, family="qp", nsim=10000, post_end=None):
    pre = m[(m.t >= start) & (m.t < t_p)]
    rhs = ("t + " if trend else "") + "s1 + c1 + s2 + c2"
    cols = (["t"] if trend else []) + ["s1", "c1", "s2", "c2"]; XC = ["Intercept"] + cols
    if family == "nb":
        mp = smf.glm(f"y ~ {rhs}", data=pre, family=sm.families.NegativeBinomial(alpha=0.01), offset=pre.off).fit()
        # estima alpha por momentos
        mu = mp.fittedvalues; a = max(((pre.y - mu) ** 2 - mu).sum() / (mu ** 2).sum(), 1e-4)
        mp = smf.glm(f"y ~ {rhs}", data=pre, family=sm.families.NegativeBinomial(alpha=a), offset=pre.off).fit()
        b = mp.params[XC].values; V = mp.cov_params().loc[XC, XC].values
        X = sm.add_constant(m[cols], has_constant="add").values; MU = np.exp(rng.multivariate_normal(b, V, nsim) @ X.T + m.off.values)
        YS = rng.poisson(rng.gamma(1 / a, MU * a))
    else:
        mp = smf.glm(f"y ~ {rhs}", data=pre, family=sm.families.Poisson(), offset=pre.off).fit(scale="X2")
        b = mp.params[XC].values; V = mp.cov_params().loc[XC, XC].values; ph = max(mp.scale, 1.0001)
        X = sm.add_constant(m[cols], has_constant="add").values; MU = np.exp(rng.multivariate_normal(b, V, nsim) @ X.T + m.off.values)
        YS = rng.poisson(rng.gamma(MU / (ph - 1), ph - 1))
    e = np.exp(X @ b + m.off.values); out = {}
    pend = post_end if post_end is not None else len(m)
    for lab, sel in [("Pandemic", (m.t >= t_p) & (m.t < t_q)), ("Post", (m.t >= t_q) & (m.t < pend)), ("Cumulative", (m.t >= t_p) & (m.t < pend))]:
        ii = np.where(sel.values)[0]; o = m.y.values[ii].sum(); ee = e[ii].sum(); s_ = YS[:, ii].sum(axis=1)
        out[lab] = dict(obs=o, exp=ee, OE=o / ee, OE_lo=o / np.percentile(s_, 97.5), OE_hi=o / np.percentile(s_, 2.5), diff=o - ee)
    return out

SUBS = {"All": cn, "Excluding PDL": cn[cn.ppl != 1], "PDL": cn[cn.ppl == 1]}
SPECS = [("Main: Poisson-quasi, trend + seasonality, Jan 2015-Feb 2020", dict()),
         ("Shorter baseline (Mar 2017-Feb 2020), trend", dict(start=26)),
         ("No trend, baseline Mar 2017-Feb 2020 (flat rate)", dict(trend=False, start=26)),
         ("No trend, baseline Mar 2018-Feb 2020 (flat rate)", dict(trend=False, start=38)),
         ("Negative binomial, trend + seasonality", dict(family="nb")),
         ("Post-pandemic window 2022-2024 (2025 data preliminary)", dict(post_end=t_25))]
rows = []
for sn, sub in SUBS.items():
    m = frame(sub)
    for lab, kw in SPECS:
        o = cf(m, **kw)
        for per in ["Pandemic", "Post", "Cumulative"]:
            rows.append(dict(series=sn, spec=lab, period=per, **o[per]))
S = pd.DataFrame(rows)

# excluindo UFs com qualidade de dados questionável (PB: série com descontinuidade)
for drop in [["PB"], ["PB", "PE", "DF"]]:
    sub = cn[(~cn.uf.isin(drop)) & (cn.ppl != 1)]
    # população sem as UFs
    P = pd.read_csv("/home/claude/work/br/ext/pop_uf.csv"); P.columns = [int(c) if c.isdigit() else c for c in P.columns]
    yrs = list(range(2015, 2026))
    pbr = P[(P.uf == "BR") & P.grupo.isin(["Homens", "Mulheres"])][yrs].sum()
    pdr = P[(P.uf.isin(drop)) & P.grupo.isin(["Homens", "Mulheres"])][yrs].sum()
    ratio = ((pbr - pdr) / pbr).values.astype(float)
    m = frame(sub); mid = np.array([pd.Period(f"{y}-07", "M").ordinal for y in yrs]); idx = np.array([p.ordinal for p in months])
    m["off"] = m["off"] + np.log(np.interp(idx, mid, ratio))
    o = cf(m)
    for per in ["Pandemic", "Post", "Cumulative"]:
        rows.append(dict(series="Excluding PDL", spec=f"Excluding {', '.join(drop)}", period=per, **o[per]))
S = pd.DataFrame(rows)
pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 70)
print(S.assign(r=lambda x: x.OE.round(3).astype(str) + " (" + x.OE_lo.round(3).astype(str) + "-" + x.OE_hi.round(3).astype(str) + ")")
      .pivot_table(index=["series", "spec"], columns="period", values="r", aggfunc="first")[["Pandemic", "Post", "Cumulative"]])
R["sens"] = S
pickle.dump(R, open("/home/claude/work/br/res_br.pkl", "wb"))
