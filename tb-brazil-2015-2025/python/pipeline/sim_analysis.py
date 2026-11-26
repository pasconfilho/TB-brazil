import re, json, pickle, warnings, numpy as np, pandas as pd, statsmodels.api as sm, statsmodels.formula.api as smf
warnings.filterwarnings("ignore"); rng = np.random.default_rng(2026)
MES = {m: i + 1 for i, m in enumerate(["Janeiro", "Fevereiro", "Março", "Abril", "Maio", "Junho", "Julho", "Agosto", "Setembro", "Outubro", "Novembro", "Dezembro"])}
lines = open("tabnet_sim_tb_mes_uf.txt").read().splitlines()
hdr = [l for l in lines if l.startswith('"Ano/mês')][0].replace('"', '').split(";")
rows = []
for l in lines:
    m = re.match(r'^"[\.\s]*([A-Za-zçÇ]+)/(\d{4})";(.*)$', l)
    if m and m.group(1) in MES:
        v = [0 if x == "-" else int(x) for x in m.group(3).split(";")]
        rows.append([int(m.group(2)), MES[m.group(1)]] + v)
D = pd.DataFrame(rows, columns=["ano", "mes"] + hdr[1:]); D = D[D.ano.between(2015, 2025)].reset_index(drop=True)
assert len(D) == 132, len(D)
D.to_csv("sim_tb_mensal_uf.csv", index=False)
print(D.groupby("ano").Total.sum().to_dict())
R = pickle.load(open("/home/claude/work/br/res_br.pkl", "rb")); mb = R["monthly"].copy()
mb["y"] = D.Total.values; t_p, t_q = 62, 84
def cf(m, nsim=10000, start=0, trend=True, post_end=None):
    pre = m[(m.t >= start) & (m.t < t_p)]; cols = (["t"] if trend else []) + ["s1", "c1", "s2", "c2"]; XC = ["Intercept"] + cols
    mp = smf.glm("y ~ " + " + ".join(cols), data=pre, family=sm.families.Poisson(), offset=pre.off).fit(scale="X2")
    b = mp.params[XC].values; V = mp.cov_params().loc[XC, XC].values; ph = max(mp.scale, 1.0001)
    X = sm.add_constant(m[cols], has_constant="add").values; MU = np.exp(rng.multivariate_normal(b, V, nsim) @ X.T + m.off.values)
    YS = rng.poisson(rng.gamma(MU / (ph - 1), ph - 1)); e = np.exp(X @ b + m.off.values); out = {}
    pend = post_end or len(m)
    sels = [("Pandemic", (m.t >= t_p) & (m.t < t_q)), ("Post", (m.t >= t_q) & (m.t < pend)), ("Cumulative", (m.t >= t_p) & (m.t < pend))]
    sels += [(str(y), (m.t >= max(t_p, (y - 2015) * 12)) & (m.t < (y - 2014) * 12)) for y in range(2020, 2026)]
    for lab, sel in sels:
        ii = np.where(sel.values)[0]; o = m.y.values[ii].sum(); s_ = YS[:, ii].sum(axis=1)
        out[lab] = dict(obs=int(o), exp=e[ii].sum(), OE=o / e[ii].sum(), lo=o / np.percentile(s_, 97.5), hi=o / np.percentile(s_, 2.5),
                        diff=o - e[ii].sum(), dlo=o - np.percentile(s_, 97.5), dhi=o - np.percentile(s_, 2.5))
    return out, e, np.percentile(YS, 2.5, axis=0), np.percentile(YS, 97.5, axis=0)
o, e, lo, hi = cf(mb); mb["exp_d"], mb["lo_d"], mb["hi_d"] = e, lo, hi
for k, v in o.items(): print(k, {a: round(float(b), 3) for a, b in v.items()})
sens = {"flat_2017": cf(mb, start=26, trend=False)[0], "flat_2018": cf(mb, start=38, trend=False)[0], "to2024": cf(mb, post_end=120)[0]}
for k, v in sens.items(): print(k, {p: round(v[p]["OE"], 3) for p in ["Pandemic", "Post", "Cumulative"]}, {p: (round(v[p]["lo"], 3), round(v[p]["hi"], 3)) for p in ["Post"]})
# taxa de mortalidade e razão óbitos/casos
ann = D.groupby("ano").Total.sum(); pop = R["pop_br"]
tb = pd.read_parquet("/home/claude/work/br/tb_br.parquet", columns=["DT_DIAG", "TRATAMENTO", "ID_MN_RESI"])
tb["ano"] = pd.to_datetime(tb.DT_DIAG, errors="coerce").dt.year
newall = tb[(tb.TRATAMENTO == "1") & tb.ID_MN_RESI.notna()].groupby("ano").size().reindex(range(2015, 2026))
A = pd.DataFrame({"deaths": ann, "mort_rate": ann / pop.values * 1e5, "new_cases_all_forms": newall.values})
A["deaths_per_100_cases"] = A.deaths / A.new_cases_all_forms * 100
print(A.round(2))
cat = json.load(open("tabnet_sim_tb_extra.json"))
json.dump({"cf": o, "sens": sens, "annual": A.reset_index().to_dict("list")}, open("sim_results.json", "w"), default=float, indent=1)
mb[["y", "exp_d", "lo_d", "hi_d"]].to_csv("sim_monthly_cf.csv")
