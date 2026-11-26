"""Análises nacionais: TB pulmonar com confirmação laboratorial, casos novos (TRATAMENTO=1), Brasil 2015-2025."""
import pickle, warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

warnings.filterwarnings("ignore")
rng = np.random.default_rng(2026)
YEARS = list(range(2015, 2026))
UF = {'11': 'RO', '12': 'AC', '13': 'AM', '14': 'RR', '15': 'PA', '16': 'AP', '17': 'TO', '21': 'MA', '22': 'PI', '23': 'CE',
      '24': 'RN', '25': 'PB', '26': 'PE', '27': 'AL', '28': 'SE', '29': 'BA', '31': 'MG', '32': 'ES', '33': 'RJ', '35': 'SP',
      '41': 'PR', '42': 'SC', '43': 'RS', '50': 'MS', '51': 'MT', '52': 'GO', '53': 'DF'}
REG = {'1': 'North', '2': 'Northeast', '3': 'Southeast', '4': 'South', '5': 'Centre-West'}
AGES = ['0-14', '15-24', '25-34', '35-44', '45-54', '55-64', '65+']
R = {}

# ---------------- dados ----------------
d = pd.read_parquet("/home/claude/work/br/tb_br.parquet")
d["dtd"] = pd.to_datetime(d.DT_DIAG, errors="coerce"); d["ano"] = d.dtd.dt.year; d["mes"] = d.dtd.dt.to_period("M")
d["lab"] = (d.BACILOSC_E == "1") | (d.CULTURA_ES == "1") | d.TEST_MOLEC.isin(["1", "2"])
d["cuf"] = d.ID_MN_RESI.fillna("").str[:2]; d["uf"] = d.cuf.map(UF); d["regiao"] = d.cuf.str[0].map(REG)
d["mun"] = d.ID_MN_RESI.fillna("").str[:6]
def idade(s):
    if not isinstance(s, str) or len(s) != 4 or not s.isdigit(): return np.nan
    u, v = int(s[0]), int(s[1:]); return v if u == 4 else (0 if u in (1, 2, 3) else np.nan)
d["idade"] = d.NU_IDADE_N.map(idade)
d["faixa"] = pd.cut(d.idade, [-1, 14, 24, 34, 44, 54, 64, 200], labels=AGES).astype("object")
d["sexo"] = d.CS_SEXO.map({"M": "Male", "F": "Female"})
d["raca"] = d.CS_RACA.map({"1": "White", "2": "Black", "3": "Asian", "4": "Brown (pardo)", "5": "Indigenous"})
d["se"] = pd.to_numeric(d.SITUA_ENCE, errors="coerce")
d["ppl"] = d.POP_LIBER.map({"1": 1, "2": 0}); d["rua"] = d.POP_RUA.map({"1": 1, "2": 0})
for v, c in [("alcool", "AGRAVALCOO"), ("drogas", "AGRAVDROGA"), ("tabaco", "AGRAVTABAC"), ("diabetes", "AGRAVDIABE")]:
    d[v] = d[c].map({"1": 1, "2": 0})
d["hiv"] = d.HIV.map({"1": "Positive", "2": "Negative", "3": "Pending", "4": "Not done"})

flow = [("TB records in national Sinan files 2014-2025", len(d))]
s = d[d.ano.between(2015, 2025)]; flow.append(("Diagnosed 2015-2025", len(s)))
s = s[s.uf.notna()]; flow.append(("Resident in Brazil (valid municipality)", len(s)))
s = s[s.FORMA == "1"]; flow.append(("Pulmonary only", len(s))); allp = s.copy()
s = s[s.lab]; flow.append(("Laboratory-confirmed", len(s))); notif = s.copy()
s = s[s.TRATAMENTO == "1"]; flow.append(("New cases", len(s))); cn = s.copy()
R["flow"] = flow
R["tabnet_check"] = int((d.ano.between(2015, 2025) & (d.FORMA == "1") & d.lab).sum())
R["ms2025_check"] = int(((d.ano == 2025) & (d.TRATAMENTO == "1")).sum())
cn["per"] = np.where(cn.ano <= 2019, "Pre", np.where(cn.ano <= 2021, "Pandemic", "Post"))

# ---------------- populações ----------------
P = pd.read_csv("/home/claude/work/br/ext/pop_uf.csv"); P.columns = [int(c) if c.isdigit() else c for c in P.columns]
def popget(uf, grupo):
    return P[(P.uf == uf) & (P.grupo == grupo)][YEARS].iloc[0].astype(float)
pop_br = popget("BR", "Homens") + popget("BR", "Mulheres")
pop_uf = pd.DataFrame({u: popget(u, "Homens") + popget(u, "Mulheres") for u in UF.values()})
pop_age_br = pd.DataFrame({a: popget("BR", a) for a in AGES})
pop_sex_br = pd.DataFrame({"Male": popget("BR", "Homens"), "Female": popget("BR", "Mulheres")})
pop_age_uf = {u: pd.DataFrame({a: popget(u, a) for a in AGES}) for u in UF.values()}
W = pop_age_br.loc[2022] / pop_age_br.loc[2022].sum()
RACE = pd.read_csv("/home/claude/work/br/ext/race_uf_2022.csv", dtype={"cuf": str})
race_br = RACE[["Branca", "Preta", "Amarela", "Parda", "Indígena"]].sum()
race_br.index = ["White", "Black", "Asian", "Brown (pardo)", "Indigenous"]
reg_of = {u: REG[c[0]] for c, u in UF.items()}
pop_reg = pd.DataFrame({r: pop_uf[[u for u in UF.values() if reg_of[u] == r]].sum(axis=1) for r in REG.values()})

def by_year(df): return df.groupby("ano").size().reindex(YEARS, fill_value=0)
def asr(df, pa):
    c = df.groupby(["ano", "faixa"]).size().unstack().reindex(index=YEARS, columns=AGES).fillna(0)
    return (c / pa * 1e5).mul(W, axis=1).sum(axis=1)

# ---------------- série anual Brasil ----------------
ann = pd.DataFrame({"cases": by_year(cn), "pop": pop_br, "notif": by_year(notif), "pulm_new_all": by_year(allp[allp.TRATAMENTO == "1"])}, index=YEARS)
ann["rate"] = ann.cases / ann["pop"] * 1e5; ann["asr"] = asr(cn, pop_age_br)
ann["rate_pulm_new_all"] = ann.pulm_new_all / ann["pop"] * 1e5
ann["pct_lab"] = by_year(allp[(allp.TRATAMENTO == "1") & allp.lab]) / ann.pulm_new_all * 100
ann["pct_retreat"] = (1 - ann.cases / ann.notif) * 100
for sx in ["Male", "Female"]:
    ann["rate_" + sx] = by_year(cn[cn.sexo == sx]) / pop_sex_br[sx] * 1e5
ann["ppl"] = by_year(cn[cn.ppl == 1]); ann["ppl_pct"] = ann.ppl / ann.cases * 100
ann["rate_nonppl"] = (ann.cases - ann.ppl) / ann["pop"] * 1e5
ann["pct_trm"] = cn.groupby("ano").TEST_MOLEC.apply(lambda x: x.isin(["1", "2", "3", "4"]).mean() * 100).reindex(YEARS)
ann["pct_only_trm"] = cn.groupby("ano").apply(lambda g: (g.TEST_MOLEC.isin(["1", "2"]) & (g.BACILOSC_E != "1") & (g.CULTURA_ES != "1")).mean() * 100).reindex(YEARS)
R["ann"] = ann

# taxas por UF e região
uf_rate = pd.DataFrame({u: by_year(cn[cn.uf == u]) / pop_uf[u] * 1e5 for u in UF.values()})
uf_asr = pd.DataFrame({u: asr(cn[cn.uf == u], pop_age_uf[u]) for u in UF.values()})
uf_ppl = pd.DataFrame({u: by_year(cn[(cn.uf == u) & (cn.ppl == 1)]) / by_year(cn[cn.uf == u]) * 100 for u in UF.values()})
uf_nonppl = pd.DataFrame({u: by_year(cn[(cn.uf == u) & (cn.ppl != 1)]) / pop_uf[u] * 1e5 for u in UF.values()})
reg_rate = pd.DataFrame({r: by_year(cn[cn.regiao == r]) / pop_reg[r] * 1e5 for r in REG.values()})
reg_ppl = pd.DataFrame({r: by_year(cn[(cn.regiao == r) & (cn.ppl == 1)]) / by_year(cn[cn.regiao == r]) * 100 for r in REG.values()})
R.update(uf_rate=uf_rate, uf_asr=uf_asr, uf_ppl=uf_ppl, uf_nonppl=uf_nonppl, reg_rate=reg_rate, reg_ppl=reg_ppl)

def apc(years, rates):
    x = sm.add_constant(np.array(years, float) - years[0]); y = np.log10(rates)
    res = sm.GLSAR(y, x, rho=1).iterative_fit(maxiter=50)
    b, se = res.params[1], res.bse[1]; t = stats.t.ppf(0.975, res.df_resid); f = lambda v: (10**v - 1) * 100
    lo, hi = f(b - t * se), f(b + t * se)
    return dict(APC=f(b), lo=lo, hi=hi, p=2 * (1 - stats.t.cdf(abs(b / se), res.df_resid)),
                trend="Increasing" if lo > 0 else ("Decreasing" if hi < 0 else "Stationary"))
ap = []
for nm, ser in [("Brazil, crude", ann.rate), ("Brazil, age-standardised", ann.asr), ("Brazil, excluding PDL", ann.rate_nonppl),
                ("Brazil, male", ann.rate_Male), ("Brazil, female", ann.rate_Female)] + [(r, reg_rate[r]) for r in REG.values()]:
    for lab, ys in [("2015-2019", range(2015, 2020)), ("2021-2025", range(2021, 2026)), ("2015-2025", range(2015, 2026))]:
        ys = list(ys); ap.append(dict(series=nm, period=lab, **apc(ys, ser.loc[ys].values)))
R["apc"] = pd.DataFrame(ap)

# ---------------- série mensal ----------------
months = pd.period_range("2015-01", "2025-12", freq="M")
t_p = list(months).index(pd.Period("2020-03", "M")); t_q = list(months).index(pd.Period("2022-01", "M"))
def monthly_frame(sub, popser):
    m = sub.groupby("mes").size().reindex(months, fill_value=0).rename("y").to_frame()
    mid = np.array([pd.Period(f"{y}-07", "M").ordinal for y in YEARS]); idx = np.array([p.ordinal for p in months])
    pv = popser.values.astype(float); s0 = (pv[1] - pv[0]) / 12; s1 = (pv[-1] - pv[-2]) / 12
    pop_m = np.interp(idx, mid, pv); pop_m = np.where(idx < mid[0], pv[0] + (idx - mid[0]) * s0, pop_m)
    pop_m = np.where(idx > mid[-1], pv[-1] + (idx - mid[-1]) * s1, pop_m)
    m["pop"] = pop_m; m["t"] = np.arange(len(m)); mo = np.array([p.month for p in months])
    m["dias"] = [p.days_in_month for p in months]
    for k in (1, 2):
        m[f"s{k}"] = np.sin(2 * np.pi * k * mo / 12); m[f"c{k}"] = np.cos(2 * np.pi * k * mo / 12)
    m["pand"] = ((m.t >= t_p) & (m.t < t_q)).astype(int); m["post"] = (m.t >= t_q).astype(int)
    m["t_pand"] = np.where(m.t >= t_p, m.t - t_p, 0); m["t_post"] = np.where(m.t >= t_q, m.t - t_q, 0)
    m["off"] = np.log(m["pop"]) + np.log(m.dias)
    return m
FOUR = "s1 + c1 + s2 + c2"; XC = ["Intercept", "t", "s1", "c1", "s2", "c2"]
def counterfactual(m, nsim=10000):
    pre = m[m.t < t_p]
    mp = smf.glm(f"y ~ t + {FOUR}", data=pre, family=sm.families.Poisson(), offset=pre.off).fit(scale="X2")
    X = sm.add_constant(m[["t", "s1", "c1", "s2", "c2"]]).values
    b = mp.params[XC].values; V = mp.cov_params().loc[XC, XC].values; ph = max(mp.scale, 1.0001)
    B = rng.multivariate_normal(b, V, nsim); MU = np.exp(B @ X.T + m.off.values)
    YS = rng.poisson(rng.gamma(MU / (ph - 1), ph - 1))
    e = np.exp(X @ b + m.off.values)
    out = {}
    for lab, sel in [("Pandemic", (m.t >= t_p) & (m.t < t_q)), ("Post", m.t >= t_q), ("Cumulative", m.t >= t_p)]:
        ii = np.where(sel.values)[0]; o = m.y.values[ii].sum(); ee = e[ii].sum(); s_ = YS[:, ii].sum(axis=1)
        out[lab] = dict(obs=o, exp=ee, diff=o - ee, diff_lo=o - np.percentile(s_, 97.5), diff_hi=o - np.percentile(s_, 2.5),
                        OE=o / ee, OE_lo=o / np.percentile(s_, 97.5), OE_hi=o / np.percentile(s_, 2.5))
    return out, e, np.percentile(YS, 2.5, axis=0), np.percentile(YS, 97.5, axis=0)

mb = monthly_frame(cn, pop_br)
its = smf.glm(f"y ~ t + pand + t_pand + post + t_post + {FOUR}", data=mb, family=sm.families.Poisson(), offset=mb.off).fit(
    scale="X2", cov_type="HAC", cov_kwds={"maxlags": 12})
R["its_params"] = pd.DataFrame({"b": its.params, "se": its.bse, "p": its.pvalues})
cov = its.cov_params()
def seg(names):
    b = sum(its.params[n] for n in names); se = np.sqrt(sum(cov.loc[a, c] for a in names for c in names))
    f = lambda x: (np.exp(12 * x) - 1) * 100; return dict(pct=f(b), lo=f(b - 1.96 * se), hi=f(b + 1.96 * se))
R["its_irr"] = {k: dict(IRR=np.exp(its.params[k]), lo=np.exp(its.params[k] - 1.96 * its.bse[k]), hi=np.exp(its.params[k] + 1.96 * its.bse[k]), p=its.pvalues[k]) for k in ["pand", "t_pand", "post", "t_post"]}
R["its_trends"] = {"Pre": seg(["t"]), "Pandemic": seg(["t", "t_pand"]), "Post": seg(["t", "t_pand", "t_post"])}
cf, e, lo, hi = counterfactual(mb); mb["exp"], mb["pi_lo"], mb["pi_hi"] = e, lo, hi; R["cf_br"] = cf
for nm, sub in [("nonppl", cn[cn.ppl != 1]), ("ppl", cn[cn.ppl == 1])]:
    mm = monthly_frame(sub, pop_br); o, e, lo, hi = counterfactual(mm)
    R[f"cf_{nm}"] = o; mb[f"y_{nm}"], mb[f"exp_{nm}"], mb[f"lo_{nm}"], mb[f"hi_{nm}"] = mm.y.values, e, lo, hi
R["monthly"] = mb
# por UF
cfu = {}
for u in UF.values():
    sub = cn[cn.uf == u]
    o1, *_ = counterfactual(monthly_frame(sub, pop_uf[u]), 5000)
    o2, *_ = counterfactual(monthly_frame(sub[sub.ppl != 1], pop_uf[u]), 5000)
    cfu[u] = dict(all=o1, nonppl=o2)
R["cf_uf"] = cfu
cfr = {}
for r in REG.values():
    sub = cn[cn.regiao == r]
    o1, *_ = counterfactual(monthly_frame(sub, pop_reg[r]), 5000); o2, *_ = counterfactual(monthly_frame(sub[sub.ppl != 1], pop_reg[r]), 5000)
    cfr[r] = dict(all=o1, nonppl=o2)
R["cf_reg"] = cfr

# ---------------- RT pós x pré ----------------
PER = {"Pre": range(2015, 2020), "Pandemic": range(2020, 2022), "Post": range(2022, 2026)}
def rt(a, pa, b, pb):
    r = (a / pa) / (b / pb); se = np.sqrt(1 / a + 1 / b); return r, r * np.exp(-1.96 * se), r * np.exp(1.96 * se)
rows = []
def addrow(g, c, cnt, popy):
    o = dict(group=g, cat=c)
    for p, ys in PER.items():
        ys = list(ys); o[f"n_{p}"] = cnt.loc[ys].sum(); o[f"rate_{p}"] = cnt.loc[ys].sum() / popy.loc[ys].sum() * 1e5
    for p in ["Pandemic", "Post"]:
        r_, l_, h_ = rt(o[f"n_{p}"], popy.loc[list(PER[p])].sum(), o["n_Pre"], popy.loc[list(PER["Pre"])].sum())
        o[f"RT_{p}"], o[f"RT_{p}_lo"], o[f"RT_{p}_hi"] = r_, l_, h_
    rows.append(o)
addrow("Total", "Brazil", by_year(cn), pop_br)
addrow("Total", "Excluding PDL", by_year(cn[cn.ppl != 1]), pop_br)
for sx in ["Male", "Female"]: addrow("Sex", sx, by_year(cn[cn.sexo == sx]), pop_sex_br[sx])
for a in AGES: addrow("Age group (years)", a, by_year(cn[cn.faixa == a]), pop_age_br[a])
for rc, pp in race_br.items(): addrow("Race/skin colour*", rc, by_year(cn[cn.raca == rc]), pd.Series(pp, index=YEARS))
for r in REG.values(): addrow("Region", r, by_year(cn[cn.regiao == r]), pop_reg[r])
R["rt"] = pd.DataFrame(rows)
rows = []
for u in UF.values():
    addrow("UF", u, by_year(cn[cn.uf == u]), pop_uf[u])
    addrow("UF_nonppl", u, by_year(cn[(cn.uf == u) & (cn.ppl != 1)]), pop_uf[u])
R["rt_uf"] = pd.DataFrame(rows[-54:]) if False else pd.DataFrame(rows)
rr = []
kb = (cn.raca == "White").sum()
for rc, pp in race_br.items():
    k = (cn.raca == rc).sum(); r_, l_, h_ = rt(k, pp, kb, race_br["White"])
    rr.append(dict(race=rc, n=k, rate=k / (pp * 11) * 1e5, RT=r_, lo=l_, hi=h_))
R["race"] = pd.DataFrame(rr).set_index("race"); R["race_ign"] = cn.raca.isna().mean() * 100

# ---------------- perfil por período ----------------
def dist(col, cats):
    t = pd.crosstab(cn[col].fillna("Missing"), cn.per)[["Pre", "Pandemic", "Post"]]
    t = t.reindex([c for c in cats if c in t.index] + [c for c in t.index if c not in cats])
    tk = t.drop(index=[i for i in t.index if i == "Missing"])
    chi2, p, dof, _ = stats.chi2_contingency(tk); V = np.sqrt(chi2 / (tk.values.sum() * (min(tk.shape) - 1)))
    return t, dict(p=p, V=V)
cn["ppl_s"] = cn.ppl.map({1: "Yes", 0: "No"}); cn["rua_s"] = cn.rua.map({1: "Yes", 0: "No"})
R["prof"] = {nm: dist(col, cats) for nm, col, cats in [
    ("Sex", "sexo", ["Male", "Female"]), ("Age group (years)", "faixa", AGES),
    ("Race/skin colour", "raca", ["White", "Black", "Brown (pardo)", "Asian", "Indigenous"]),
    ("Region", "regiao", list(REG.values())),
    ("Deprived of liberty", "ppl_s", ["Yes", "No"]), ("Homeless", "rua_s", ["Yes", "No"])]}
R["comorb_per"] = cn.groupby("per")[["alcool", "drogas", "tabaco", "diabetes"]].mean() * 100
R["ppl_male"] = (cn[cn.ppl == 1].sexo == "Male").mean() * 100
R["ppl_15_34"] = cn[cn.ppl == 1].faixa.isin(["15-24", "25-34"]).mean() * 100
R["ppl_share_m1534"] = cn[(cn.sexo == "Male") & cn.faixa.isin(["15-24", "25-34"])].ppl.mean() * 100

# ---------------- desfechos ----------------
def wilson(k, n):
    z = 1.96; ph = k / n; den = 1 + z**2 / n; c = (ph + z**2 / (2 * n)) / den
    h = z * np.sqrt(ph * (1 - ph) / n + z**2 / (4 * n**2)) / den; return ph * 100, (c - h) * 100, (c + h) * 100
cn["closed"] = cn.se.between(1, 10); cn["aband"] = cn.se.isin([2, 10])
oy = cn.groupby("ano").apply(lambda g: pd.Series({"n": len(g), "unclosed": (~g.closed).mean() * 100,
    "cure": (g.se == 1).sum() / g.closed.sum() * 100, "ltfu": g.aband.sum() / g.closed.sum() * 100,
    "death_tb": (g.se == 3).sum() / g.closed.sum() * 100, "death_other": (g.se == 4).sum() / g.closed.sum() * 100,
    "transfer": (g.se == 5).sum() / g.closed.sum() * 100}))
R["out_year"] = oy
R["cure_apc"] = apc(list(range(2015, 2025)), oy.loc[2015:2024, "cure"].values)
R["ltfu_apc"] = apc(list(range(2015, 2025)), oy.loc[2015:2024, "ltfu"].values)
coh = cn[(cn.ano <= 2024) & cn.closed].copy(); coh["perc"] = np.where(coh.ano <= 2019, "Pre", np.where(coh.ano <= 2021, "Pandemic", "Post"))
R["out_per"] = {p: dict(n=len(g), cure=wilson((g.se == 1).sum(), len(g)), ltfu=wilson(g.aband.sum(), len(g)),
                        death=wilson((g.se == 3).sum(), len(g))) for p, g in coh.groupby("perc")}
R["out_uf"] = coh.groupby(["uf", "perc"]).apply(lambda g: pd.Series({"cure": (g.se == 1).mean() * 100, "ltfu": g.aband.mean() * 100, "n": len(g)})).unstack()
# sub-grupos
R["out_sub"] = {nm: dict(n=int(msk.sum()), cure=wilson((coh[msk].se == 1).sum(), msk.sum()), ltfu=wilson(coh[msk].aband.sum(), msk.sum()))
                for nm, msk in [("PDL", coh.ppl == 1), ("Homeless", coh.rua == 1), ("HIV positive", coh.hiv == "Positive"),
                                ("Illicit drugs", coh.drogas == 1), ("Others", (coh.ppl != 1) & (coh.rua != 1) & (coh.hiv != "Positive"))]}
# modelo
rg = coh[coh.se != 5].copy(); rg["y"] = rg.aband.astype(int)
rg["raca_r"] = rg.raca.fillna("Missing").replace({"Asian": "Asian/Indigenous", "Indigenous": "Asian/Indigenous"})
rg["faixa_r"] = rg.faixa.fillna("25-34").replace({"0-14": "0-24", "15-24": "0-24"})
rg["hiv_r"] = rg.hiv.fillna("Not done/missing").replace({"Not done": "Not done/missing", "Pending": "Not done/missing"})
for v in ["ppl", "rua", "alcool", "drogas", "tabaco", "diabetes"]: rg[v + "_r"] = rg[v].fillna(0).astype(int)
rg["tdo_r"] = rg.TRAT_SUPER.map({"1": "Yes", "2": "No"}).fillna("Missing")
rg = rg.dropna(subset=["sexo"])
form = ("y ~ C(sexo, Treatment('Female')) + C(faixa_r, Treatment('35-44')) + C(raca_r, Treatment('White')) + ppl_r + rua_r + "
        "C(hiv_r, Treatment('Negative')) + alcool_r + drogas_r + tabaco_r + diabetes_r + C(tdo_r, Treatment('Yes')) + "
        "C(perc, Treatment('Pre')) + C(uf)")
mod = smf.glm(form, data=rg, family=sm.families.Poisson()).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(rg.uf)[0]})
res = pd.DataFrame({"PR": np.exp(mod.params), "lo": np.exp(mod.conf_int()[0]), "hi": np.exp(mod.conf_int()[1]), "p": mod.pvalues})
R["ltfu_model"] = res[~res.index.str.startswith("C(uf)")]; R["ltfu_n"] = len(rg); R["ltfu_prev"] = rg.y.mean() * 100

# ---------------- laboratório e HIV ----------------
R["lab_per"] = {p: dict(trm=wilson(g.TEST_MOLEC.isin(["1", "2", "3", "4"]).sum(), len(g)),
                        cult=wilson(g.CULTURA_ES.isin(["1", "2", "3"]).sum(), len(g)),
                        dst=wilson(g[g.CULTURA_ES == "1"].TEST_SENSI.isin(["1", "2", "3", "4", "5"]).sum(), (g.CULTURA_ES == "1").sum()),
                        rr=wilson(g.TEST_SENSI.isin(["2", "3"]).sum(), g.TEST_SENSI.isin(["1", "2", "3", "4", "5"]).sum()),
                        hr=wilson((g.TEST_SENSI == "1").sum(), g.TEST_SENSI.isin(["1", "2", "3", "4", "5"]).sum()),
                        mdr=wilson((g.TEST_SENSI == "3").sum(), g.TEST_SENSI.isin(["1", "2", "3", "4", "5"]).sum()),
                        rr_xpert=wilson((g.TEST_MOLEC == "2").sum(), g.TEST_MOLEC.isin(["1", "2"]).sum()),
                        hiv_test=wilson(g.hiv.isin(["Positive", "Negative"]).sum(), len(g)),
                        hiv_pos=wilson((g.hiv == "Positive").sum(), g.hiv.isin(["Positive", "Negative"]).sum()))
                for p, g in cn.groupby("per")}

# ---------------- municípios e regiões imediatas ----------------
MP = pd.read_csv("/home/claude/work/br/ext/mun_pop_rgi.csv", dtype={"cod7": str}); MP["mun"] = MP.cod7.str[:6]
mun = MP.set_index("mun")
for nm, sub in [("cases", cn), ("cases_pre", cn[cn.ano <= 2019]), ("cases_post", cn[cn.ano >= 2022]), ("cases_ppl", cn[cn.ppl == 1]),
                ("cases_pre_np", cn[(cn.ano <= 2019) & (cn.ppl != 1)]), ("cases_post_np", cn[(cn.ano >= 2022) & (cn.ppl != 1)])]:
    mun[nm] = sub.groupby("mun").size().reindex(mun.index).fillna(0)
mun["cases_np"] = mun.cases - mun.cases_ppl
mun["pct_ppl"] = np.where(mun.cases > 0, mun.cases_ppl / mun.cases * 100, np.nan)
R["mun"] = mun; R["mun_unmatched"] = int((~cn.mun.isin(mun.index)).sum())
rgi = mun.groupby("rgi")[["pop2022", "cases", "cases_np", "cases_pre_np", "cases_post_np", "cases_ppl"]].sum()
rgi["rate_np"] = rgi.cases_np / (rgi.pop2022 * 11) * 1e5
r_, l_, h_ = rt(rgi.cases_post_np.clip(lower=0.5), rgi.pop2022 * 4, rgi.cases_pre_np.clip(lower=0.5), rgi.pop2022 * 5)
rgi["RT_np"], rgi["RT_np_lo"], rgi["RT_np_hi"] = r_, l_, h_
names = pd.read_csv("/home/claude/work/br/ext/rgi.csv").set_index("rgi"); rgi = rgi.join(names)
R["rgi"] = rgi
R["cn_n"] = len(cn); R["notif_n"] = len(notif)
R["pop_br"] = pop_br; R["pop_uf"] = pop_uf
pickle.dump(R, open("/home/claude/work/br/res_br.pkl", "wb"))
cn.to_parquet("/home/claude/work/br/cn_br.parquet")
print("ok", len(cn))
