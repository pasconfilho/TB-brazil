"""Nova Figura 1 (Brasil): razão observado/esperado trimestral para casos fora das prisões, casos na PPL e óbitos por TB,
com painéis menores de contagens observadas e esperadas (com e sem tendência)."""
import pickle, warnings, numpy as np, pandas as pd, statsmodels.api as sm, statsmodels.formula.api as smf
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
warnings.filterwarnings("ignore"); rng = np.random.default_rng(2026)
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID, SHADE, SURF = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#efeeea", "#ffffff"
mpl.rcParams.update({"font.family": "Liberation Sans", "font.size": 7.5, "axes.titlesize": 8.5, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID,
    "grid.linewidth": 0.5, "legend.frameon": False, "axes.axisbelow": True, "xtick.major.size": 2.5, "ytick.major.size": 2.5})
R = pickle.load(open("/home/claude/work/br/res_br.pkl", "rb")); mb = R["monthly"].copy(); mb["t"] = np.arange(len(mb))
deaths = pd.read_csv("/home/claude/work/br/sim/sim_tb_mensal_uf.csv").Total.values
t_p, t_q = 62, 84
def sims(y, trend=True, start=0, nsim=5000):
    m = mb.copy(); m["y"] = y; pre = m[(m.t >= start) & (m.t < t_p)]
    cols = (["t"] if trend else []) + ["s1", "c1", "s2", "c2"]; XC = ["Intercept"] + cols
    mp = smf.glm("y ~ " + " + ".join(cols), data=pre, family=sm.families.Poisson(), offset=pre.off).fit(scale="X2")
    b = mp.params[XC].values; V = mp.cov_params().loc[XC, XC].values; ph = max(mp.scale, 1.0001)
    X = sm.add_constant(m[cols], has_constant="add").values
    MU = np.exp(rng.multivariate_normal(b, V, nsim) @ X.T + m.off.values); YS = rng.poisson(rng.gamma(MU / (ph - 1), ph - 1))
    return np.exp(X @ b + m.off.values), YS
q = np.arange(132) // 3; nq = 44
def quarterly(y, e, YS):
    oq = np.bincount(q, y); eq = np.bincount(q, e); Yq = np.stack([np.bincount(q, r) for r in YS])
    return oq, eq, oq / eq, oq / np.percentile(Yq, 97.5, axis=0), oq / np.percentile(Yq, 2.5, axis=0)
S = {}
for key, y in [("np", mb.y_nonppl.values), ("ppl", mb.y_ppl.values), ("deaths", deaths)]:
    e, YS = sims(y); ef, _ = sims(y, trend=False, start=26, nsim=10)
    S[key] = dict(y=y, q=quarterly(y, e, YS), ef=np.bincount(q, ef))
xq = 2015 + (np.arange(nq) * 3 + 1.5) / 12  # meio do trimestre
SER = [("np", "Cases outside prisons", BLUE, "o"), ("ppl", "Cases among PDL", ORANGE, "s"), ("deaths", "Tuberculosis deaths", AQUA, "^")]

fig = plt.figure(figsize=(7.2, 5.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1.55, 1], hspace=0.5, wspace=0.32)
ax = fig.add_subplot(gs[0, :])
ax.axvspan(2020 + 2 / 12, 2022, color=SHADE, lw=0, zorder=0)
ax.text(2021.08, 1.62, "Pandemic", ha="center", va="top", fontsize=7, color=INK2)
ax.axhline(1, color=INK2, lw=0.8, zorder=1)
for key, lab, c, mk in SER:
    oq, eq, r, lo, hi = S[key]["q"]
    ax.fill_between(xq, lo, hi, color=c, alpha=0.12, lw=0, zorder=2)
    ax.plot(xq, r, color=c, lw=1.6, zorder=3)
    ax.plot(xq, r, mk, color=c, ms=3.6, mec=SURF, mew=0.6, zorder=4)
ax.set_yscale("log"); ax.set_ylim(0.4, 1.7)
yt = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6]
ax.yaxis.set_major_locator(FixedLocator(yt)); ax.yaxis.set_major_formatter(FixedFormatter([str(v) for v in yt])); ax.yaxis.set_minor_locator(NullLocator())
ax.set_xlim(2015, 2026.75); ax.set_xticks(range(2015, 2026)); ax.grid(axis="y", color=GRID, lw=0.5)
ax.set_ylabel("Observed / expected (log scale)")
ax.set_title("A. Observed relative to the pre-pandemic trajectory, by quarter")
# rótulos diretos com O/E 2022-2025
oe_post = {"np": R["cf_nonppl"]["Post"], "ppl": R["cf_ppl"]["Post"]}
import json
sim = json.load(open("/home/claude/work/br/sim/sim_results.json"))["cf"]["Post"]
oe_post["deaths"] = dict(OE=sim["OE"], OE_lo=sim["lo"], OE_hi=sim["hi"])
lab_y = {"np": 1.08, "ppl": 0.5, "deaths": 1.42}
for key, lab, c, mk in SER:
    o = oe_post[key]
    ax.text(2026.05, lab_y[key], f"{lab}\n2022-25: {o['OE']:.2f} ({o['OE_lo']:.2f}-{o['OE_hi']:.2f})", color=c if False else INK, fontsize=6.8, va="center", ha="left")
ax.spines["bottom"].set_bounds(2015, 2026)
# painéis B: contagens trimestrais
titles = {"np": "B. Cases outside prisons", "ppl": "C. Cases among PDL", "deaths": "D. Tuberculosis deaths"}
for i, (key, lab, c, mk) in enumerate(SER):
    a = fig.add_subplot(gs[1, i]); oq, eq, *_ = S[key]["q"]
    a.axvspan(2020 + 2 / 12, 2022, color=SHADE, lw=0, zorder=0)
    a.plot(xq, eq, color=INK2, lw=1.0, ls=(0, (3, 1.5)), label="Expected, with trend")
    a.plot(xq[xq > 2020.1], S[key]["ef"][xq > 2020.1], color=INK2, lw=1.0, ls=(0, (1, 1.3)), label="Expected, flat 2017-19")
    a.plot(xq, oq, color=c, lw=1.6, label="Observed")
    a.set_title(titles[key], fontsize=8); a.set_xticks([2015, 2018, 2021, 2024]); a.set_xlim(2015, 2026)
    lo_, hi_ = min(oq.min(), eq.min()), max(oq.max(), eq.max(), S[key]["ef"].max())
    a.set_ylim(lo_ * 0.85, hi_ * 1.08)
    a.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, p: f"{v/1000:.1f}k" if v >= 1000 else f"{v:.0f}"))
    if i == 0: a.set_ylabel("Count per quarter")
    if i == 1: a.legend(loc="upper left", fontsize=6.2, handlelength=2.2, borderaxespad=0.2)
for ext, kw in [("png", {"dpi": 300}), ("pdf", {}), ("tif", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})]:
    fig.savefig(f"/home/claude/work/br/figv2/Figure1_v2.{ext}", bbox_inches="tight", facecolor="white", **kw)
print({k: (round(S[k]["q"][2][-1], 2)) for k in S})
