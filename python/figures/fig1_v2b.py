"""Versão alternativa da Figura 1: observado x esperado (mensal) para cada série, com faixa de predição e diferença sombreada."""
import json, pickle, warnings, numpy as np, pandas as pd, statsmodels.api as sm, statsmodels.formula.api as smf
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore"); rng = np.random.default_rng(2026)
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID, SHADE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#efeeea"
mpl.rcParams.update({"font.family": "Liberation Sans", "font.size": 7.5, "axes.titlesize": 8.5, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID,
    "grid.linewidth": 0.5, "legend.frameon": False, "axes.axisbelow": True, "xtick.major.size": 2.5, "ytick.major.size": 2.5})
R = pickle.load(open("/home/claude/work/br/res_br.pkl", "rb")); mb = R["monthly"].copy(); mb["t"] = np.arange(len(mb))
deaths = pd.read_csv("/home/claude/work/br/sim/sim_tb_mensal_uf.csv").Total.values
sim = json.load(open("/home/claude/work/br/sim/sim_results.json"))["cf"]
t_p, t_q = 62, 84
def sims(y, trend=True, start=0, nsim=5000):
    m = mb.copy(); m["y"] = y; pre = m[(m.t >= start) & (m.t < t_p)]
    cols = (["t"] if trend else []) + ["s1", "c1", "s2", "c2"]; XC = ["Intercept"] + cols
    mp = smf.glm("y ~ " + " + ".join(cols), data=pre, family=sm.families.Poisson(), offset=pre.off).fit(scale="X2")
    b = mp.params[XC].values; V = mp.cov_params().loc[XC, XC].values; ph = max(mp.scale, 1.0001)
    X = sm.add_constant(m[cols], has_constant="add").values
    MU = np.exp(rng.multivariate_normal(b, V, nsim) @ X.T + m.off.values); YS = rng.poisson(rng.gamma(MU / (ph - 1), ph - 1))
    return np.exp(X @ b + m.off.values), np.percentile(YS, 2.5, axis=0), np.percentile(YS, 97.5, axis=0)
x = 2015 + (np.arange(132) + 0.5) / 12
oe = {"np": (R["cf_nonppl"]["Pandemic"], R["cf_nonppl"]["Post"]), "ppl": (R["cf_ppl"]["Pandemic"], R["cf_ppl"]["Post"]),
      "deaths": ({"OE": sim["Pandemic"]["OE"], "OE_lo": sim["Pandemic"]["lo"], "OE_hi": sim["Pandemic"]["hi"]},
                 {"OE": sim["Post"]["OE"], "OE_lo": sim["Post"]["lo"], "OE_hi": sim["Post"]["hi"]})}
SER = [("np", mb.y_nonppl.values, "A. New laboratory-confirmed pulmonary cases outside prisons", BLUE, "Cases per month"),
       ("ppl", mb.y_ppl.values, "B. New laboratory-confirmed pulmonary cases among people deprived of liberty", ORANGE, "Cases per month"),
       ("deaths", deaths, "C. Deaths with tuberculosis as underlying cause", AQUA, "Deaths per month")]
fig, axs = plt.subplots(3, 1, figsize=(7.2, 7.4), sharex=True, gridspec_kw={"hspace": 0.32})
for ax, (key, y, ttl, c, yl) in zip(axs, SER):
    e, lo, hi = sims(y); f = x >= 2020 + 2 / 12
    ax.axvspan(2020 + 2 / 12, 2022, color=SHADE, lw=0, zorder=0)
    ax.fill_between(x[f], lo[f], hi[f], color=INK2, alpha=0.10, lw=0, zorder=1, label="95% prediction interval")
    ax.fill_between(x[f], y[f], e[f], where=y[f] < e[f], color=c, alpha=0.22, lw=0, interpolate=True, zorder=2)
    ax.fill_between(x[f], y[f], e[f], where=y[f] >= e[f], color=c, alpha=0.45, lw=0, interpolate=True, zorder=2)
    ax.plot(x, e, color=INK2, lw=1.1, ls=(0, (4, 2)), zorder=3, label="Expected from Jan 2015-Feb 2020 trend")
    if key == "ppl":
        ef, _, _ = sims(y, trend=False, start=26, nsim=10)
        ax.plot(x[f], ef[f], color=INK2, lw=1.0, ls=(0, (1, 1.4)), zorder=3, label="Expected, flat 2017-2019 level")
        ax.text(2026.05, ef[-1], "Flat\nbaseline", fontsize=6.5, color=INK2, va="center")
        ax.text(2026.05, e[-1], "With\ntrend", fontsize=6.5, color=INK2, va="center")
    ax.plot(x, y, color=c, lw=1.3, zorder=4, label="Observed")
    ax.set_title(ttl); ax.set_ylabel(yl)
    top = max(y.max(), hi.max()); bot = min(y.min(), lo[f].min())
    ax.set_ylim(bot * 0.8, top * 1.28)
    p, q = oe[key]
    ax.text(2021.08, top * 1.16, f"Mar 2020-Dec 2021\nO/E {p['OE']:.2f} ({p['OE_lo']:.2f}-{p['OE_hi']:.2f})", ha="center", va="top", fontsize=6.8, color=INK)
    ax.text(2024.0, top * 1.16, f"2022-2025\nO/E {q['OE']:.2f} ({q['OE_lo']:.2f}-{q['OE_hi']:.2f})", ha="center", va="top", fontsize=6.8, color=INK)
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, p: f"{v:,.0f}"))
axs[0].text(2021.08, axs[0].get_ylim()[0] + 0.03 * np.diff(axs[0].get_ylim())[0], "Pandemic", ha="center", fontsize=6.8, color=INK2)
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
H = [Line2D([], [], color=INK, lw=1.3, label="Observed (series colour)"), Line2D([], [], color=INK2, lw=1.1, ls=(0, (4, 2)), label="Expected, Jan 2015-Feb 2020 trend"),
     Line2D([], [], color=INK2, lw=1.0, ls=(0, (1, 1.4)), label="Expected, flat 2017-2019 level (B)"), Patch(color=INK2, alpha=0.10, lw=0, label="95% prediction interval"),
     Patch(color=MUTED, alpha=0.30, lw=0, label="Fewer than expected"), Patch(color=MUTED, alpha=0.65, lw=0, label="More than expected")]
fig.legend(handles=H, loc="lower center", ncol=3, fontsize=6.6, bbox_to_anchor=(0.5, 0.0), handlelength=2.6, columnspacing=1.6)
axs[-1].set_xlim(2015, 2026); axs[-1].set_xticks(range(2015, 2027))
for ext, kw in [("png", {"dpi": 300}), ("pdf", {}), ("tif", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})]:
    fig.savefig(f"/home/claude/work/br/figv2/Figure1_v2b.{ext}", bbox_inches="tight", facecolor="white", **kw)
