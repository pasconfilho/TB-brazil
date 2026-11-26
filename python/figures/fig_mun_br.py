"""Figura S4 (funil + municípios mais afetados) e Tabela S6, Brasil."""
import json, numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

OUT = "/home/claude/work/fig_br"
BLUE, ORANGE, INK, INK2, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
mpl.rcParams.update({"font.family": "Liberation Sans", "font.size": 8, "axes.titlesize": 9, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "axes.axisbelow": True})
UF = {'11': 'RO', '12': 'AC', '13': 'AM', '14': 'RR', '15': 'PA', '16': 'AP', '17': 'TO', '21': 'MA', '22': 'PI', '23': 'CE',
      '24': 'RN', '25': 'PB', '26': 'PE', '27': 'AL', '28': 'SE', '29': 'BA', '31': 'MG', '32': 'ES', '33': 'RJ', '35': 'SP',
      '41': 'PR', '42': 'SC', '43': 'RS', '50': 'MS', '51': 'MT', '52': 'GO', '53': 'DF'}

g = pd.read_parquet("/home/claude/work/br/mun_spatial.parquet")
g["uf"] = g.cod7.str[:2].map(UF); g["lab"] = g.name + " (" + g.uf + ")"
g["r"] = g.cases_np / (g.pop2022 * 11) * 1e5
g["pre"] = g.cases_pre_np / (g.pop2022 * 5) * 1e5; g["post"] = g.cases_post_np / (g.pop2022 * 4) * 1e5
br = g.cases_np.sum() / (g.pop2022.sum() * 11)
pre_br = g.cases_pre_np.sum() / (g.pop2022.sum() * 5) * 1e5; post_br = g.cases_post_np.sum() / (g.pop2022.sum() * 4) * 1e5
big = g[g.pop2022 >= 100000].sort_values("r", ascending=False)
top = big.head(25).copy()
rr = lambda a, b: (a / 4) / (b / 5)
top["RR"] = rr(top.cases_post_np, top.cases_pre_np)
se = np.sqrt(1 / top.cases_post_np + 1 / top.cases_pre_np)
top["RR_lo"] = np.exp(np.log(top.RR) - 1.96 * se); top["RR_hi"] = np.exp(np.log(top.RR) + 1.96 * se)

fig, axs = plt.subplots(1, 2, figsize=(7.4, 5.4), gridspec_kw={"width_ratios": [1.15, 1], "wspace": 0.55})
# A: funnel
ax = axs[0]
P = np.logspace(np.log10(g.pop2022.min()), np.log10(g.pop2022.max()), 300)
# 99.8% limits: Poisson (dotted) and adjusted for overdispersion (solid; additive random effects, Spiegelhalter 2005)
FS = pd.read_csv("/home/claude/work/gh/tb-brazil-2015-2025/output/tables/funnel_summary.csv").set_index("measure").value
tau2 = (FS["tau (per 100,000)"] / 1e5) ** 2
for v, ls, lab_ in [(br / (P * 11), (0, (1, 1.5)), "Poisson 99.8% limit"), (br / (P * 11) + tau2, "-", "Overdispersion-adjusted 99.8% limit")]:
    ax.plot(P, (br + 3.09 * np.sqrt(v)) * 1e5, color=INK2, lw=0.7, ls=ls, label=lab_)
    ax.plot(P, np.clip(br - 3.09 * np.sqrt(v), 0, None) * 1e5, color=INK2, lw=0.7, ls=ls)
hi = g.pop2022 >= 100000
ax.scatter(g.pop2022[~hi], g.r[~hi], s=4, color=MUTED, alpha=0.35, lw=0, label="<100,000 inhabitants")
ax.scatter(g.pop2022[hi], g.r[hi], s=9, color=BLUE, alpha=0.75, lw=0, label="≥100,000 inhabitants")
ax.axhline(br * 1e5, color=ORANGE, lw=1)
lab = {"Guarujá": (5, 0, "left"), "Manaus": (5, 0, "left"), "São Vicente": (-5, 0, "right"), "Belém": (5, 0, "left"),
       "Rio de Janeiro": (0, 8, "center"), "São Paulo": (0, 8, "right"), "Duque de Caxias": (-5, 0, "right"), "Cubatão": (-5, 0, "right")}
for nm, (dx, dy, ha) in lab.items():
    r = g[(g.name == nm) & hi].iloc[0]
    ax.annotate(nm, (r.pop2022, r.r), xytext=(dx, dy), textcoords="offset points", fontsize=6.6, color=INK, va="center", ha=ha)
ax.set_xscale("log"); ax.set_ylim(0, 100); ax.grid(axis="y", color=GRID, lw=0.6)
ax.set_xlabel("Resident population (2022 Census, log scale)"); ax.set_ylabel("Mean annual rate per 100,000, excluding PDL")
ax.set_title("A. Municipal incidence, 2015-2025")
ax.legend(loc="upper left", fontsize=6.8, markerscale=1.6, handletextpad=0.2)
ax.text(1.4e7, br * 1e5 + 3, f"Brazil {br * 1e5:.1f}", ha="right", fontsize=6.6, color=ORANGE)

# B: dumbbell top 25
ax = axs[1]; y = np.arange(len(top))[::-1]
for yi, (_, r) in zip(y, top.iterrows()):
    c = BLUE if r.RR_lo > 1 else (ORANGE if r.RR_hi < 1 else MUTED)
    ax.plot([r.pre, r.post], [yi, yi], color=c, lw=1.3, alpha=0.8)
ax.scatter(top.pre, y, s=16, facecolor="white", edgecolor=INK2, lw=0.9, zorder=3, label="2015-2019")
ax.scatter(top.post, y, s=18, color=INK, zorder=3, label="2022-2025")
ax.axvline(pre_br, color=MUTED, lw=0.7, ls=(0, (2, 2))); ax.axvline(post_br, color=INK2, lw=0.7, ls=(0, (2, 2)))
ax.set_yticks(y); ax.set_yticklabels(top.lab, fontsize=6.8); ax.set_xlim(0, 100)
ax.grid(axis="x", color=GRID, lw=0.6); ax.tick_params(axis="y", length=0)
ax.set_xlabel("Rate per 100,000 per year, excluding PDL"); ax.set_title("B. Municipalities ≥100,000 with highest rates")
ax.legend(loc="lower right", fontsize=6.8, handletextpad=0.2)
for ext, kw in [("png", {"dpi": 300}), ("pdf", {}), ("tif", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})]:
    fig.savefig(f"{OUT}/FigureS4.{ext}", bbox_inches="tight", facecolor="white", **kw)

# Tabela S6
f0 = lambda x: f"{x:,.0f}"; f1 = lambda x: f"{x:.1f}"
rows = []
for _, r in top.iterrows():
    rows.append(dict(kind="", cells=[r.lab, f0(r.pop2022), f0(r.cases), f1(r.pct_ppl), f1(r.r), f1(r.pre), f1(r.post),
                                     f"{r.RR:.2f} ({r.RR_lo:.2f}-{r.RR_hi:.2f})", "Yes" if r.lisa_np == "High-high" else "No"]))
RRb = rr(g.cases_post_np.sum(), g.cases_pre_np.sum()); seb = np.sqrt(1 / g.cases_post_np.sum() + 1 / g.cases_pre_np.sum())
rows.append(dict(kind="bold", cells=["Brazil", f0(g.pop2022.sum()), f0(g.cases.sum()), f1(g.cases_ppl.sum() / g.cases.sum() * 100), f1(br * 1e5),
                                     f1(pre_br), f1(post_br), f"{RRb:.2f} ({np.exp(np.log(RRb) - 1.96 * seb):.2f}-{np.exp(np.log(RRb) + 1.96 * seb):.2f})", ""]))
T = json.load(open("/home/claude/work/br/tables_br.json"))
T["s6"] = dict(header=["Municipality (state)", "Population\n2022", "New cases", "PDL %", "Rate excl.\nPDL", "Rate\n2015-19", "Rate\n2022-25",
                       "RR post/pre\n(95% CI)", "High-high\ncluster"], rows=rows)
json.dump(T, open("/home/claude/work/br/tables_br.json", "w"), ensure_ascii=False, indent=1)
for r in rows: print(" | ".join(r["cells"]))
print("n big", len(big), "above 99.8% limit:", ((big.cases_np > br * big.pop2022 * 11 + 3.09 * np.sqrt(br * big.pop2022 * 11))).sum())
