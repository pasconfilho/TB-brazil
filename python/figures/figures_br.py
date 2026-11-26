import pickle, os
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects
from matplotlib.colors import BoundaryNorm, ListedColormap

R = pickle.load(open("/home/claude/work/br/res_br.pkl", "rb"))
OUT = "/home/claude/work/fig_br"; os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, AQUA, YEL = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, MUTED, GRID, SHADE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#efeeea"
mpl.rcParams.update({"font.family": "Liberation Sans", "font.size": 8.5, "axes.titlesize": 9.5, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID,
    "grid.linewidth": 0.6, "legend.frameon": False, "lines.linewidth": 2, "savefig.dpi": 300, "axes.axisbelow": True})
def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(f"{OUT}/{name}.tif", bbox_inches="tight", facecolor="white", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight", facecolor="white"); plt.close(fig)

m = R["monthly"]; x = np.arange(len(m)); tp, tq = 62, 84
ticks = [i for i in range(len(m)) if i % 12 == 0]; tl = [str(2015 + i // 12) for i in ticks]
def band(ax): ax.axvspan(tp - 0.5, tq - 0.5, color=SHADE, lw=0, zorder=0)

# Figure 1
fig, axs = plt.subplots(2, 1, figsize=(7.2, 5.6), gridspec_kw={"hspace": 0.45})
ax = axs[0]; band(ax); fut = x >= tp
ax.fill_between(x[fut], m.pi_lo[fut], m.pi_hi[fut], color=BLUE, alpha=0.14, lw=0, label="95% prediction interval")
ax.plot(x, m["exp"], color=BLUE, lw=1.4, ls=(0, (4, 2)), label="Expected (Jan 2015-Feb 2020 trend and seasonality)")
ax.plot(x, m.y, color=INK, lw=1.1, label="Observed")
ax.set_xticks(ticks); ax.set_xticklabels(tl); ax.set_xlim(-1, len(m)); ax.set_ylim(2500, 6500)
ax.set_ylabel("New cases per month"); ax.set_title("A. New laboratory-confirmed pulmonary TB cases, Brazil")
ax.legend(loc="upper left", fontsize=7, ncol=2, handlelength=1.8)
cf = R["cf_br"]
ax.text(73, 3000, f"O/E {cf['Pandemic']['OE']:.2f}", ha="center", fontsize=7.2, color=INK2)
ax.text(108, 3000, f"O/E {cf['Post']['OE']:.2f}", ha="center", fontsize=7.2, color=INK2)
ax.text((tp + tq) / 2, 2650, "Pandemic", ha="center", fontsize=7.5, color=INK2)
ax = axs[1]; band(ax)
for col, c, lab in [("nonppl", BLUE, "Not deprived of liberty"), ("ppl", ORANGE, "Deprived of liberty")]:
    ax.fill_between(x[fut], m[f"lo_{col}"][fut], m[f"hi_{col}"][fut], color=c, alpha=0.13, lw=0)
    ax.plot(x, m[f"exp_{col}"], color=c, lw=1.3, ls=(0, (4, 2))); ax.plot(x, m[f"y_{col}"], color=c, lw=1.1, label=lab)
ax.set_xticks(ticks); ax.set_xticklabels(tl); ax.set_xlim(-1, len(m)); ax.set_ylim(0, 5600)
ax.set_ylabel("New cases per month"); ax.set_title("B. By deprivation of liberty")
ax.text(131, 4900, f"O/E 2022-2025: {R['cf_nonppl']['Post']['OE']:.2f}", ha="right", fontsize=7.2, color=INK2)
ax.text(131, 1150, f"O/E 2022-2025: {R['cf_ppl']['Post']['OE']:.2f}", ha="right", fontsize=7.2, color=INK2)
ax.legend(loc="upper left", fontsize=7.5, ncol=2)
save(fig, "Figure1")

# Figure 2: O/E por UF (não PPL), pandemia e pós
REGORD = {"North": ["RO", "AC", "AM", "RR", "PA", "AP", "TO"], "Northeast": ["MA", "PI", "CE", "RN", "PB", "PE", "AL", "SE", "BA"],
          "Southeast": ["MG", "ES", "RJ", "SP"], "South": ["PR", "SC", "RS"], "Centre-West": ["MS", "MT", "GO", "DF"]}
fig, axs = plt.subplots(1, 2, figsize=(7.2, 6.6), sharey=True, gridspec_kw={"wspace": 0.08})
rows = []
for reg, ufs in REGORD.items():
    rows.append(("reg", reg, R["cf_reg"][reg]["nonppl"]))
    for u in ufs: rows.append(("uf", u, R["cf_uf"][u]["nonppl"]))
rows.append(("br", "Brazil", R["cf_nonppl"]))
yy = np.arange(len(rows))[::-1]
for ax, per, ttl in [(axs[0], "Pandemic", "A. Pandemic (Mar 2020-Dec 2021)"), (axs[1], "Post", "B. Post-pandemic (2022-2025)")]:
    for y_, (kind, nm, o) in zip(yy, rows):
        v = o[per]; c = ORANGE if v["OE_lo"] > 1 else (BLUE if v["OE_hi"] < 1 else MUTED)
        ms = 6 if kind != "uf" else 4.5; mk = "D" if kind != "uf" else "o"
        ax.plot([v["OE_lo"], v["OE_hi"]], [y_, y_], color=c, lw=1.4 if kind == "uf" else 2)
        ax.plot(v["OE"], y_, mk, color=c if kind == "uf" else INK, ms=ms, mec="white", mew=0.8, zorder=3)
    ax.axvline(1, color=INK2, lw=0.9); ax.set_xscale("log"); ax.set_xticks([0.5, 0.75, 1, 1.5, 2]); ax.set_xticklabels(["0.5", "0.75", "1", "1.5", "2"])
    ax.minorticks_off(); ax.set_xlim(0.45, 2.3); ax.grid(axis="x", color=GRID); ax.grid(axis="y", visible=False)
    ax.set_title(ttl, fontsize=9); ax.set_xlabel("Observed/expected ratio (95% PI)")
axs[0].set_yticks(yy); axs[0].set_yticklabels([nm for _, nm, _ in rows], fontsize=7.4)
for lab_ in axs[0].get_yticklabels():
    if lab_.get_text() in REGORD or lab_.get_text() == "Brazil": lab_.set_fontweight("bold")
save(fig, "Figure2")

# Figure 3: forest RT
rt = R["rt"]
groups = [("Total", rt[rt.group == "Total"]), ("Sex", rt[rt.group == "Sex"]), ("Age group (years)", rt[rt.group == "Age group (years)"]),
          ("Race/skin colour*", rt[rt.group == "Race/skin colour*"]), ("Region", rt[rt.group == "Region"])]
fig, ax = plt.subplots(figsize=(7.2, 6.6)); yv = 0
for title, sub in groups:
    ax.text(-0.02, yv, title, transform=ax.get_yaxis_transform(), ha="right", va="center", fontweight="bold", fontsize=8.5); yv -= 1
    for r in sub.itertuples():
        c = ORANGE if r.RT_Post_lo > 1 else (BLUE if r.RT_Post_hi < 1 else MUTED)
        ax.plot([r.RT_Post_lo, r.RT_Post_hi], [yv, yv], color=c, lw=1.6); ax.plot(r.RT_Post, yv, "o", color=c, ms=5.5, mec="white", mew=1, zorder=3)
        lab = r.cat.replace("Excluding PDL", "Excluding people deprived of liberty").replace("65+", "≥65")
        ax.text(-0.02, yv, lab, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8, color=INK2)
        ax.text(1.02, yv, f"{r.rate_Pre:.1f}  →  {r.rate_Post:.1f}", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=7.6, color=INK2)
        ax.text(1.30, yv, f"{r.RT_Post:.2f} ({r.RT_Post_lo:.2f}–{r.RT_Post_hi:.2f})", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=7.6, color=INK2)
        yv -= 1
    yv -= 0.4
ax.axvline(1, color=INK2, lw=0.9); ax.set_xscale("log"); tk = [0.9, 1, 1.2, 1.4, 1.6, 1.8]
ax.set_xticks(tk); ax.set_xticklabels([str(t) for t in tk]); ax.minorticks_off(); ax.set_xlim(0.88, 1.9); ax.set_ylim(yv + 0.4, 0.8)
ax.set_yticks([]); ax.grid(axis="x", color=GRID); ax.grid(axis="y", visible=False); ax.spines["left"].set_visible(False)
ax.text(1.02, 0.9, "Rate pre → post\n(per 100,000/year)", transform=ax.get_yaxis_transform(), fontsize=7.6, fontweight="bold", va="bottom")
ax.text(1.30, 0.9, "RR (95% CI)", transform=ax.get_yaxis_transform(), fontsize=7.6, fontweight="bold", va="bottom")
ax.set_xlabel("Incidence rate ratio, post-pandemic (2022-2025) vs pre-pandemic (2015-2019)")
save(fig, "Figure3")

# Figure 4: maps
g = gpd.read_file("/home/claude/work/br/mun_spatial.gpkg"); g["geometry"] = g.geometry.make_valid().buffer(0)
uf_b = g.dissolve(by=g.cod7.str[:2]); uf_b["geometry"] = uf_b.buffer(0.03).buffer(-0.03)
rgi = R["rgi"]
g["rgi"] = g["rgi"].astype(int)
g = g.join(rgi[["RT_np"]], on="rgi"); grgi = g
fig, axs = plt.subplots(1, 3, figsize=(10.5, 4.2), gridspec_kw={"wspace": 0.02})
bounds = [0, 5, 10, 20, 30, 45, 60, 1e9]; ramp = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
cmap = ListedColormap(ramp); norm = BoundaryNorm(bounds, cmap.N)
ax = axs[0]; g.plot(column="seb_np", cmap=cmap, norm=norm, ax=ax, linewidth=0)
uf_b.boundary.plot(ax=ax, color="white", lw=0.4); ax.set_axis_off(); ax.set_title("A. Incidence excluding PDL\n(spatial empirical Bayes)", fontsize=9)
labs = ["<5", "5–9", "10–19", "20–29", "30–44", "45–59", "≥60"]
ax.legend(handles=[mpl.patches.Patch(fc=c, label=l) for c, l in zip(ramp, labs)], loc="lower left", fontsize=6.4, title="per 100,000/yr", title_fontsize=6.6)
CITY = {"1302603": ("Manaus", -6, 4, "right"), "1501402": ("Belém", 4, 4, "left"), "2304400": ("Fortaleza", 4, 3, "left"), "2611606": ("Recife", 5, 0, "left"),
        "2927408": ("Salvador", 5, -2, "left"), "3304557": ("Rio de Janeiro", 5, -4, "left"), "3550308": ("São Paulo", -5, -5, "right"), "4314902": ("Porto Alegre", 5, -3, "left")}
gp = g.to_crs(5880)
for cod, (nm, dx, dy, ha) in CITY.items():
    pt = gp[gp.cod7 == cod].geometry.representative_point().to_crs(4326).iloc[0]
    ax.plot(pt.x, pt.y, "o", ms=2.6, color=INK, mec="white", mew=0.5, zorder=5)
    ax.annotate(nm, (pt.x, pt.y), xytext=(dx, dy), textcoords="offset points", fontsize=6.2, color=INK, ha=ha, va="center", zorder=6,
                path_effects=[mpl.patheffects.withStroke(linewidth=1.8, foreground="white")])
ax = axs[1]; colors = {"High-high": "#c0392b", "Low-low": "#2a78d6", "High-low": "#f19b75", "Low-high": "#9ec5f4", "Not significant": "#efeeea"}
g.plot(color=g.lisa_np.map(colors), ax=ax, linewidth=0); uf_b.boundary.plot(ax=ax, color="#8a8984", lw=0.3)
ax.set_axis_off(); ax.set_title("B. Local clusters (LISA)\nexcluding PDL", fontsize=9)
ax.legend(handles=[mpl.patches.Patch(fc=c, label=l) for l, c in colors.items()], loc="lower left", fontsize=6.4)
ax = axs[2]
b2 = [0, 0.8, 0.95, 1.05, 1.2, 1.5, 100]; r2 = ["#2a78d6", "#9ec5f4", "#efeeea", "#fbdccf", "#f19b75", "#c0392b"]
grgi.plot(column="RT_np", cmap=ListedColormap(r2), norm=BoundaryNorm(b2, 6), ax=ax, linewidth=0)
uf_b.boundary.plot(ax=ax, color="#52514e", lw=0.35); ax.set_axis_off()
ax.set_title("C. Change post vs pre-pandemic excluding PDL\n(immediate geographic regions)", fontsize=9)
ax.legend(handles=[mpl.patches.Patch(fc=c, label=l) for c, l in zip(r2, ["<0.80", "0.80–0.94", "0.95–1.04", "1.05–1.19", "1.20–1.49", "≥1.50"])],
          loc="lower left", fontsize=6.4, title="Rate ratio", title_fontsize=6.6)
save(fig, "Figure4")

# Figure 5: outcomes
oy = R["out_year"].loc[2015:2024]; yrs = oy.index.values
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.9), gridspec_kw={"wspace": 0.85, "width_ratios": [1.15, 1.25]})
ax = axs[0]; ax.axvspan(2019.5, 2021.5, color=SHADE, lw=0, zorder=0)
for col, c, mk, lab in [("cure", BLUE, "o", "Cure"), ("ltfu", ORANGE, "s", "Loss to follow-up"), ("death_tb", AQUA, "^", "TB death")]:
    ax.plot(yrs, oy[col], color=c, marker=mk, ms=4, mec="white", mew=1); ax.text(2024.4, oy[col].iloc[-1] + (1.5 if col == "ltfu" else (-1.5 if col == "death_tb" else 0)), lab, va="center", fontsize=7.4)
ax.axhline(85, color=INK2, lw=0.8, ls=(0, (3, 2))); ax.axhline(5, color=INK2, lw=0.8, ls=(0, (3, 2)))
ax.text(2015, 86.5, "Cure target ≥85%", fontsize=6.8, color=INK2); ax.text(2015, 6.3, "Loss to follow-up target <5%", fontsize=6.8, color=INK2)
ax.set_xlim(2014.6, 2028.2); ax.set_ylim(0, 95); ax.set_xticks(yrs[::2]); ax.set_xticklabels([str(y) for y in yrs[::2]], fontsize=7.4)
ax.set_ylabel("% of closed cases"); ax.set_xlabel("Year of diagnosis"); ax.set_title("A. Treatment outcomes")
ax = axs[1]; rp = R["ltfu_model"]
lab = {"C(sexo, Treatment('Female'))[T.Male]": "Male sex", "C(faixa_r, Treatment('35-44'))[T.0-24]": "Age <25 years¹",
       "C(faixa_r, Treatment('35-44'))[T.65+]": "Age ≥65 years¹", "C(raca_r, Treatment('White'))[T.Black]": "Black²",
       "C(raca_r, Treatment('White'))[T.Brown (pardo)]": "Brown (pardo)²", "rua_r": "Homelessness", "drogas_r": "Illicit drug use",
       "alcool_r": "Alcohol use", "tabaco_r": "Smoking", "C(hiv_r, Treatment('Negative'))[T.Positive]": "HIV positive³",
       "C(hiv_r, Treatment('Negative'))[T.Not done/missing]": "HIV not tested³", "C(tdo_r, Treatment('Yes'))[T.No]": "No DOT planned⁴",
       "ppl_r": "Deprived of liberty", "C(perc, Treatment('Pre'))[T.Pandemic]": "Pandemic period⁵", "C(perc, Treatment('Pre'))[T.Post]": "Post-pandemic period⁵"}
items = [(v, rp.loc[k]) for k, v in lab.items() if k in rp.index]; yy = np.arange(len(items))[::-1]
for y_, (nm, r) in zip(yy, items):
    c = ORANGE if r.lo > 1 else (BLUE if r.hi < 1 else MUTED)
    ax.plot([r.lo, r.hi], [y_, y_], color=c, lw=1.6); ax.plot(r.PR, y_, "o", color=c, ms=5, mec="white", mew=1, zorder=3)
    ax.text(3.6, y_, f"{r.PR:.2f} ({r.lo:.2f}–{r.hi:.2f})", va="center", fontsize=7, color=INK2)
ax.set_yticks(yy); ax.set_yticklabels([i[0] for i in items], fontsize=7.6); ax.axvline(1, color=INK2, lw=0.9)
ax.set_xscale("log"); ax.set_xticks([0.25, 0.5, 1, 2]); ax.set_xticklabels(["0.25", "0.5", "1", "2"]); ax.minorticks_off(); ax.set_xlim(0.3, 3.4)
ax.grid(axis="x", color=GRID); ax.grid(axis="y", visible=False); ax.set_xlabel("Adjusted prevalence ratio (95% CI)"); ax.set_title("B. Factors associated with loss to follow-up")
save(fig, "Figure5")

# Supplementary: PDL share by region; lab
fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.9), gridspec_kw={"wspace": 0.5})
ax = axs[0]; rp_ = R["reg_ppl"]; cols = [BLUE, ORANGE, AQUA, YEL, "#6250d6"]
for (r, c) in zip(rp_.columns, cols): ax.plot(rp_.index, rp_[r], color=c, marker="o", ms=3.5, mec="white", mew=0.8, label=r, lw=1.6)
ax.set_ylabel("% of new cases in PDL"); ax.set_title("A. Share of cases deprived of liberty"); ax.legend(fontsize=6.8, ncol=2, loc="upper left")
ax.set_ylim(0, 32); ax.set_xticks(rp_.index[::2])
ax = axs[1]; ann = R["ann"]
ax.plot(ann.index, ann.pct_trm, color=BLUE, marker="o", ms=3.5, mec="white", mew=0.8); ax.plot(ann.index, ann.pct_only_trm, color=ORANGE, marker="s", ms=3.5, mec="white", mew=0.8)
ax.plot(ann.index, ann.pct_lab, color=AQUA, marker="^", ms=3.5, mec="white", mew=0.8)
ax.text(2025.3, ann.pct_trm.iloc[-1] - 3, "Xpert done", va="center", fontsize=7); ax.text(2025.3, ann.pct_only_trm.iloc[-1], "Confirmed\nonly by Xpert", va="center", fontsize=7)
ax.text(2025.3, ann.pct_lab.iloc[-1] + 4, "Lab-confirmed\n(all new pulm.)", va="center", fontsize=7)
ax.set_xlim(2014.6, 2029); ax.set_ylim(0, 100); ax.set_xticks(ann.index[::2]); ax.set_ylabel("%"); ax.set_title("B. Laboratory diagnosis")
save(fig, "FigureS2")
print(os.listdir(OUT))
