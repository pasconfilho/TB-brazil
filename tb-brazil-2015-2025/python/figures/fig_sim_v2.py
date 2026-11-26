import json, pandas as pd, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
BLUE, ORANGE, INK, INK2, MUTED, GRID, SHADE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#efeeea"
mpl.rcParams.update({"font.family": "Liberation Sans", "font.size": 7.5, "axes.titlesize": 8.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "grid.linewidth": 0.5})
A = pd.DataFrame(json.load(open("/home/claude/work/br/sim/sim_results.json"))["annual"]).set_index("ano"); yrs = A.index.values
M = pd.read_csv("/home/claude/work/gh/tb-brazil-2015-2025/output/tables/mortality_annual.csv").set_index("year")
A["asmr"] = M.loc[yrs, "asmr"].values
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw={"wspace": 0.3})
for ax, col, c, mk, ttl, yl, fmt in [(axs[0], "mort_rate", BLUE, "o", "A. Tuberculosis mortality", "Deaths per 100,000", "{:.1f}"),
                                     (axs[1], "deaths_per_100_cases", ORANGE, "s", "B. Deaths per 100 new notified cases", "Deaths per 100 cases", "{:.1f}")]:
    ax.axvspan(2019.5, 2021.5, color=SHADE, lw=0, zorder=0)
    ax.plot(yrs, A[col], color=c, lw=1.8, marker=mk, ms=4.5, mec="white", mew=0.8)
    for y in (2015, 2019, 2025):
        ax.annotate(fmt.format(A.loc[y, col]), (y, A.loc[y, col]), xytext=(0, 7 if (ax is axs[1] or y == 2025) else -13), textcoords="offset points", ha="center", fontsize=6.8, color=INK)
    lo, hi = A[col].min(), A[col].max(); ax.set_ylim(lo - (hi - lo) * 0.35, hi + (hi - lo) * 0.3)
    ax.set_title(ttl); ax.set_ylabel(yl); ax.set_xticks(range(2015, 2026, 2)); ax.set_xlim(2014.5, 2025.5)
ax = axs[0]; ax.plot(yrs, A["asmr"], color=INK2, lw=1.2, ls=(0, (3, 2)), marker="o", ms=3, mfc="white", mec=INK2, mew=0.8)
lo = min(A["asmr"].min(), A["mort_rate"].min()); hi = max(A["asmr"].max(), A["mort_rate"].max()); ax.set_ylim(lo - (hi - lo) * 0.35, hi + (hi - lo) * 0.3)
ax.annotate("Crude", (2017, A.loc[2017, "mort_rate"]), xytext=(0, -12), textcoords="offset points", ha="center", fontsize=6.6, color=BLUE)
ax.annotate("Age-standardised", (2017, A.loc[2017, "asmr"]), xytext=(4, 7), textcoords="offset points", ha="center", fontsize=6.6, color=INK2)
axs[0].text(2020.5, axs[0].get_ylim()[0] + 0.02, "Pandemic", ha="center", va="bottom", fontsize=6.6, color=INK2)
for ext, kw in [("png", {"dpi": 300}), ("pdf", {}), ("tif", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})]:
    fig.savefig(f"/home/claude/work/fig_br/FigureS5.{ext}", bbox_inches="tight", facecolor="white", **kw)
