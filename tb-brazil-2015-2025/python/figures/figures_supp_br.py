"""Figuras suplementares em inglês: S1 fluxograma, S3 qualidade por UF."""
import pickle
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

R = pickle.load(open("/home/claude/work/br/res_br.pkl", "rb"))
OUT = "/home/claude/work/fig_br"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
mpl.rcParams.update({"font.family": "Liberation Sans", "font.size": 8})
def save(fig, name):
    for ext, kw in [("png", {}), ("pdf", {}), ("tif", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})]:
        fig.savefig(f"{OUT}/{name}.{ext}", bbox_inches="tight", facecolor="white", **({"dpi": 300} | kw))
    plt.close(fig)

# ---------- S1 flowchart ----------
fl = R["flow"]
steps = [(f"Tuberculosis records in national Sinan files, 2014-2025\nn = {fl[0][1]:,}", f"Diagnosis outside 2015-2025: {fl[0][1]-fl[1][1]:,}"),
         (f"Diagnosed 2015-2025\nn = {fl[1][1]:,}", f"Residence outside Brazil or invalid municipality: {fl[1][1]-fl[2][1]:,}"),
         (f"Resident in Brazil\nn = {fl[2][1]:,}", f"Extrapulmonary or mixed forms: {fl[2][1]-fl[3][1]:,}"),
         (f"Pulmonary tuberculosis only\nn = {fl[3][1]:,}", f"No laboratory confirmation*: {fl[3][1]-fl[4][1]:,}"),
         (f"Laboratory-confirmed pulmonary tuberculosis\nn = {fl[4][1]:,}", f"Retreatment, transfers and other entry types**: {fl[4][1]-fl[5][1]:,}"),
         (f"New laboratory-confirmed pulmonary cases\nn = {fl[5][1]:,}", None)]
fig, ax = plt.subplots(figsize=(6.6, 6.2)); ax.set_xlim(0, 10); ax.set_ylim(0, len(steps) * 1.6); ax.axis("off")
for i, (t, ex) in enumerate(steps):
    y = (len(steps) - i - 0.5) * 1.6
    ax.text(3.0, y, t, ha="center", va="center", fontsize=8.2, bbox=dict(boxstyle="square,pad=0.6", fc="white", ec=INK, lw=0.8),
            fontweight="bold" if i == len(steps) - 1 else "normal")
    if ex:
        ax.annotate("", xy=(3.0, y - 1.6 + 0.42), xytext=(3.0, y - 0.42), arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.8))
        ax.plot([3.0, 5.6], [y - 0.8, y - 0.8], color=INK, lw=0.8)
        ax.text(5.7, y - 0.8, "Excluded: " + ex, ha="left", va="center", fontsize=7.4, color=INK2,
                bbox=dict(boxstyle="square,pad=0.45", fc="#f7f7f5", ec=MUTED, lw=0.6), wrap=True)
save(fig, "FigureS1")

# ---------- S3 data quality heatmap ----------
Q = pd.read_csv("/home/claude/work/br/qualidade_uf.csv", index_col=0).rename(index={"Brasil": "Brazil"})
cols = [("casos", "New\ncases", None, "n"), ("sem_encerr_%", "No\noutcome*", False, "%"),
        ("raca_ign_%", "Race\nmissing", False, "%"), ("escol_ign_%", "Schooling\nmissing", False, "%"),
        ("ppl_ign_%", "PDL\nmissing", False, "%"), ("hiv_result_%", "HIV\nresult", True, "%"),
        ("trm_%", "Xpert\ndone", True, "%"), ("cultura_%", "Culture\ndone", True, "%"),
        ("transf_%", "Closed as\ntransfer*", False, "%"), ("ppl_%", "Cases\nin PDL", None, "%"),
        ("cura_%", "Cure*", True, "%"), ("abandono_%", "Loss to\nfollow-up*", False, "%")]
bad = LinearSegmentedColormap.from_list("b", ["#f7f7f5", "#fbdccf", "#f19b75", "#d24d1c"])
neu = LinearSegmentedColormap.from_list("n", ["#f7f7f5", "#e4e3df", "#b5b4ae"])
fig, ax = plt.subplots(figsize=(8.6, 7.6))
ax.set_xlim(0, len(cols)); ax.set_ylim(0, len(Q)); ax.invert_yaxis()
for j, (c, lab, hb, u) in enumerate(cols):
    v = Q[c].values[:-1]; lo, hi = np.nanmin(v), np.nanmax(v)
    for i, (uf, val) in enumerate(Q[c].items()):
        if hb is None: cm, x = neu, (val - lo) / (hi - lo + 1e-9)
        elif hb: cm, x = bad, (hi - val) / (hi - lo + 1e-9)
        else: cm, x = bad, (val - lo) / (hi - lo + 1e-9)
        x = min(max(x, 0), 1)
        ax.add_patch(plt.Rectangle((j, i), 1, 1, fc=cm(x * 0.95), ec="white", lw=1.2))
        txt = f"{val:,.0f}" if u == "n" else f"{val:.1f}"
        ax.text(j + 0.5, i + 0.5, txt, ha="center", va="center", fontsize=6.8,
                color="white" if x > 0.7 else INK, fontweight="bold" if uf == "Brazil" else "normal")
ax.set_xticks(np.arange(len(cols)) + 0.5); ax.set_xticklabels([c[1] for c in cols], fontsize=6.9); ax.xaxis.tick_top()
ax.set_yticks(np.arange(len(Q)) + 0.5); ax.set_yticklabels(Q.index, fontsize=7.5)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0); ax.axhline(len(Q) - 1, color=INK, lw=0.8)
save(fig, "FigureS3")
print("ok")
