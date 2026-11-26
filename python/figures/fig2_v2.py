"""Figura 2 nova: mapa em grade das UFs com O/E (excluindo PPL), pandemia e pós-pandemia."""
import pickle, numpy as np, matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
mpl.rcParams.update({"font.family": "Liberation Sans", "font.size": 7.5, "axes.titlesize": 8.5, "axes.titleweight": "bold", "axes.titlelocation": "left"})
R = pickle.load(open("/home/claude/work/br/res_br.pkl", "rb"))
POS = {"RR": (0, 2), "AP": (0, 4), "AM": (1, 1), "PA": (1, 3), "MA": (1, 4), "CE": (1, 5), "RN": (1, 6),
       "AC": (2, 0), "RO": (2, 1), "MT": (2, 2), "TO": (2, 3), "PI": (2, 4), "PE": (2, 5), "PB": (2, 6),
       "MS": (3, 2), "GO": (3, 3), "BA": (3, 4), "SE": (3, 5), "AL": (3, 6),
       "PR": (4, 2), "DF": (4, 3), "MG": (4, 4), "ES": (4, 5), "SC": (5, 2), "SP": (5, 3), "RJ": (5, 4), "RS": (6, 2)}
BINS = [0, 0.75, 0.85, 0.95, 1.05, 1.15, 1.30, 9]
COLS = ["#1c5cab", "#5598e7", "#9ec5f4", "#f0efec", "#f4b3ae", "#e86c66", "#b8302f"]
LAB = ["<0.75", "0.75-0.84", "0.85-0.94", "0.95-1.04", "1.05-1.14", "1.15-1.29", "≥1.30"]
def col(v): return COLS[np.searchsorted(BINS, round(v, 2) + 1e-9, side="right") - 1]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 4.4), gridspec_kw={"wspace": 0.06})
for ax, per, ttl in [(axs[0], "Pandemic", "A. Pandemic (March 2020-December 2021)"), (axs[1], "Post", "B. Post-pandemic (2022-2025)")]:
    for uf, (r, c) in POS.items():
        o = R["cf_uf"][uf]["nonppl"][per]; v = o["OE"]; sig = o["OE_lo"] > 1 or o["OE_hi"] < 1
        fc = col(v); dark = fc in ("#b8302f", "#1c5cab")
        ax.add_patch(FancyBboxPatch((c + 0.05, -r + 0.05), 0.9, 0.9, boxstyle="round,pad=0,rounding_size=0.08", fc=fc,
                                    ec=INK if sig else "white", lw=1.3 if sig else 0.8))
        tc = "white" if dark else INK
        ax.text(c + 0.5, -r + 0.62, uf, ha="center", va="center", fontsize=7.6, fontweight="bold", color=tc)
        ax.text(c + 0.5, -r + 0.32, f"{v:.2f}", ha="center", va="center", fontsize=6.6, color=tc)
    b = R["cf_nonppl"][per]
    ax.set_title(ttl, fontsize=8.5)
    ax.text(0, -6.0 + 0.05, f"Brazil\nO/E {b['OE']:.2f}\n({b['OE_lo']:.2f}-{b['OE_hi']:.2f})", fontsize=7, va="bottom", color=INK)
    ax.set_xlim(-0.1, 7.1); ax.set_ylim(-6.1, 1.05); ax.set_aspect("equal"); ax.axis("off")
# legenda
lx = fig.add_axes([0.2, 0.07, 0.6, 0.05]); lx.axis("off"); lx.set_xlim(0, 7); lx.set_ylim(0, 1)
for i, (c, l) in enumerate(zip(COLS, LAB)):
    lx.add_patch(Rectangle((i + 0.04, 0.45), 0.92, 0.45, fc=c, ec="white")); lx.text(i + 0.5, 0.2, l, ha="center", va="center", fontsize=6.6, color=INK2)
lx.text(3.5, 1.25, "Observed/expected new cases outside prisons (black outline: 95% prediction interval excludes 1)", ha="center", fontsize=6.8, color=INK2)
for ext, kw in [("png", {"dpi": 300}), ("pdf", {}), ("tif", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})]:
    fig.savefig(f"/home/claude/work/br/figv2/Figure2_v2.{ext}", bbox_inches="tight", facecolor="white", **kw)
