import pandas as pd, numpy as np, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

Q = pd.read_csv("/home/claude/work/br/qualidade_uf.csv", index_col=0)
cols = [("casos", "Casos\nnovos", None, "n"), ("sem_encerr_%", "Sem\nencerra-\nmento*", False, "%"),
        ("raca_ign_%", "Raça/cor\nignorada", False, "%"), ("escol_ign_%", "Escola-\nridade\nignorada", False, "%"),
        ("ppl_ign_%", "PPL\nignorada", False, "%"), ("hiv_result_%", "HIV com\nresultado", True, "%"),
        ("trm_%", "TRM-TB\nrealizado", True, "%"), ("cultura_%", "Cultura\nrealizada", True, "%"),
        ("transf_%", "Encerr.\ncomo trans-\nferência*", False, "%"), ("ppl_%", "Casos\nem PPL", None, "%"),
        ("cura_%", "Cura*", True, "%"), ("abandono_%", "Aban-\ndono*", False, "%")]
mpl.rcParams.update({"font.family": "Liberation Sans", "font.size": 7.5})
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
        txt = f"{val:,.0f}".replace(",", ".") if u == "n" else f"{val:.1f}".replace(".", ",")
        ax.text(j + 0.5, i + 0.5, txt, ha="center", va="center", fontsize=6.8,
                color="white" if x > 0.7 else "#0b0b0b", fontweight="bold" if uf == "Brasil" else "normal")
ax.set_xticks(np.arange(len(cols)) + 0.5); ax.set_xticklabels([c[1] for c in cols], fontsize=6.9); ax.xaxis.tick_top()
ax.set_yticks(np.arange(len(Q)) + 0.5); ax.set_yticklabels(Q.index, fontsize=7.5)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0); ax.axhline(len(Q) - 1, color="#0b0b0b", lw=0.8)
ax.text(0, len(Q) + 0.9, "Casos novos de TB pulmonar com confirmação laboratorial, 2015-2025, por UF de residência (%). Laranja mais intenso: pior qualidade\n"
        "ou pior desfecho; cinza: variáveis descritivas. Escala relativa entre UFs em cada coluna. *Coorte 2015-2024. Fonte: Sinan/DATASUS, extração em 24/11/2026.",
        fontsize=6.6, color="#52514e", va="top")
fig.savefig("/home/claude/work/fig_br/qualidade_uf.png", dpi=200, bbox_inches="tight", facecolor="white")
