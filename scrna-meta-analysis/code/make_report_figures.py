# Publication figures for the galectin-7 OSCC scRNA-seq report.
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT = r"F:/oral_cancer_galectin-7"
RES = os.path.join(ROOT, "results", "scrna_galectin7")
OUT = os.path.join(RES, "report", "figures")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

TYPE_ORDER = [
    "Malignant", "Epithelial", "Fibroblast", "Endothelial", "Pericyte",
    "T cell", "B cell", "B/Plasma", "Macrophage", "Dendritic", "Myeloid",
    "Mast", "Myocyte", "Unassigned",
]
TYPE_COLOR = {
    "Malignant": "#9B2335",
    "Epithelial": "#E07A5F",
    "Fibroblast": "#2A9D8F",
    "Endothelial": "#457B9D",
    "Pericyte": "#1D3557",
    "T cell": "#6D597A",
    "B cell": "#4C956C",
    "B/Plasma": "#4C956C",
    "Macrophage": "#E9C46A",
    "Dendritic": "#F4A261",
    "Myeloid": "#E9C46A",
    "Mast": "#BC6C25",
    "Myocyte": "#6C757D",
    "Unassigned": "#ADB5BD",
}
GENES = ["LGALS1", "LGALS2", "LGALS3", "LGALS7", "LGALS7B", "LGALS8", "LGALS9"]
DATASETS = ["GSE181919", "GSE103322", "GSE172577"]
DS_LABEL = {
    "GSE181919": "GSE181919\nChoi, 10x",
    "GSE103322": "GSE103322\nPuram, Smart-seq2",
    "GSE172577": "GSE172577\nPeng, 10x",
}


def load_primary():
    frames = []
    for ds in DATASETS:
        df = pd.read_csv(os.path.join(RES, ds, "galectin_by_sample_celltype.csv"))
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    return df[df["tissue"] == "primary_tumor"].copy()


def weighted_pct(df):
    g = df.groupby(["dataset", "gene", "cell_type"], as_index=False).agg(
        n_cells=("n_cells", "sum"), n_detected=("n_detected", "sum")
    )
    g["pct"] = 100.0 * g["n_detected"] / g["n_cells"]
    return g


primary = load_primary()
pct = weighted_pct(primary)
lg7 = pct[pct["gene"] == "LGALS7"].copy()

# Figure 1. LGALS7 detection by cell type
fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.4), sharey=False)
for ax, ds in zip(axes, DATASETS):
    sub = lg7[lg7["dataset"] == ds].copy()
    sub = sub[sub["n_cells"] >= 50]
    sub["cell_type"] = pd.Categorical(sub["cell_type"], TYPE_ORDER, ordered=True)
    sub = sub.sort_values("pct")
    colors = [TYPE_COLOR.get(ct, "#888888") for ct in sub["cell_type"]]
    ax.barh(sub["cell_type"].astype(str), sub["pct"], color=colors, height=0.72)
    for y, (p, n) in enumerate(zip(sub["pct"], sub["n_cells"])):
        ax.text(p + 1.2, y, f"{p:.1f}%", va="center", fontsize=8, color="#333333")
    ax.set_xlim(0, 100)
    ax.set_title(DS_LABEL[ds], fontsize=11, loc="left", fontweight="medium")
    ax.set_xlabel("LGALS7 detected (%)")
    ax.axvline(0, color="#222", linewidth=0.6)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig1_lgals7_detection.png"), dpi=180, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig1_lgals7_detection.pdf"), bbox_inches="tight")
plt.close()

# Figure 2. Choi epithelial states
choi = pd.read_csv(os.path.join(RES, "GSE181919", "galectin_by_sample_celltype.csv"))
choi = choi[choi["gene"] == "LGALS7"]
states = [
    ("normal", "Epithelial", "Normal\nepithelium"),
    ("leukoplakia", "Epithelial", "Leukoplakia\nepithelium"),
    ("primary_tumor", "Malignant", "OSCC\nmalignant"),
]
rows = []
for tissue, ct, label in states:
    sub = choi[(choi["tissue"] == tissue) & (choi["cell_type"] == ct)]
    rows.append({
        "label": label,
        "pct": 100 * sub["n_detected"].sum() / sub["n_cells"].sum(),
        "n": int(sub["n_cells"].sum()),
        "n_samples": sub["sample_id"].nunique(),
    })
st = pd.DataFrame(rows)
fig, ax = plt.subplots(figsize=(5.2, 4.0))
cols = ["#4C6A7A", "#C4A35A", "#9B2335"]
bars = ax.bar(st["label"], st["pct"], color=cols, width=0.68)
for b, r in zip(bars, st.itertuples()):
    ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 2.2,
            f"{r.pct:.1f}%\n{r.n:,} cells\n{r.n_samples} samples",
            ha="center", va="bottom", fontsize=8, color="#333")
ax.set_ylim(0, 100)
ax.set_ylabel("LGALS7 detected (%)")
ax.set_title("GSE181919 oral cavity, epithelial compartment", loc="left", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig2_choi_epithelial_states.png"), dpi=180, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig2_choi_epithelial_states.pdf"), bbox_inches="tight")
plt.close()

# Figure 3. Galectin family heatmap
show_types = {
    "GSE181919": ["Malignant", "Fibroblast", "Endothelial", "T cell", "B/Plasma", "Macrophage"],
    "GSE103322": ["Malignant", "Fibroblast", "Endothelial", "T cell", "B cell", "Macrophage"],
    "GSE172577": ["Epithelial", "Fibroblast", "Endothelial", "T cell", "B/Plasma", "Myeloid"],
}
fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.6))
for ax, ds in zip(axes, DATASETS):
    types = show_types[ds]
    sub = pct[(pct["dataset"] == ds) & (pct["gene"].isin(GENES)) & (pct["cell_type"].isin(types))]
    mat = sub.pivot(index="gene", columns="cell_type", values="pct").reindex(index=GENES, columns=types)
    im = ax.imshow(mat.values, cmap="YlOrRd", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(types)))
    ax.set_xticklabels(types, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(GENES)))
    ax.set_yticklabels(GENES, fontsize=8)
    ax.set_title(DS_LABEL[ds].replace("\n", " · "), fontsize=10, loc="left")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.values[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=7,
                        color="white" if val > 55 else "#222")
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
fig.subplots_adjust(right=0.90, wspace=0.45)
cax = fig.add_axes([0.92, 0.18, 0.012, 0.64])
cb = fig.colorbar(im, cax=cax)
cb.set_label("Detected (%)", fontsize=9)
cb.outline.set_linewidth(0.4)
fig.savefig(os.path.join(OUT, "fig3_galectin_heatmap.png"), dpi=180, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig3_galectin_heatmap.pdf"), bbox_inches="tight")
plt.close()

# Figure 4. Program and selected-gene correlations
cor = pd.read_csv(os.path.join(RES, "trajectory", "LGALS7_cancer_program_correlations.csv"))
cor = cor[cor["subset"] == "all_tumor_epithelial"]
modules = [
    "basal", "differentiated", "epithelial_adhesion", "proliferation",
    "emt", "invasion", "oncogenic_signaling", "LGALS7B",
]
genes = ["KRT5", "KRT14", "KRT13", "IVL", "PKP1", "DSG3", "EPCAM",
         "MKI67", "LAMC2", "EGFR", "MYC", "CCND1", "CDKN2A", "LGALS3", "LGALS9"]
module_lab = {
    "basal": "Basal program",
    "differentiated": "Differentiated program",
    "epithelial_adhesion": "Epithelial adhesion",
    "proliferation": "Proliferation",
    "emt": "EMT",
    "invasion": "Invasion",
    "oncogenic_signaling": "Oncogenic signaling",
    "LGALS7B": "LGALS7B",
}

def rho_matrix(features):
    sub = cor[cor["feature"].isin(features)]
    mat = sub.pivot(index="feature", columns="dataset", values="rho").reindex(index=features, columns=DATASETS)
    return mat

def draw_rho(ax, mat, labels):
    norm = TwoSlopeNorm(vmin=-0.7, vcenter=0, vmax=0.7)
    im = ax.imshow(mat.values, cmap="RdBu_r", norm=norm, aspect="auto")
    ax.set_xticks(range(3))
    ax.set_xticklabels(["GSE181919", "GSE103322", "GSE172577"], fontsize=8)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.values[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if abs(val) > 0.38 else "#222")
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    return im

fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.6), gridspec_kw={"width_ratios": [1.05, 1.15]})
m1 = rho_matrix(modules)
im = draw_rho(axes[0], m1, [module_lab[f] for f in modules])
axes[0].set_title("Gene programs", loc="left", fontsize=11)
m2 = rho_matrix(genes)
draw_rho(axes[1], m2, genes)
axes[1].set_title("Selected genes", loc="left", fontsize=11)
fig.subplots_adjust(right=0.90, wspace=0.55)
cax = fig.add_axes([0.92, 0.18, 0.015, 0.64])
cb = fig.colorbar(im, cax=cax)
cb.set_label("Spearman rho with LGALS7", fontsize=9)
cb.outline.set_linewidth(0.4)
fig.savefig(os.path.join(OUT, "fig4_lgals7_correlations.png"), dpi=180, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig4_lgals7_correlations.pdf"), bbox_inches="tight")
plt.close()

# Figure 5. Pseudotime
pt_tests = pd.read_csv(os.path.join(RES, "trajectory", "LGALS7_pseudotime_tests.csv"))
fig, axes = plt.subplots(1, 3, figsize=(11.6, 3.9), sharey=False)
for ax, ds in zip(axes, DATASETS):
    df = pd.read_csv(os.path.join(RES, "trajectory", f"{ds}_pseudotime_cells.csv"))
    d = df[np.isfinite(df["pseudotime"])].copy()
    score = d["differentiated"] - d["basal"]
    # subsample for drawing
    rng = np.random.default_rng(7)
    if len(d) > 2500:
        idx = rng.choice(len(d), 2500, replace=False)
        d = d.iloc[idx]
        score = score.iloc[idx]
    sc = ax.scatter(d["pseudotime"], d["LGALS7"], c=score, cmap="RdBu_r",
                    s=8, alpha=0.45, linewidths=0, vmin=-1.2, vmax=1.2)
    # binned mean
    bins = np.linspace(d["pseudotime"].min(), d["pseudotime"].max(), 12)
    centers, means = [], []
    for a, b in zip(bins[:-1], bins[1:]):
        m = (d["pseudotime"] >= a) & (d["pseudotime"] < b)
        if m.sum() >= 15:
            centers.append((a + b) / 2)
            means.append(d.loc[m, "LGALS7"].mean())
    ax.plot(centers, means, color="#1A1A1A", linewidth=1.6)
    rho = pt_tests[(pt_tests.dataset == ds) & (pt_tests.feature == "LGALS7_vs_pseudotime")]["rho"].iloc[0]
    n = int(pt_tests[(pt_tests.dataset == ds) & (pt_tests.feature == "LGALS7_vs_pseudotime")]["n_cells"].iloc[0])
    ax.set_title(f"{ds}\nrho = {rho:.2f}, n = {n:,}", fontsize=10, loc="left")
    ax.set_xlabel("Slingshot pseudotime")
    ax.set_ylabel("LGALS7")
fig.subplots_adjust(right=0.90, wspace=0.32)
cax = fig.add_axes([0.92, 0.18, 0.012, 0.64])
cb = fig.colorbar(sc, cax=cax)
cb.set_label("Differentiated − basal", fontsize=8)
cb.outline.set_linewidth(0.4)
fig.savefig(os.path.join(OUT, "fig5_pseudotime.png"), dpi=180, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig5_pseudotime.pdf"), bbox_inches="tight")
plt.close()
print("figures written to", OUT)
