"""QC-only PCA: look for outlier samples and for structure beyond diagnosis."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA

from parser import load_config, load_expression
from metadata import load_metadata

cfg = load_config()
M, probe_info = load_expression(cfg)
meta = load_metadata(cfg)

# QC-only preparation (no labels used): drop probes with >10% missing,
# fill the remaining gaps with each probe's median, keep the 5000 most variable probes
keep = M.isna().mean(axis=1) <= 0.10
Xs = M.loc[keep].T
Xs = Xs.fillna(Xs.median())
Xs = Xs[Xs.var().nlargest(5000).index]

pca = PCA(n_components=10, random_state=cfg["seed"])
pcs = pd.DataFrame(
    pca.fit_transform(Xs), index=Xs.index, columns=[f"PC{i}" for i in range(1, 11)]
)
pcs.to_csv("reports/pca_scores.csv")
print("variance explained by PC1-PC5 (%):", (pca.explained_variance_ratio_[:5] * 100).round(1))

# which samples sit furthest from the crowd in the first 5 PCs?
five = pcs[["PC1", "PC2", "PC3", "PC4", "PC5"]]
centered = five - five.median()
z = centered / (1.4826 * centered.abs().median())
extreme = z.abs().max(axis=1).sort_values(ascending=False).head(10)
out = pd.concat(
    [extreme.rename("max_robust_z"), meta.loc[extreme.index, ["diagnosis", "sex", "age"]]],
    axis=1,
)
print("10 most extreme samples:")
print(out.round(1).to_string())


def scatter(by, palette, path):
    fig, ax = plt.subplots(figsize=(6, 5))
    for label, color in palette.items():
        idx = meta.index[meta[by] == label].intersection(pcs.index)
        ax.scatter(pcs.loc[idx, "PC1"], pcs.loc[idx, "PC2"], s=14, c=color, label=label, alpha=0.8)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


scatter(
    "diagnosis",
    {"Alzheimer's disease": "#D95F02", "Huntington's disease": "#7570B3", "non-demented": "#1B9E77"},
    "reports/figures/pca_diagnosis.png",
)
scatter(
    "sex",
    {"female": "#CC79A7", "male": "#0072B2"},
    "reports/figures/pca_sex.png",
)