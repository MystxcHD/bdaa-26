"""Is sample number (processing order?) tangled up with diagnosis?"""
import pandas as pd

from parser import load_config
from metadata import load_metadata

cfg = load_config()
meta = load_metadata(cfg)
pcs = pd.read_csv("reports/pca_scores.csv", index_col=0)

df = meta.join(pcs)
df["sample_number"] = df.index.str.extract(r"(\d+)", expand=False).astype(int)

# 1. where do the diagnoses sit along the sample numbers?
df["block"] = pd.cut(df["sample_number"], bins=[0, 100, 200, 300, 400, 500, 624])
print(pd.crosstab(df["block"], df["diagnosis"]))

# 2. does PC2 still track sample number within each diagnosis?
for dx, g in df.groupby("diagnosis"):
    print(dx, "PC2 vs sample_number r =", round(g["PC2"].corr(g["sample_number"]), 2))
