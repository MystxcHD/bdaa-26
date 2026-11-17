"""What do the top principal components line up with?"""
import pandas as pd

from parser import load_config
from metadata import load_metadata

cfg = load_config()
meta = load_metadata(cfg)
pcs = pd.read_csv("reports/pca_scores.csv", index_col=0)
qc = pd.read_csv("reports/sample_qc.csv", index_col=0)

df = meta.join(pcs).join(qc)
df["sex_male"] = (df["sex"] == "male").astype(int)
df["sample_number"] = df.index.str.extract(r"(\d+)", expand=False).astype(int)

traits = ["age", "sex_male", "sample_number", "missing_frac", "median_M"]
pc_cols = ["PC1", "PC2", "PC3", "PC4", "PC5"]

print("correlation of each PC with each trait:")
print(df[pc_cols + traits].corr().loc[pc_cols, traits].round(2).to_string())

print("mean PC score by diagnosis:")
print(df.groupby("diagnosis")[pc_cols].mean().round(1).to_string())