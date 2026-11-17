"""Join metadata with sample QC and check for confounding."""
import pandas as pd

from parser import load_config
from metadata import load_metadata

cfg = load_config()
meta = load_metadata(cfg)
qc = pd.read_csv("reports/sample_qc.csv", index_col=0)

# 1. do the sample IDs line up exactly?
print("IDs only in metadata:", sorted(set(meta.index) - set(qc.index))[:5])
print("IDs only in expression:", sorted(set(qc.index) - set(meta.index))[:5])

df = meta.join(qc, how="inner")
print("joined samples:", len(df))

# 2. sex by diagnosis
print(pd.crosstab(df["diagnosis"], df["sex"]))

# 3. age by diagnosis (full summary)
print(df.groupby("diagnosis")["age"].describe().round(1).to_string())

# 4. technical QC by diagnosis
print(df.groupby("diagnosis")[["median_M", "missing_frac"]].mean().round(4))