"""Check that the sex in the metadata agrees with sex-linked gene expression."""
import numpy as np
import pandas as pd

from parser import load_config, load_expression
from metadata import load_metadata

cfg = load_config()
M, probe_info = load_expression(cfg)
meta = load_metadata(cfg)


def probes_for(genes):
    pattern = r"(?:^|,)(?:%s)(?:,|$)" % "|".join(genes)
    return probe_info.index[probe_info["Gene"].fillna("").str.contains(pattern, regex=True)]


xist = probes_for(["XIST"])
ygenes = probes_for(["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY"])
print("XIST probes:", len(xist), " Y-gene probes:", len(ygenes))

# the pool mixes both sexes, so M is positive for the sex that has more of the gene
score = M.loc[xist].mean(axis=0) - M.loc[ygenes].mean(axis=0)  # > 0 looks female

df = meta.join(score.rename("female_score"), how="inner")
scored = df["female_score"].notna()
df["predicted"] = np.where(
    ~scored, "unscored", np.where(df["female_score"] > 0, "female", "male")
)
df["sex_mismatch"] = scored & (df["predicted"] != df["sex"])

mismatch = df[df["sex_mismatch"]]
print("sex mismatches:", len(mismatch), "of", len(df))
print("unscored samples:", (~scored).sum())
print(mismatch[["gsm", "diagnosis", "sex", "female_score"]].to_string())

df[["gsm", "diagnosis", "sex", "age", "female_score", "predicted", "sex_mismatch"]].to_csv(
    "reports/sex_check.csv"
)