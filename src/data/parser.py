"""Load the GSE33000 expression file and compute log2 ratios (sample vs pool)."""
import numpy as np
import pandas as pd
import yaml


def load_config(path="configs/data_config.yaml"):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def load_expression(cfg):
    """Return (M, probe_info).

    M: DataFrame, rows = probes (reporterID), columns = PFC_1 ... PFC_N,
       values = log2(sample / pool). Failed spots stay as NaN.
    probe_info: DataFrame indexed by reporterID with Gene and Transcript.
    """
    path = cfg["paths"]["expression_file"]

    header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    pool_cols = [c for c in header if c.startswith("HBTRC_PF_Pool_")]
    pfc_cols = [c for c in header if c.startswith("PFC_")]
    n = len(pfc_cols)
    assert len(pool_cols) == n, "pool and sample column counts differ"

    dtypes = {c: np.float32 for c in pool_cols + pfc_cols}
    dtypes["reporterID"] = str
    df = pd.read_csv(path, sep="\t", dtype=dtypes)

    # pair PFC_i with HBTRC_PF_Pool_i by number
    pool = df[[f"HBTRC_PF_Pool_{i}" for i in range(1, n + 1)]].to_numpy()
    pfc = df[[f"PFC_{i}" for i in range(1, n + 1)]].to_numpy()
    M = pd.DataFrame(
        np.log2(pfc / pool),
        index=df["reporterID"],
        columns=[f"PFC_{i}" for i in range(1, n + 1)],
    )
    probe_info = df[["reporterID", "Gene", "Transcript"]].set_index("reporterID")
    return M, probe_info


if __name__ == "__main__":
    cfg = load_config()
    M, probe_info = load_expression(cfg)
    print("probes x samples:", M.shape)
    print("probe IDs unique:", M.index.is_unique)
    print("missing fraction:", float(np.isnan(M.to_numpy()).mean()))
    print("overall median:", float(np.nanmedian(M.to_numpy())))
    print(probe_info.head())