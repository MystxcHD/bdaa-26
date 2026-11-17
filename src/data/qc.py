"""Per-sample and per-probe quality tables for the log-ratio matrix."""
import numpy as np
import pandas as pd

from parser import load_config, load_expression


def sample_qc(M):
    vals = M.to_numpy()
    return pd.DataFrame(
        {
            "median_M": np.nanmedian(vals, axis=0),
            "missing_frac": np.isnan(vals).mean(axis=0),
        },
        index=M.columns,
    )


def probe_qc(M):
    return pd.DataFrame(
        {"missing_frac": np.isnan(M.to_numpy()).mean(axis=1)},
        index=M.index,
    )


if __name__ == "__main__":
    cfg = load_config()
    M, probe_info = load_expression(cfg)

    s = sample_qc(M)
    p = probe_qc(M)
    s.to_csv("reports/sample_qc.csv")
    p.to_csv("reports/probe_qc.csv")

    print(s.describe())
    order = s["median_M"].abs().sort_values(ascending=False).index
    print("10 most off-center samples:")
    print(s.loc[order].head(10))

    print("probes above each missing-value level:")
    for t in [0, 0.01, 0.05, 0.10, 0.20, 0.50]:
        print(f"  more than {t:.0%} missing: {(p['missing_frac'] > t).sum()} probes")