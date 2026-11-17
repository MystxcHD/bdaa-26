"""Turn the probe-level log-ratio matrix into a samples x genes matrix.

Every step here is label-free (no diagnosis information is used), so it is safe to run
once on all samples. Imputation of the remaining gaps and any supervised feature selection
must still happen inside the training folds.
"""
import numpy as np
import pandas as pd

from parser import load_config, load_expression
from metadata import load_metadata


def load_annotation(cfg):
    """Probe -> gene symbol and Entrez Gene ID, from GEO's GPL4372 annotation table."""
    p = cfg["platform"]
    a = pd.read_csv(
        cfg["paths"]["platform_annot"],
        sep="\t",
        skiprows=p["annot_header_lines"],
        dtype=str,
        quoting=3,
    )
    a = a[~a[p["id_col"]].str.startswith("!")]
    a = a.rename(columns={p["id_col"]: "probe", p["symbol_col"]: "symbol", p["entrez_col"]: "entrez_id"})
    return a[["probe", "symbol", "entrez_id"]].set_index("probe")


def probe_table(M, cfg):
    """One row per probe: its gene, its quality numbers, and why it was dropped (if it was)."""
    pc = cfg["probes"]
    t = load_annotation(cfg).reindex(M.index)
    t["symbol"] = t["symbol"].fillna("").str.strip()
    t["entrez_id"] = t["entrez_id"].fillna("").str.strip()
    t["missing_frac"] = M.isna().mean(axis=1).to_numpy()
    q = np.nanpercentile(M.to_numpy(), [25, 75], axis=1)
    t["iqr"] = q[1] - q[0]
    t["reason_dropped"] = ""

    def drop(mask, why):
        t.loc[mask & (t["reason_dropped"] == ""), "reason_dropped"] = why

    drop(t["missing_frac"] > pc["max_missing_frac"], "too_many_missing")
    drop(t["symbol"] == "", "no_gene_symbol")
    if pc["drop_ambiguous"]:
        drop(t["symbol"].str.contains("///") | t["entrez_id"].str.contains("///"), "ambiguous_gene")
    if pc["drop_loc_genes"]:
        drop(t["symbol"].str.startswith("LOC"), "loc_gene")
    return t


def collapse_to_genes(M, t):
    """Keep one probe per gene (Entrez ID): the one with the largest IQR across samples."""
    assert t is not None
    keep = t[t["reason_dropped"] == ""].copy()
    keep = keep.sort_values(["entrez_id", "iqr"], ascending=[True, False])
    n_probes = keep.groupby("entrez_id").size()
    chosen = keep.drop_duplicates("entrez_id", keep="first").copy()
    chosen["n_probes"] = chosen["entrez_id"].map(n_probes)
    # column names are gene symbols; if one symbol covers two Entrez IDs, add the ID to tell them apart
    dup = chosen["symbol"].duplicated(keep=False)
    chosen["gene"] = np.where(dup, chosen["symbol"] + "_" + chosen["entrez_id"], chosen["symbol"])
    X = M.loc[chosen.index].T
    X.columns = chosen["gene"].to_numpy()
    return X, chosen


def build_metadata(sample_ids):
    """Diagnosis, age, sex, plus the QC flags from earlier steps (flags are recorded, not applied)."""
    cfg = load_config()
    meta = load_metadata(cfg)
    qc = pd.read_csv("reports/sample_qc.csv", index_col=0)
    sx = pd.read_csv("reports/sex_check.csv", index_col=0)[["female_score", "predicted", "sex_mismatch"]]
    sx = sx.rename(columns={"predicted": "sex_from_expression"})
    out = meta.join(qc).join(sx)
    out["sex_mismatch"] = out["sex_mismatch"].fillna(False).astype(bool)
    out["sex_unscored"] = out["sex_from_expression"].eq("unscored")
    return out.loc[sample_ids]


def build(cfg):
    assert cfg["probes"]["collapse_method"] == "max_iqr", "only max_iqr is implemented"
    M, _ = load_expression(cfg)
    t = probe_table(M, cfg)
    X, chosen = collapse_to_genes(M, t)
    t["kept"] = t.index.isin(chosen.index)
    t.loc[(t["reason_dropped"] == "") & ~t["kept"], "reason_dropped"] = "duplicate_probe_for_gene"
    meta = build_metadata(X.index)
    excluded = cfg["samples"]["exclude_samples"]
    if excluded:
        X = X.drop(index=excluded)
        meta = meta.drop(index=excluded)
    return X, meta, chosen, t


def save(cfg, X, meta, chosen, t):
    o = cfg["output"]
    X.astype("float32").to_parquet(o["matrix_file"])
    meta.to_parquet(o["metadata_file"])
    info = chosen[["gene", "symbol", "entrez_id", "n_probes", "missing_frac", "iqr"]].copy()
    info.insert(0, "probe", chosen.index)
    info.to_csv(o["gene_info_file"], index=False)
    t.to_csv(o["probe_decisions_file"])


if __name__ == "__main__":
    cfg = load_config()
    X, meta, chosen, t = build(cfg)
    save(cfg, X, meta, chosen, t)

    print("probes by outcome:")
    print(t["reason_dropped"].replace("", "kept (one per gene)").value_counts().to_string())
    print()
    print("matrix (samples x genes):", X.shape)
    print("missing fraction:", round(float(np.isnan(X.to_numpy()).mean()), 4))
    print("column names unique:", X.columns.is_unique)
    print(meta["diagnosis"].value_counts().to_string())
    print("sex mismatches flagged:", int(meta["sex_mismatch"].sum()), " unscored:", int(meta["sex_unscored"].sum()))
    for gname in ["XIST", "HTT", "APP", "MAPT", "GFAP", "SNCA"]:
        hit = chosen[chosen["symbol"] == gname]
        print(gname, "->", "present, probes for it: " + str(int(hit["n_probes"].iloc[0])) if len(hit) else "not in matrix")
