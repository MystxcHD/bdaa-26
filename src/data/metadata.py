"""Load sample metadata (diagnosis, age, sex) from the ArrayExpress SDRF."""
import pandas as pd

from parser import load_config


def load_metadata(cfg):
    sdrf = pd.read_csv(cfg["paths"]["sdrf"], sep="\t", dtype=str)
    # each GSM has two rows: Cy3 = pool, Cy5 = the brain sample. Keep the brain rows.
    sdrf = sdrf[sdrf["Label"] == cfg["sdrf"]["keep_label"]].copy()

    meta = pd.DataFrame(
        {
            "sample_id": sdrf["Comment [Sample_title]"].str.strip(),
            "gsm": sdrf["Source Name"].str.split().str[0],
            "diagnosis": sdrf["Characteristics [disease status]"].str.strip(),
            "sex": sdrf["Characteristics [sex]"].str.strip(),
            "age_raw": sdrf["Characteristics [age]"].str.strip(),
        }
    )
    meta["age"] = pd.to_numeric(meta["age_raw"].str.extract(r"(\d+)")[0])
    return meta.set_index("sample_id")


if __name__ == "__main__":
    cfg = load_config()
    meta = load_metadata(cfg)

    print("samples:", len(meta), " IDs unique:", meta.index.is_unique)
    print(meta["diagnosis"].value_counts(dropna=False))
    print(meta["sex"].value_counts(dropna=False))

    odd = meta[~meta["age_raw"].fillna("").str.match(r"^\d+ yrs$")]
    print("ages not written like '57 yrs':", len(odd))
    print(odd["age_raw"].value_counts(dropna=False).head(10))

    print(meta.groupby("diagnosis")["age"].describe())