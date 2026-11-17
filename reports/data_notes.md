# GSE33000 data notes (Member 1, Sprint 1)

Status: processed dataset v1 built and verified. Several decisions are open for the team (bottom).

## What the data is
- 624 post-mortem dorsolateral prefrontal cortex (BA9) samples: 310 Alzheimer's (AD), 157 Huntington's (HD), 157 non-demented controls.
- Platform: Rosetta/Merck Human 44k 1.1 (GPL4372), two-color. Each sample was hybridized against a shared pool made from the same samples (the pool mixes AD, HD and control tissue, so 0 means "cohort average", not "healthy").
- `GSE33000_raw_data.txt` holds the per-channel intensities (pool and sample). The expression value used here is the log2 ratio sample / pool, computed per probe per array. It is not normalized beyond that; per-array medians are already close to 0 (range -0.07 to +0.10).
- Only age, sex and diagnosis are available as sample information (confirmed against the SDRF and the investigation design file). No RIN, post-mortem interval, pH or batch.

## Processed files (data/processed, not tracked by git)
- `expression_samples_x_genes.parquet`: 624 samples x 19,031 genes, float32, 1.23% of values missing (left as NaN on purpose).
- `sample_metadata.parquet`: diagnosis, sex, age, GSM id, per-sample QC numbers and flags.
- `gene_info.csv`: which probe was kept for each gene, how many probes the gene had, symbol, Entrez ID.
- `reports/probe_decisions.csv`: every one of the 39,280 probes and why it was kept or dropped.
- Rebuild everything with `python src/data/preprocess.py`. Reading the parquet needs `python -m pip install pyarrow`.

## What was done (all steps use no diagnosis information)
1. Probes missing in more than 10% of samples dropped: 686.
2. Probes with no gene symbol in GEO's GPL4372 annotation dropped: 15,029.
3. Probes matching several genes ("///" in the annotation) dropped: 125.
4. One probe kept per gene (Entrez ID), the one with the largest interquartile range: 4,409 duplicate probes dropped.
5. Result: 19,031 genes (15,345 backed by one probe, 3,686 by two or more; 1,072 are LOC-named placeholder genes).
Checked: 18 values (6 genes x 3 samples) match the raw file exactly; XIST is high in recorded females and low in recorded males.

## Findings the team should know
- Age differs a lot by diagnosis (mean, range): AD 80.6 (53-100), HD 55.9 (18-85), controls 63.5 (22-106).
- Sex differs a lot by diagnosis (female/male): AD 175/135, HD 74/83, controls 34/123 (controls are 78% male).
- Both are confounds: a model can partly separate groups using age or sex-linked genes. Suggested: an age-only and a sex-only baseline classifier as benchmarks.
- Sex check from XIST and Y-gene expression: 9 samples disagree with recorded sex (PFC_601, 561, 395, 285, 260, 172, 123, 76, 52). One male sample could not be scored. These may be mislabeled or swapped samples, and the diagnosis label could be wrong for them too. They are flagged in the metadata (`sex_mismatch`), not removed.
- Missing values per sample: mean 1.2%, max 8.7%. By diagnosis: AD 1.3%, HD 0.8%, controls 1.5%.
- PCA (QC only): PC1 (33.7%) separates controls from AD/HD, PC2 (10.1%) separates HD from the rest. No isolated outlier samples (largest robust z about 4.5).
- Sample numbers are sorted by diagnosis (AD 1-310, controls 311-467, HD 468-624). Never use them as a feature or a batch proxy.
- Within each diagnosis, PC2 still trends with sample number (r = -0.65 in AD, -0.63 in controls, +0.37 in HD). It is not explained by age, sex or missingness. Cause unknown (possibly processing order), and it cannot be checked with the available files.

## Known limitations of v1
- The GEO annotation is from 2007, so some symbols are old (the Huntington gene appears as `HD`, not `HTT`). Use the Entrez ID in gene_info.csv when converting to current symbols for enrichment.
- Probes without a symbol in that annotation were dropped even where the array's own gene label has one. Example: PSEN1 (3 probes labeled PSEN1 in the data file) is missing from v1. About 8,500 such probes look like normal gene symbols (some are accession-number junk).

## Open decisions for the team
1. Exclude the 9 sex-mismatch samples (and the unscored one)? Member 1 suggests excluding them for modeling and noting it.
2. Keep or drop the LOC placeholder genes (config: probes.drop_loc_genes).
3. Add a flagged second gene source (the array's own labels) to recover genes like PSEN1? Would add about 5,000 less certain genes, marked by source.
4. Center each sample at its own median (label-free)? It would remove a small technical axis (PC3) without touching diagnosis.
5. How to handle age and sex (Members 2 and 3). Any adjustment must be fit inside the training folds. The GPL4372 annotation has chromosome locations if X/Y genes should be flagged.
6. Fill the remaining 1.2% gaps inside each training fold only.
