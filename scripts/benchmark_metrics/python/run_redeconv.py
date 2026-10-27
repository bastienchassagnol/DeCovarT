"""Run the official ReDeconv implementation on a batch of scenarios.

Usage: python run_redeconv.py <in_dir> <out_dir>

<in_dir>/manifest.csv lists `scenario_id`. Each <in_dir>/<scenario_id>/ holds
`bulk.csv` (genes x bulk samples, first column gene ids), `celltypes.csv`
(column `celltype`) and the plug-in reference moments `signature_mean.csv`
and `signature_sd.csv` (genes x cell types).

The moments are written in the signature format read by
`redeconv.__ReDeconv_P.read_sig_gene_mean_std()` (tab-separated, one
`<celltype>_mean` / `<celltype>_std` column pair per type) and the bulk in the
format read by `read_bulk_RNA_seq_column_for_sample()`. The default
`ReDeconv()` entry point (MLMD, six reweighting loops) is then called. Writing
the signature from the given moments bypasses ReDeconv's signature-gene
selection (t-tests, fold change, top-k by coefficient of variation): every
gene is used. ReDeconv reports no confidence interval.

Writes <out_dir>/<scenario_id>/estimates.csv (cell types x bulk columns) and
<out_dir>/status.csv (scenario_id, status, message, elapsed_sec).
"""

import argparse
import contextlib
import importlib
import io
import os
import time
import traceback
from pathlib import Path

import pandas as pd

os.environ.setdefault("MPLBACKEND", "Agg")
redeconv_p = importlib.import_module("redeconv.__ReDeconv_P")


def write_signature(mean, sd, path):
    columns = {}
    for celltype in mean.columns:
        columns[f"{celltype}_mean"] = mean[celltype]
        columns[f"{celltype}_std"] = sd[celltype]
    signature = pd.DataFrame(columns, index=mean.index)
    signature.index.name = "gene"
    signature.to_csv(path, sep="\t", float_format="%.17g")


def fit_scenario(scenario_dir, work_dir):
    bulk = pd.read_csv(scenario_dir / "bulk.csv", index_col=0)
    celltypes = pd.read_csv(scenario_dir / "celltypes.csv")["celltype"].tolist()
    mean = pd.read_csv(scenario_dir / "signature_mean.csv", index_col=0)
    sd = pd.read_csv(scenario_dir / "signature_sd.csv", index_col=0)
    if list(mean.index) != list(bulk.index) or list(sd.index) != list(bulk.index):
        raise ValueError("signature genes differ from bulk.csv genes")
    if list(mean.columns) != celltypes or list(sd.columns) != celltypes:
        raise ValueError("signature columns differ from celltypes.csv")

    signature_file = work_dir / "signature.tsv"
    bulk_file = work_dir / "bulk.tsv"
    result_file = work_dir / "percentages.tsv"
    write_signature(mean, sd, signature_file)
    bulk_out = bulk.copy()
    bulk_out.index.name = "gene"
    bulk_out.to_csv(bulk_file, sep="\t", float_format="%.17g")

    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        redeconv_p.ReDeconv(str(signature_file), str(bulk_file), str(result_file))

    result = pd.read_csv(result_file, sep="\t", index_col=0)
    estimates = result.T.loc[celltypes, list(bulk.columns)]
    return estimates


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("in_dir", type=Path)
    parser.add_argument("out_dir", type=Path)
    args = parser.parse_args()

    manifest = pd.read_csv(args.in_dir / "manifest.csv")
    status = []
    for scenario_id in manifest["scenario_id"].astype(str):
        out_dir = args.out_dir / scenario_id
        out_dir.mkdir(parents=True, exist_ok=True)
        start = time.perf_counter()
        try:
            estimates = fit_scenario(args.in_dir / scenario_id, out_dir)
            estimates.to_csv(out_dir / "estimates.csv")
            status.append((scenario_id, "ok", "", time.perf_counter() - start))
        except Exception as error:  # record the failure, no fallback fit
            message = f"{type(error).__name__}: {error}"
            traceback.print_exc()
            status.append((scenario_id, "error", message, time.perf_counter() - start))

    pd.DataFrame(
        status, columns=["scenario_id", "status", "message", "elapsed_sec"]
    ).to_csv(args.out_dir / "status.csv", index=False)


if __name__ == "__main__":
    main()
