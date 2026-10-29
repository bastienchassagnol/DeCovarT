"""Run the official RNA-Sieve implementation on a batch of scenarios.

Usage: python run_rnasieve.py <in_dir> <out_dir> [--level 0.95] [--workers 1]

<in_dir>/manifest.csv lists `scenario_id`. Each <in_dir>/<scenario_id>/ holds
`bulk.csv` (genes x bulk samples, first column gene ids), `celltypes.csv`
(column `celltype`, in output order) and one `reference_<k>.csv` per cell type
k = 1..J (genes x reference cells, same gene order as the bulk).

The fit uses `rnasieve.preprocessing.model_from_raw_counts()` with its
preprocessing switched off through its own arguments (trim_percent = 0: no
cell trimming; normalization = False: no library-size rescaling;
gene_thresh = 0: no gene filter), then `RNASieveModel.predict()` jointly over
all bulk columns of the scenario (the package default).

Three marginal intervals per cell type and bulk column:

* `lower.csv`/`upper.csv`: the released
  `RNASieveModel.compute_marginal_confidence_intervals()`. When the averaged
  Hessian has a condition number of 1e6 or more, rnasieve 0.1.4 subsamples
  genes and returns H B H' instead of H^-1 B H^-T (stats.py:224), and the
  interval collapses to [0, 1].
* `lower_godambe.csv`/`upper_godambe.csv`: the Godambe (sandwich) interval of
  the paper's "Confidence intervals" section, H^-1 B H^-T / G with H the
  per-gene Hessian averaged over all genes and B the covariance of the
  per-gene scores, both from rnasieve's own functions, without subsampling.
* `lower_fisher.csv`/`upper_fisher.csv`: the Wald interval from the inverse
  empirical Fisher information, B^-1 / G, which the Godambe matrix reduces to
  without misspecification (H = -B). The averaged Hessian itself is not used
  here: at an estimate on a face of the simplex it is indefinite.

The parameters are (alpha_1, ..., alpha_{K-1}, n), with alpha_K = 1 - sum of
the others; the variance of alpha_K is therefore 1' V 1 over the first K - 1
coordinates. Intervals are alpha +/- z sqrt(var), clipped to [0, 1] as in
rnasieve.

Writes <out_dir>/<scenario_id>/{estimates,lower,upper,lower_fisher,
upper_fisher,lower_godambe,upper_godambe}.csv (cell types x bulk columns) and
<out_dir>/status.csv (scenario_id, status, message, elapsed_sec).
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from rnasieve.preprocessing import model_from_raw_counts
from rnasieve.stats import (
    _godambe_hessian_log_likelihood_alpha,
    _godambe_partial_log_likelihood_alpha,
)
from scipy.stats import norm

from batch import run_batch


def marginal_variances(alpha, phi, sigma, psi, n):
    """Fisher and Godambe variances of (alpha_1, ..., alpha_K) for one bulk."""
    rows = range(phi.shape[0])
    scores = np.vstack(
        [
            _godambe_partial_log_likelihood_alpha(
                alpha, phi[i : i + 1], sigma[i : i + 1], psi[i : i + 1], n
            )
            for i in rows
        ]
    )
    hessian = np.mean(
        [
            _godambe_hessian_log_likelihood_alpha(
                alpha, phi[i : i + 1], sigma[i : i + 1], psi[i : i + 1], n
            )
            for i in rows
        ],
        axis=0,
    )
    n_genes = phi.shape[0]
    score_cov = np.cov(scores.T)
    inv_hessian = np.linalg.inv(hessian)
    godambe = inv_hessian @ score_cov @ inv_hessian.T / n_genes
    fisher = np.linalg.inv(score_cov) / n_genes

    def by_type(v):
        v = v[:-1, :-1]
        return np.append(np.diag(v), np.sum(v))

    return by_type(fisher), by_type(godambe)


def fit_scenario(scenario_dir, level):
    bulk = pd.read_csv(scenario_dir / "bulk.csv", index_col=0)
    celltypes = pd.read_csv(scenario_dir / "celltypes.csv")["celltype"].tolist()
    raw_counts = {}
    for k, celltype in enumerate(celltypes, start=1):
        reference = pd.read_csv(scenario_dir / f"reference_{k}.csv", index_col=0)
        if list(reference.index) != list(bulk.index):
            raise ValueError(f"reference_{k}.csv genes differ from bulk.csv genes")
        raw_counts[celltype] = reference.to_numpy(dtype=float)

    model, psis = model_from_raw_counts(
        raw_counts,
        bulk.to_numpy(dtype=float),
        trim_percent=0.0,
        gene_thresh=0.0,
        normalization=False,
    )
    if psis.shape[0] != bulk.shape[0]:
        raise ValueError(
            f"RNA-Sieve dropped {bulk.shape[0] - psis.shape[0]} all-zero genes"
        )
    alpha_hats = model.predict(psis)
    cis = model.compute_marginal_confidence_intervals(sig=1.0 - level)

    labels = list(alpha_hats.columns)
    z = norm.isf((1.0 - level) / 2)
    phi = model.phi_hat.to_numpy()
    sigma = model.observed_sigma.to_numpy()[model.filter_idxs]
    psi_all = model.psis.to_numpy()[model.filter_idxs]
    tables = {"lower": [], "upper": []}
    for kind in ("fisher", "godambe"):
        tables[f"lower_{kind}"] = []
        tables[f"upper_{kind}"] = []
    for i in range(alpha_hats.shape[0]):
        alpha = alpha_hats.to_numpy()[i].reshape(1, -1)
        tables["lower"].append([ci[0] for ci in cis[i]])
        tables["upper"].append([ci[1] for ci in cis[i]])
        try:
            variances = marginal_variances(
                alpha,
                phi,
                sigma,
                psi_all[:, i].reshape(-1, 1),
                model.n_hats.to_numpy()[i],
            )
        except np.linalg.LinAlgError:
            variances = (np.full(len(labels), np.nan),) * 2
        for kind, var in zip(("fisher", "godambe"), variances):
            half = z * np.sqrt(np.where(var >= 0, var, np.nan))
            tables[f"lower_{kind}"].append(np.clip(alpha[0] - half, 0.0, 1.0))
            tables[f"upper_{kind}"].append(np.clip(alpha[0] + half, 0.0, 1.0))

    def frame(values):
        out = pd.DataFrame(np.asarray(values).T, index=labels, columns=bulk.columns)
        return out.loc[celltypes]

    estimates = pd.DataFrame(
        alpha_hats.to_numpy().T, index=labels, columns=bulk.columns
    )
    return {"estimates": estimates.loc[celltypes]} | {
        name: frame(values) for name, values in tables.items()
    }


def fit_and_write(scenario_dir, out_dir, level):
    for name, table in fit_scenario(scenario_dir, level).items():
        table.to_csv(out_dir / f"{name}.csv")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("in_dir", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--level", type=float, default=0.95)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    run_batch(
        "RNA-Sieve",
        fit_and_write,
        args.in_dir,
        args.out_dir,
        workers=args.workers,
        level=args.level,
    )


if __name__ == "__main__":
    main()
