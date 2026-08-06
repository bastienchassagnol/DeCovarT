# Overlap Measures in Gaussian Mixtures

For two multivariate Gaussians (components $i$ and $j$ with densities $f_i,f_j$ and mixing weights $\pi_i,\pi_j$), the **theoretical overlap** can be defined as the Bayes classification error:
\[
\eta_{ij} \;=\; \int \min\{\pi_i f_i(x),\,\pi_j f_j(x)\}\,dx.
\]
This is equivalently the expected misclassification rate of the optimal (likelihood‐ratio) classifier.  In practice, especially with unequal covariance matrices or higher dimension, this integral has no closed form and is *intractable* to compute exactly.  For example, for $d>1$ and distinct $\Sigma_i,\Sigma_j$, one would need to integrate a Gaussian density over complicated regions (the decision boundary of the classifier).  As a result, modern tools use *pairwise* overlap approximations.

**Pairwise overlap (misclassification).**  The R package **mixsim** defines the *pairwise overlap* between two components as the sum of their misclassification probabilities:
\[
\text{Overlap}_{ij} \;=\; P(\text{class }i\to j) + P(\text{class }j\to i),
\]
where $P(i\to j)$ is the probability a point from component $i$ is (Bayes) classified as $j$.  Equivalently, if priors are equal, $\text{Overlap}_{ij}=2\int \min\{f_i(x),f_j(x)\}dx$.  MixSim’s `overlap()` function computes the full $K\times K$ **OmegaMap** of these misclassification probabilities given any means `Mu` and covariance arrays `S` (so it handles distinct covariances).  In the output of `overlap(Pi,Mu,S)`, **BarOmega** is defined as the average of all pairwise overlaps, and **MaxOmega** as the largest overlap between any pair.  (The mixsim *MixSim()* function similarly returns `BarOmega` and `MaxOmega` along with the parameters, which can be used to simulate mixtures at a desired overlap.)  For unequal covariance Gaussians, mixsim internally uses numerical integration (via the Davies algorithm for quadratic forms) to evaluate each $P(i\to j)$.

**Global vs pairwise overlap.**  In principle one could define a *global overlap* for a $K$-component mixture as the overall Bayes error of distinguishing all $K$ classes: e.g.\
\[
\eta_{\rm global} \;=\; 1 - \int \max_{1\le k\le K}\{\pi_k f_k(x)\}\,dx.
\]
However, this multi-class integral is highly complicated in $d>1$, especially with heteroscedastic covariances.  There is no simple analytic formula for $\eta_{\rm global}$ when $K>2$, so standard practice (as in mixsim) is to use *pairwise* overlaps as proxies.  Hence mixsim reports the average or maximum pairwise overlap rather than a single global value.  In effect, the average BarOmega is a crude measure of overall class separation; the maximum overlap identifies the two “closest” clusters.

**Mixture (latent) vs complete (labeled) case.**  The above discussion assumes a mixture model where class labels are unobserved.  If instead one has a *complete dataset with known labels* for every point, one can directly compute classification error.  In that case the overlap is just the misclassification rate of a classifier (e.g.\ Gaussian Bayes or linear discriminant).  The theoretical overlap $\eta_{ij}=\int \min\{\pi_i f_i,\pi_j f_j\}$ still applies for each pair, and the overall error can be measured by how many true labels are mis-assigned.  But in the latent-mixture setting (as in mixsim simulation), one typically cannot compute $\eta_{\rm global}$ in closed form, so one resorts to the pairwise (`overlap()`) approach.

**R functions:** In R, mixsim provides these tools.  `overlap(Pi,Mu,S)` returns the pairwise OmegaMap and the average (`BarOmega`) and maximum (`MaxOmega`) overlaps.  A related function `overlapGOM(Pi,Mu,S)` computes the “generalized overlap” measure of Maitra (often denoted *goMega*) for $K$ components.  (MixSim also provides `MixSim()` for sampling mixture parameters at fixed Bar/MaxOmega, as well as `simdataset()` to generate data.)

## Comparing Clustering Solutions: Matching and Metrics

When comparing a “true” clustering or parameters to an inferred Gaussian mixture, one must account for label-switching (component permutations).  The clusters are only identified *up to permutation*, so one must optimally match the inferred labels to the true labels before comparing.  A common approach is to solve an assignment problem (Hungarian algorithm) between the two label sets.  For example, one can compute a cost between inferred and true components (e.g.\ Euclidean distance between means, or divergence between distributions) and use the Hungarian (Kuhn–Munkres) algorithm to find the best one-to-one matching that minimizes total cost.  In R, this can be done with the `clue::solve_LSAP` function (from the **clue** package) or similar routines.  This optimal pairing ensures that each inferred cluster is matched to the true cluster it most closely resembles (in location or distribution) – an idea reminiscent of optimal transport on discrete labels.

After aligning labels, one can compare parameters or assignments.  A widely used measure is the **adjusted Rand index** (ARI), which quantifies agreement between two labelings.  In mixsim this is available via `RandIndex(id_true, id_est)`, which returns the Rand index and ARI.  Another is the **Variation of Information** (VI), an information-theoretic distance between partitions.  Mixsim provides `VarInf(id1,id2)` for VI.  Both ARI and VI are *label-invariant* (they give the same value under any permutation of labels), so they avoid having to choose a matching at all.  In practice, one might report ARI/VI to assess clustering accuracy regardless of label permutation.

In summary, best practice is to first align components (e.g. by Hungarian matching on means or overlap) and then compare either the numeric parameters (means, covariances) or the assignments.  Metrics like ARI or VI (implemented by **MixSim** and many other packages) provide an overall similarity score between true and estimated clusterings.  This avoids lexicographic ordering hacks.  If one wants to compare parameters (means, covariances) directly, matching by minimizing mean‐distance or Kullback–Leibler divergence and then computing parameter errors is common.  Recent work also explores using optimal-transport distances (Wasserstein) between component distributions for matching, but the Hungarian approach remains a simple standard.

**Sources:** MixSim package documentation and literature on mixture overlap; theoretical treatment of Gaussian overlap and Bayes error; clustering-comparison metrics.
