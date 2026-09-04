#!/usr/bin/env bash
# Cloud Agent setup for the DeCovarT R package.
#
# Self-contained and idempotent: it can bootstrap a bare Ubuntu image
# (installing R, the r2u binary-CRAN apt repository, Pandoc, Quarto and
# pre-commit) and, on a warm base where those already exist, becomes a
# fast dependency refresh. Everything is keyed off "is it already there?"
# so it is safe to run repeatedly.
#
# Layers:
#   1. System toolchain  - R + r2u/bspm (binary CRAN), Pandoc, Quarto,
#                          pre-commit.
#   2. R package deps     - Imports, Suggests, dev/check tooling, the
#                          optional Bioconductor package ComplexHeatmap,
#                          and pkgcheck (rOpenSci), installed as apt
#                          binaries via bspm where possible.
#   3. Repository state   - install the package so library(DeCovarT)
#                          works and pre-build the pre-commit hooks.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

export DEBIAN_FRONTEND=noninteractive
export PATH="$HOME/.local/bin:$PATH"
NCPUS="$(nproc 2>/dev/null || echo 2)"

# ---------------------------------------------------------------------------
# 1. System toolchain: R + r2u (binary CRAN via bspm)
# ---------------------------------------------------------------------------
if ! command -v R >/dev/null 2>&1; then
  echo "==> Installing R + r2u (binary CRAN) for Ubuntu"
  . /etc/os-release
  CODENAME="${VERSION_CODENAME:-noble}"
  case "$CODENAME" in
    noble|jammy) : ;;
    *)
      echo "WARN: untested Ubuntu codename '$CODENAME'; defaulting to noble." >&2
      CODENAME="noble"
      ;;
  esac

  sudo apt-get update -qq
  sudo apt-get install -y --no-install-recommends \
    ca-certificates gnupg wget curl locales make

  sudo locale-gen en_GB.UTF-8 en_US.UTF-8 || true

  # CRAN R apt repo (marutter) and the r2u binary-CRAN repo, each with its
  # own key file (a shared file would let one key clobber the other).
  wget -q -O- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
    | sudo tee /etc/apt/trusted.gpg.d/cran_r_key.asc >/dev/null
  echo "deb [arch=amd64] https://cloud.r-project.org/bin/linux/ubuntu ${CODENAME}-cran40/" \
    | sudo tee /etc/apt/sources.list.d/cran_r.list >/dev/null

  wget -q -O- https://eddelbuettel.github.io/r2u/assets/dirk_eddelbuettel_key.asc \
    | sudo tee /etc/apt/trusted.gpg.d/r2u_key.asc >/dev/null
  echo "deb [arch=amd64] https://r2u.stat.illinois.edu/ubuntu ${CODENAME} main" \
    | sudo tee /etc/apt/sources.list.d/cranapt.list >/dev/null

  sudo tee /etc/apt/preferences.d/99cranapt >/dev/null <<'EOF'
Package: *
Pin: release o=CRAN-Apt Project
Pin: release l=CRAN-Apt Packages
Pin-Priority: 700
EOF

  sudo apt-get update -qq
  sudo apt-get install -y --no-install-recommends \
    r-base-core r-base-dev python3-dbus python3-gi python3-apt r-cran-bspm

  # Configure bspm: set options BEFORE enable() so the D-Bus backend is
  # skipped (sudo mode) and no warning is emitted for non-root users.
  RPROFILE_SITE="$(R RHOME)/etc/Rprofile.site"
  if ! sudo grep -q 'bspm::enable' "$RPROFILE_SITE" 2>/dev/null; then
    sudo tee -a "$RPROFILE_SITE" >/dev/null <<'EOF'

## r2u / bspm: install CRAN (and Bioconductor r-bioc-*) packages as apt
## binaries. Options must precede enable() so non-root users use sudo mode.
options(bspm.version.check = FALSE)
options(bspm.sudo = TRUE)
suppressMessages(bspm::enable())
EOF
  fi
  # Make the local site-library writable so source-installs / user installs
  # (e.g. pkgcheck from GitHub) do not need root.
  sudo mkdir -p /usr/local/lib/R/site-library
  sudo chown -R "$(id -u):$(id -g)" /usr/local/lib/R/site-library
else
  echo "==> R already present: $(R --version | head -1)"
fi

# ---------------------------------------------------------------------------
# Pandoc, Quarto, pre-commit
# ---------------------------------------------------------------------------
if ! command -v pandoc >/dev/null 2>&1; then
  echo "==> Installing Pandoc"
  sudo apt-get install -y --no-install-recommends pandoc
fi

if ! command -v quarto >/dev/null 2>&1; then
  echo "==> Installing Quarto CLI"
  QVER="$(curl -sL https://quarto.org/docs/download/_download.json \
    | grep -oP '"version":\s*"\K[0-9.]+' | head -1)"
  curl -sL -o /tmp/quarto.deb \
    "https://github.com/quarto-dev/quarto-cli/releases/download/v${QVER}/quarto-${QVER}-linux-amd64.deb"
  sudo dpkg -i /tmp/quarto.deb
  rm -f /tmp/quarto.deb
fi

if ! command -v pre-commit >/dev/null 2>&1; then
  echo "==> Installing pre-commit (pipx)"
  sudo apt-get install -y --no-install-recommends pipx
  pipx install pre-commit
  pipx ensurepath >/dev/null 2>&1 || true
fi
# Ensure interactive agent shells find pipx tools.
if ! grep -q 'HOME/.local/bin' "$HOME/.bashrc" 2>/dev/null; then
  printf '\n# Expose pipx-installed CLI tools (pre-commit)\nexport PATH="$HOME/.local/bin:$PATH"\n' \
    >> "$HOME/.bashrc"
fi

# ---------------------------------------------------------------------------
# 2. R package dependencies (binary via r2u/bspm)
# ---------------------------------------------------------------------------
echo "==> Installing R dependencies + dev/check tooling"
NCPUS="$NCPUS" Rscript --no-init-file - <<'RS'
options(bspm.sudo = TRUE, warn = 1L,
        Ncpus = as.integer(Sys.getenv("NCPUS", "2")))

dev_tools <- c(
  "devtools", "roxygen2", "rcmdcheck", "testthat", "covr", "lintr",
  "desc", "docopt", "cli", "remotes", "pkgdown", "visNetwork"
)
to_install <- setdiff(dev_tools, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

# DESCRIPTION Imports + Suggests (CRAN resolved to apt binaries by bspm).
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_deps(".", dependencies = TRUE, upgrade = "never")

# Optional Bioconductor visualisation dependency (r2u ships r-bioc-*).
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  try(install.packages("ComplexHeatmap"), silent = TRUE)
}

# rOpenSci pkgcheck: not on CRAN for current R; install from GitHub
# (best-effort, non-fatal).
if (!requireNamespace("pkgcheck", quietly = TRUE)) {
  try(
    remotes::install_github(
      "ropensci-review-tools/pkgcheck", upgrade = "never"
    ),
    silent = TRUE
  )
}
RS

# ---------------------------------------------------------------------------
# 3. Repository state: install the package + pre-commit hook environments
# ---------------------------------------------------------------------------
echo "==> Installing the DeCovarT package (pure R, no compilation)"
R CMD INSTALL --no-multiarch --with-keep.source .

echo "==> Pre-building pre-commit hook environments"
# install-hooks builds the hook environments only; `pre-commit install`
# would try to wire .git/hooks and fails when core.hooksPath is set.
if command -v pre-commit >/dev/null 2>&1; then
  pre-commit install-hooks || \
    echo "WARN: pre-commit hook environment build failed (non-fatal)."
fi

echo "==> DeCovarT environment setup complete"
