###############################################################################
# Figure for vignette 03-1, #fig-dea-design: D-, E- and A-optimal choice
# between two gene panels when genes are correlated. Three ggplot2 panels,
# assembled with the notation in vignettes/figures/src/dea_design/dea_design.tex.
#
# Run from the repository root (system R, not renv):
#   Rscript --no-save --no-restore scripts/figures/fig_03_1_dea_design.R
#   cd vignettes/figures/src/dea_design && latexmk -lualatex dea_design.tex
#   pdftocairo -svg dea_design.pdf ../../fig_03_1_dea_design.svg
#
# Model (03-1, #sec-optimal-design; ledger V18): single-cell bulk law with a
# common within-type covariance, so S_W(p) = Sigma_W and the information for
# the orthonormal (Helmert) coordinates z of the simplex is
#   I = n K' M' Sigma_W^{-1} M K,   Cov(z-hat) = I^{-1}.
# Example of ledger V18d (J = 3 types, G = 3 genes per panel, unit variances,
# n = 25 for the scale of the ellipses). Deterministic: no random numbers
# are drawn (the seed is set for completeness).
###############################################################################

library(ggplot2)
library(DeCovarT) # theme_decovart_facets(), the repository theme

set.seed(20261101)

out_dir <- file.path("vignettes", "figures", "src", "dea_design")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

col_a <- "#1d4e89"
col_b <- "#c45c26"
col_guide <- "#1c2b24"

n_i <- 25
K <- cbind(c(1, -1, 0) / sqrt(2), c(1, 1, -2) / sqrt(6)) # Helmert, J x (J - 1)
M_A <- rbind(c(4, 2, 2), c(4, 2, 0), c(2, 1, 4)) # genes x types
M_B <- rbind(c(2, 4, 2), c(0, 2, 4), c(0, 2, 3))
rho <- 0.8
S_rho <- diag(3)
S_rho[1, 3] <- S_rho[3, 1] <- rho

info <- function(M, S) n_i * t(K) %*% t(M) %*% solve(S) %*% M %*% K

panels <- list(
  "A, uncorrelated" = info(M_A, diag(3)),
  "B, uncorrelated" = info(M_B, diag(3)),
  "A, rho = 0.8" = info(M_A, S_rho)
)

criteria <- do.call(
  rbind,
  lapply(names(panels), function(nm) {
    I <- panels[[nm]]
    ev <- eigen(I, symmetric = TRUE)$values
    data.frame(
      panel = nm,
      D = det(I),
      E = min(ev),
      A = sum(diag(solve(I))),
      kappa = max(ev) / min(ev)
    )
  })
)
print(criteria) # D, E, A divided by n_i^2, n_i, 1/n_i reproduce ledger V18d

# ---- ellipses on the simplex plane -------------------------------------------
vertices <- t(K) %*% (diag(3) - 1 / 3) # columns: types 1..3 in z-coordinates
tri <- data.frame(
  z1 = vertices[1, ],
  z2 = vertices[2, ],
  type = paste("type", 1:3)
)
ellipse <- function(I, level = 0.95, n = 200) {
  S <- solve(I)
  r <- sqrt(stats::qchisq(level, df = 2))
  ang <- seq(0, 2 * pi, length.out = n)
  e <- eigen(S, symmetric = TRUE)
  pts <- e$vectors %*% diag(sqrt(e$values)) %*% rbind(cos(ang), sin(ang)) * r
  data.frame(z1 = pts[1, ], z2 = pts[2, ])
}
ell_df <- function(nms, scenario) {
  do.call(
    rbind,
    lapply(nms, function(nm) {
      d <- ellipse(panels[[nm]])
      d$panel <- nm
      d$scenario <- scenario
      d
    })
  )
}

theme_fig <- theme_decovart_facets(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(4, 4, 4, 4)
  )

simplex_plot <- function(d, colours, linetypes) {
  ggplot() +
    geom_polygon(
      data = tri,
      aes(z1, z2),
      fill = NA,
      colour = "grey55",
      linewidth = 0.4
    ) +
    geom_text(
      data = tri,
      aes(z1 * 1.12, z2 * 1.12, label = type),
      size = 3,
      colour = col_guide
    ) +
    geom_point(aes(0, 0), colour = col_guide, size = 0.8) +
    geom_path(
      data = d,
      aes(z1, z2, colour = panel, linetype = panel),
      linewidth = 0.8
    ) +
    scale_colour_manual(values = colours) +
    scale_linetype_manual(values = linetypes) +
    coord_equal(xlim = c(-0.9, 0.9), ylim = c(-0.95, 0.6)) +
    theme_fig +
    guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
}

save_panel <- function(p, name, width = 2.9, height = 3.0) {
  ggsave(
    filename = file.path(out_dir, paste0(name, ".pdf")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )
}

pa <- simplex_plot(
  ell_df(c("A, uncorrelated", "B, uncorrelated"), "independent genes"),
  colours = c(col_a, col_b),
  linetypes = c("solid", "solid")
)
save_panel(pa, "panel_a_independent")

pb <- simplex_plot(
  ell_df(c("A, uncorrelated", "A, rho = 0.8", "B, uncorrelated"), "correlated"),
  colours = c(
    "A, uncorrelated" = "grey60",
    "A, rho = 0.8" = col_a,
    "B, uncorrelated" = col_b
  ),
  linetypes = c(
    "A, uncorrelated" = "22",
    "A, rho = 0.8" = "solid",
    "B, uncorrelated" = "solid"
  )
)
save_panel(pb, "panel_b_correlated")

# ---- which panel each criterion prefers ---------------------------------------
adv <- function(a, b) {
  data.frame(
    criterion = factor(
      c("D: log det I", "E: log lambda_min(I)", "A: -log tr I^-1"),
      levels = c("A: -log tr I^-1", "E: log lambda_min(I)", "D: log det I")
    ),
    advantage = c(log(a$D / b$D), log(a$E / b$E), log(b$A / a$A))
  )
}
cr <- split(criteria, criteria$panel)
df_c <- rbind(
  transform(
    adv(cr[["A, uncorrelated"]], cr[["B, uncorrelated"]]),
    scenario = "genes independent in A"
  ),
  transform(
    adv(cr[["A, rho = 0.8"]], cr[["B, uncorrelated"]]),
    scenario = "genes 1 and 3 of A correlated"
  )
)
df_c$scenario <- factor(
  df_c$scenario,
  levels = c("genes independent in A", "genes 1 and 3 of A correlated")
)
pc <- ggplot(
  df_c,
  aes(advantage, criterion, colour = scenario, shape = scenario)
) +
  geom_vline(xintercept = 0, colour = col_guide, linewidth = 0.4) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c("grey45", col_a)) +
  scale_shape_manual(values = c(1, 16)) +
  scale_y_discrete(labels = c("A", "E", "D")) +
  labs(x = "log advantage of panel A over B", y = "criterion") +
  theme_decovart_facets(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.margin = margin(4, 6, 4, 4)
  ) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
save_panel(pc, "panel_c_criteria", width = 3.0)

utils::sessionInfo()
