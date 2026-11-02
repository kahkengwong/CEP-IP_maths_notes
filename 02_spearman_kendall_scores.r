# 02_spearman_kendall_scores.r
# Section 2. From Spearman's rho and Kendall's tau to the per-cell scores Si
# and Ki, and the pitcher plot.
# Both statistics are sums over cells, so each cell can be given its own
# score (preprint Methods 2.3):
#   Si = 4 (Rx_i - Rbar)(Ry_i - Rbar) / (n - 1)^2
#   Ki = (C_i - D_i) / (n - 1)
# Checks, without ties:
#   rho = 3 x mean(Si) x (n - 1) / (n + 1)   (close to 3 x mean Si)
#   tau_b = mean(Ki)
# Cells with Si > 0.05 and Ki > 0.5 form the SKCS-core.
# The scores, rho, tau and the SKCS-core come from cepip_skcs(). The same
# scores are also computed by hand (cell_scores() in 00_setup.r) and checked
# against the package. Figure 02e is plot_pitcher() from cepip.
#
# Output (written to OUT_DIR):
#   02a_rho_as_pearson_of_ranks.png  v against u with the least-squares line
#   02b_Si_in_copula_space.png       copula space coloured by Si (turbo, 5th to 95th percentile)
#   02c_Ki_in_copula_space.png       copula space coloured by Ki (turbo, 5th to 95th percentile)
#   02d_Si_Ki_distributions.png      histograms of Si and Ki with the gates
#   02e_pitcher_plot.png             Ki against Si with the SKCS-core (cepip plot_pitcher)
#   02_cell_scores.csv               Si, Ki and SKCS-core label per cell
#   02_console_log.txt               console output of this run
#   02_guardrails.txt                PASS/WARN/FAIL ledger
#
# Run: setwd() to the folder that holds 00_setup.r, then source this file.

source("00_setup.r")
start_log("02")

dat <- simulate_cells()
n <- nrow(dat)
dat$u <- to_u(dat$x)
dat$v <- to_u(dat$y)
add_check("G1_no_ties", if (!anyDuplicated(dat$x) && !anyDuplicated(dat$y)) "PASS" else "FAIL",
          "The identities below are exact only without ties")

# ---------------------------------------------------------------------------
# 1. Spearman's rho and Kendall's tau
# ---------------------------------------------------------------------------
# cepip: per-cell scores and the SKCS-core (S_i here is cepip's Si_norm)
sk <- cepip_skcs(dat$x, dat$y, donor = dat$donor, id = dat$cell, fit = FALSE)
stopifnot(identical(sk$cells$id, dat$cell))
dat$Si <- sk$cells$Si_norm
dat$Ki <- sk$cells$Ki
rho  <- sk$rho
tau  <- sk$tau

# The same scores by hand (00_setup.r, cell_scores())
hand <- cell_scores(dat$x, dat$y)
d_si <- max(abs(hand$Si - dat$Si)); d_ki <- max(abs(hand$Ki - dat$Ki))
add_check("G1b_hand_equals_cepip", if (max(d_si, d_ki) < 1e-12) "PASS" else "FAIL",
          sprintf("Si and Ki by hand against cepip_skcs(): max differences %.1e and %.1e", d_si, d_ki))
d_cor <- max(abs(rho - cor(dat$x, dat$y, method = "spearman")), abs(tau - cor(dat$x, dat$y, method = "kendall")))
add_check("G1c_rho_tau", if (d_cor < 1e-12) "PASS" else "FAIL",
          sprintf("cepip rho and tau against cor(): max difference %.1e", d_cor))
rho_from_Si <- 3 * mean(dat$Si) * (n - 1) / (n + 1)
tau_from_Ki <- mean(dat$Ki)

cat(sprintf("\nSpearman rho                    = %.10f\n", rho))
cat(sprintf("3 x mean(Si) x (n-1)/(n+1)      = %.10f\n", rho_from_Si))
cat(sprintf("3 x mean(Si)                    = %.10f (differs by (n+1)/(n-1))\n", 3 * mean(dat$Si)))
cat(sprintf("Kendall tau_b                   = %.10f\n", tau))
cat(sprintf("mean(Ki)                        = %.10f\n\n", tau_from_Ki))

add_check("G2_rho_from_Si", if (abs(rho - rho_from_Si) < 1e-10) "PASS" else "FAIL",
          sprintf("|rho - 3 mean(Si)(n-1)/(n+1)| = %.2e", abs(rho - rho_from_Si)))
add_check("G3_tau_from_Ki", if (abs(tau - tau_from_Ki) < 1e-10) "PASS" else "FAIL",
          sprintf("|tau_b - mean(Ki)| = %.2e", abs(tau - tau_from_Ki)))
add_check("G4_score_ranges",
          if (all(abs(dat$Si) <= 1 + 1e-12) && all(abs(dat$Ki) <= 1 + 1e-12)) "PASS" else "FAIL",
          sprintf("Si in [%.3f, %.3f], Ki in [%.3f, %.3f]",
                  min(dat$Si), max(dat$Si), min(dat$Ki), max(dat$Ki)))

# Si from the copula coordinates (preprint Figure 4): 4(u - 0.5)(v - 0.5)
# equals Si x ((n - 1)/(n + 1))^2, because the ranks are divided by n + 1.
Si_cop <- 4 * (dat$u - 0.5) * (dat$v - 0.5)
fac <- ((n - 1) / (n + 1))^2
add_check("G5_Si_copula_form", if (max(abs(Si_cop - dat$Si * fac)) < 1e-12) "PASS" else "FAIL",
          sprintf("4(u-0.5)(v-0.5) = Si x %.4f, max difference %.2e",
                  fac, max(abs(Si_cop - dat$Si * fac))))

dat$core <- sk$cells$is_core
n_core <- sk$n_core
add_check("G5b_core_gates", if (identical(dat$core, dat$Si > SI_GATE & dat$Ki > KI_GATE)) "PASS" else "FAIL",
          "cepip SKCS-core equals Si > 0.05 and Ki > 0.5")
cat(sprintf("SKCS-core (Si > %.2f and Ki > %.1f): %d of %d cells (%.1f%%)\n",
            SI_GATE, KI_GATE, n_core, n, 100 * n_core / n))
cat("SKCS-core share by simulated source:\n")
print(round(tapply(dat$core, dat$source, mean), 3))
add_check("G6_core_size", if (n_core >= 10) "PASS" else "FAIL", sprintf("SKCS-core n = %d", n_core))

write.csv(dat[, c("cell", "x", "y", "u", "v", "Si", "Ki", "core")],
          file.path(OUT_DIR, "02_cell_scores.csv"), row.names = FALSE)
cat("Saved 02_cell_scores.csv\n")

# ---------------------------------------------------------------------------
# 2. Figures
# ---------------------------------------------------------------------------
p02a <- ggplot(dat, aes(u, v)) +
  geom_point(colour = COL_POINT, alpha = 0.35, size = 0.9) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = COL_LINE_RED, linewidth = 0.9) +
  annotate("text", x = 0.03, y = 0.97, label = sprintf("r(u, v) = %.4f", cor(dat$u, dat$v)),
           hjust = 0, vjust = 1, size = TXT(3.2), fontface = "bold") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = tx("Spearman's \\rho is Pearson's r in copula space", bold = TRUE),
       subtitle = "Least-squares line of v on u", x = "u", y = "v") +
  theme_copula()
save_plot(p02a, "02a_rho_as_pearson_of_ranks.png", w = 5.5, h = 5.5)

score_map <- function(col, title, subtitle, name) {
  ggplot(dat[order(abs(dat[[col]])), ], aes(u, v, colour = .data[[col]])) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_point(size = 0.9) +
    scale_turbo_q(dat[[col]], tx(name, bold = TRUE)) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title, subtitle = subtitle, x = "u", y = "v",
         caption = "Colour scale: 5th to 95th percentile; values outside take the end colours") +
    theme_copula()
}
p02b <- score_map("Si", tx("Per-cell Spearman score S_{i}", bold = TRUE),
                  tx("S_{i} > 0 in the lower-left and upper-right quadrants. Near 0 at the middle GOI rank"),
                  "S_{i}")
save_plot(p02b, "02b_Si_in_copula_space.png", w = 6.2, h = 5.5)
p02c <- score_map("Ki", tx("Per-cell Kendall score K_{i}", bold = TRUE),
                  tx("K_{i} = (concordant - discordant pairs of the cell) / (n - 1)"), "K_{i}")
save_plot(p02c, "02c_Ki_in_copula_space.png", w = 6.2, h = 5.5)

# Facet strips are plotmath strings (label_parsed) so that S_i and K_i carry subscripts
pan <- c(tp("S_{i}", bold = TRUE), tp("K_{i}", bold = TRUE))
hs <- hist_frames(dat$Si, seq(-1, 1, by = 0.05), lo = -1, hi = 1, panel = pan[1])
hk <- hist_frames(dat$Ki, seq(-1, 1, by = 0.05), lo = -1, hi = 1, panel = pan[2])
bars <- rbind(hs$bars, hk$bars); dens <- rbind(hs$dens, hk$dens)
bars$panel <- factor(bars$panel, pan); dens$panel <- factor(dens$panel, pan)
gates <- data.frame(panel = factor(pan, pan), g = c(SI_GATE, KI_GATE),
                    lab = c(tp(sprintf("S_{i} = %.2f", SI_GATE), bold = TRUE),
                            tp(sprintf("K_{i} = %.1f", KI_GATE), bold = TRUE)))
p02d <- ggplot() +
  hist_layers(bars, dens) +
  geom_vline(data = gates, aes(xintercept = g), colour = COL_GATE, linetype = "dashed", linewidth = 0.6) +
  geom_text(data = gates, aes(x = g + 0.04, y = Inf, label = lab), vjust = 1.5, hjust = 0, parse = TRUE,
            size = TXT(3.0), fontface = "bold", colour = COL_GATE) +
  facet_wrap(~ panel, labeller = label_parsed) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Distributions of the per-cell scores",
       subtitle = "Bin width 0.05. Dashed lines: the SKCS-core gates",
       x = "Score", y = "Cells (%)") +
  theme_copula()
save_plot(p02d, "02d_Si_Ki_distributions.png", w = 8, h = 4)

# Figure 02e: plot_pitcher() from cepip, with subscripted labels
p02e <- restyle_cepip(
  plot_pitcher(sk),
  text_map = list("Si > 0.05" = tp("S_{i} > 0.05", bold = TRUE),
                  "Ki > 0.50" = tp("K_{i} > 0.50", bold = TRUE)),
  subtitle = tx(sprintf("\\rho = %.3f = 3 x mean(S_{i}) x (n-1)/(n+1); \\tau_{b} = %.3f = mean(K_{i})", rho, tau)),
  x = tx("Spearman rank-product contribution (S_{i})", bold = TRUE),
  y = tx("Kendall concordance contribution (K_{i})", bold = TRUE))
save_plot(p02e, "02e_pitcher_plot.png", w = 6, h = 6)

finish_log("02_spearman_kendall_scores.r")
