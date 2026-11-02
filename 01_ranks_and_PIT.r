# 01_ranks_and_PIT.r
# Section 1. Ranks and the probability integral transform (PIT).
# Maps the simulated GOI-positive cells from expression space (GOI, composite)
# to copula space (u, v), shows that the copula coordinates do not change
# under strictly increasing transformations, and shows why cells with a zero
# GOI count are left out before the IP is located.
# The composite y is the mean of eight partner genes, computed with
# cepip_composite() in 00_setup.r.
#
# Output (written to OUT_DIR):
#   01a_expression_space.png    composite against GOI
#   01b_expression_margins.png  histograms of GOI and composite with density curves
#   01c_copula_space.png        v against u
#   01d_invariance.png          monotone transforms change expression space, not copula space
#   01e_ties_zero_GOI.png       cells with GOI = 0 share one u value
#   01_console_log.txt          console output of this run
#   01_guardrails.txt           PASS/WARN/FAIL ledger
#
# Run: setwd() to the folder that holds 00_setup.r, then source this file.

source("00_setup.r")
start_log("01")

dat <- simulate_cells()
N <- nrow(dat)
add_check("G1_complete", if (all(is.finite(c(dat$x, dat$y)))) "PASS" else "FAIL",
          sprintf("n = %d GOI-positive cells from %d donors, non-finite values = %d",
                  N, length(unique(dat$donor)), sum(!is.finite(c(dat$x, dat$y)))))
comp_diff <- max(abs(dat$y - rowMeans(dat[, paste0("P", 1:8)])))
add_check("G1b_composite", if (comp_diff < 1e-12) "PASS" else "FAIL",
          sprintf("cepip_composite() equals the row mean of P1 to P8, max difference %.1e", comp_diff))

# ---------------------------------------------------------------------------
# 1. Probability integral transform by ranks
# u_i = rank(x_i) / (n + 1), v_i = rank(y_i) / (n + 1).
# Dividing by n + 1 keeps every value strictly inside (0, 1).
# ---------------------------------------------------------------------------
dat$u <- to_u(dat$x)
dat$v <- to_u(dat$y)
add_check("G2_open_interval",
          if (all(dat$u > 0 & dat$u < 1 & dat$v > 0 & dat$v < 1)) "PASS" else "FAIL",
          sprintf("u in [%.4f, %.4f], v in [%.4f, %.4f]",
                  min(dat$u), max(dat$u), min(dat$v), max(dat$v)))
ks_u <- suppressWarnings(ks.test(dat$u, "punif"))$p.value
ks_v <- suppressWarnings(ks.test(dat$v, "punif"))$p.value
add_check("G3_uniform_margins", if (min(ks_u, ks_v) > 0.05) "PASS" else "WARN",
          sprintf("KS test against Uniform(0,1): p(u) = %.3f, p(v) = %.3f", ks_u, ks_v))

cat(sprintf("\nPearson r in expression space = %.4f\n", cor(dat$x, dat$y)))
cat(sprintf("Pearson r in copula space     = %.4f (this is Spearman's rho)\n", cor(dat$u, dat$v)))
cat(sprintf("Spearman rho from cor()       = %.4f\n\n", cor(dat$x, dat$y, method = "spearman")))

# ---------------------------------------------------------------------------
# 2. Figures 01a to 01c
# ---------------------------------------------------------------------------
p01a <- ggplot(dat, aes(x, y)) +
  geom_point(colour = COL_POINT, alpha = 0.35, size = 0.9) +
  labs(title = "Expression space",
       subtitle = sprintf("n = %d GOI-positive cells", N),
       x = "GOI (log1p)", y = "Composite (log1p)") +
  theme_copula()
save_plot(p01a, "01a_expression_space.png")

brk_x <- seq(0, ceiling(max(dat$x) * 2) / 2, by = 0.25)
brk_y <- seq(0, ceiling(max(dat$y) * 2) / 2, by = 0.25)
hx <- hist_frames(dat$x, brk_x, lo = 0, panel = "GOI")
hy <- hist_frames(dat$y, brk_y, lo = 0, panel = "Composite")
bars <- rbind(hx$bars, hy$bars); dens <- rbind(hx$dens, hy$dens)
bars$panel <- factor(bars$panel, c("GOI", "Composite"))
dens$panel <- factor(dens$panel, c("GOI", "Composite"))
p01b <- ggplot() +
  hist_layers(bars, dens) +
  facet_wrap(~ panel, scales = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Margins in expression space",
       subtitle = "Bin width 0.25. Curve: kernel density reflected at 0",
       x = "Expression (log1p)", y = "Cells (%)") +
  theme_copula()
save_plot(p01b, "01b_expression_margins.png", w = 8, h = 4)

p01c <- ggplot(dat, aes(u, v)) +
  geom_point(colour = COL_POINT, alpha = 0.35, size = 0.9) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Copula space",
       subtitle = "u = rank(GOI)/(n+1), v = rank(composite)/(n+1)",
       x = "u", y = "v") +
  theme_copula()
save_plot(p01c, "01c_copula_space.png", w = 5.5, h = 5.5)

# ---------------------------------------------------------------------------
# 3. Invariance under strictly increasing transformations (figure 01d)
# expm1() undoes the log1p step and returns both axes to the CP10K scale.
# The ranks, and therefore (u, v), do not change.
# ---------------------------------------------------------------------------
dat$x_t <- expm1(dat$x)
dat$y_t <- expm1(dat$y)
dat$u_t <- to_u(dat$x_t)
dat$v_t <- to_u(dat$y_t)
max_diff <- max(abs(dat$u - dat$u_t), abs(dat$v - dat$v_t))
add_check("G4_invariance", if (max_diff == 0) "PASS" else "FAIL",
          sprintf("max |u - u_t|, |v - v_t| = %g", max_diff))

inv <- rbind(
  data.frame(panel = "1. Expression space, log1p scale",  a = dat$x,   b = dat$y),
  data.frame(panel = "2. Expression space, CP10K scale",  a = dat$x_t, b = dat$y_t),
  data.frame(panel = "3. Copula space from log1p values", a = dat$u,   b = dat$v),
  data.frame(panel = "4. Copula space from CP10K values", a = dat$u_t, b = dat$v_t))
p01d <- ggplot(inv, aes(a, b)) +
  geom_point(colour = COL_POINT, alpha = 0.30, size = 0.6) +
  facet_wrap(~ panel, scales = "free") +
  labs(title = "Strictly increasing transformations do not change copula space",
       subtitle = "Panels 3 and 4 are identical point for point",
       x = NULL, y = NULL) +
  theme_copula()
save_plot(p01d, "01d_invariance.png", w = 8, h = 7)

# ---------------------------------------------------------------------------
# 4. Ties from zero GOI counts (figure 01e)
# In single-cell data many cells have a zero count for the GOI. All tied zeros
# receive the same average rank, so they share one u value and form a
# vertical stripe. The stripe has no order along u, which is why the IP is
# located on GOI-positive cells only.
# ---------------------------------------------------------------------------
set.seed(SEED + 1)
zero_frac <- 0.30
x0 <- dat$x
x0[sample(N, round(zero_frac * N))] <- 0
n0 <- sum(x0 == 0)
u0 <- to_u(x0)
u_stripe_expected <- ((n0 + 1) / 2) / (N + 1)
u_stripe_observed <- unique(u0[x0 == 0])
add_check("G5_tie_stripe",
          if (length(u_stripe_observed) == 1 &&
              abs(u_stripe_observed - u_stripe_expected) < 1e-12) "PASS" else "FAIL",
          sprintf("%d tied zeros; stripe at u = %.4f (expected (n0+1)/2/(n+1) = %.4f)",
                  n0, u_stripe_observed[1], u_stripe_expected))

lab_zero <- sprintf("GOI = 0, tied (n=%d)", n0)
lab_pos  <- sprintf("GOI > 0 (n=%d)", N - n0)
ties_df <- data.frame(u = u0, v = dat$v,
                      zero = factor(ifelse(x0 == 0, lab_zero, lab_pos), levels = c(lab_pos, lab_zero)))
ties_df <- ties_df[order(ties_df$zero), ]
p01e <- ggplot(ties_df, aes(u, v, colour = zero)) +
  geom_vline(xintercept = u_stripe_expected, colour = "black", linetype = "dashed", linewidth = 0.4) +
  geom_point(alpha = 0.5, size = 0.9) +
  annotate("text", x = u_stripe_expected + 0.02, y = 0.03,
           label = sprintf("u = %.3f", u_stripe_expected), hjust = 0, size = TXT(3.0), fontface = "bold") +
  scale_colour_manual(values = setNames(c(COL_OTHER, COL_TIE), c(lab_pos, lab_zero))) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Tied zeros share one average rank",
       subtitle = tx(sprintf("%d of %d cells set to GOI = 0. Stripe at u = (n_{0}+1)/2/(n+1)", n0, N)),
       x = "u", y = "v", colour = NULL) +
  theme_copula() +
  legend_inside(0.98, 0.02)
save_plot(p01e, "01e_ties_zero_GOI.png", w = 5.5, h = 5.5)

finish_log("01_ranks_and_PIT.r")
