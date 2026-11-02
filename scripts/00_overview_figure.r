# 00_overview_figure.r
# Overview figure: one panel for each type of plot in scripts 01 to 04, in
# the order the scripts produce them. Each panel names the figure it comes
# from. Runs on its own from the same simulated cells (00_setup.r), so it
# can be run before or after the other scripts. The numbers come from the
# cepip package; the panels are drawn here with ggplot2.
#
# Panels:
#   A  Expression space                 (01a)
#   B  Histogram with density curve     (01b)
#   C  Copula space                     (01c)
#   D  Per-cell score map               (02b)
#   E  GAM plot with control composite  (03a)
#   F  Derivatives of the fitted curve  (03d)
#   G  Residual curve R(x)              (03e)
#   H  Brent roots of the 95% CI        (03f)
#   I  IPRS components                  (03h)
#   J  CEP-IP plot by group             (04a)
#   K  Pitcher plot, score space        (04b)
#   L  Pitcher plot, copula space       (04e)
#
# Output (figures/, tables/ and logs/ of the repository folder):
#   00_overview.png               the 12-panel overview
#   00_console_log.txt            console output of this run
#   00_guardrails.txt             PASS/WARN/FAIL ledger
#
# Run: setwd() to the repository folder, then source this file from scripts/.

source("scripts/00_setup.r")
start_log("00")

BASE <- 9
# Panel titles and subtitles accept the markup of tp() (for example S_{i}).
tag_title <- function(tag, title, src) labs(title = tx(sprintf("%s. %s", tag, title), bold = TRUE),
                                            subtitle = tx(src))

# ---------------------------------------------------------------------------
# 1. Data and pipeline (cepip, as in scripts 01 to 04)
# ---------------------------------------------------------------------------
dat <- simulate_cells()
n <- nrow(dat)
dat$u <- to_u(dat$x); dat$v <- to_u(dat$y)

f_all  <- cepip_fit(dat$x, dat$y, donor = dat$donor, id = dat$cell)
f_ctrl <- cepip_fit(dat$x, dat$y_ctrl, donor = dat$donor, id = dat$cell)
ip     <- cepip_ip(dat$x, dat$y, donor = dat$donor, id = dat$cell, fit = f_all)
sk     <- cepip_skcs(dat$x, dat$y, donor = dat$donor, id = dat$cell)
dat$Si <- sk$cells$Si_norm; dat$Ki <- sk$cells$Ki; dat$core <- sk$cells$is_core
core_idx <- which(dat$core)
dat$trep_core <- FALSE; dat$trep_core[core_idx] <- sk$ep$is_TREP
dat$group <- "Non-SKCS-core"
side <- ifelse(dat$y < ip$ip_y, "Pre-IP", "Post-IP")
dat$group[dat$core] <- paste(side[dat$core], ifelse(dat$trep_core[dat$core], "TREP", "non-TREP"))
dat$group <- factor(dat$group, names(GROUP_COLS))
dat$u_SK <- to_u(dat$Si); dat$v_SK <- to_u(dat$Ki)
ugate <- sum(dat$Si <= SI_GATE) / (n + 1)
vgate <- sum(dat$Ki <= KI_GATE) / (n + 1)

ctl <- simulate_control_arm()
ipc <- cepip_ip(ctl$x, ctl$y, donor = ctl$donor, id = ctl$cell)

cat(sprintf("n = %d, SKCS-core = %d, k = %d, DE = %.3f, IPx = %.4f, IPy = %.4f, IPRS = %.3f (%s)\n\n",
            n, sum(dat$core), f_all$k, f_all$de, ip$ip_x, ip$ip_y, ip$iprs$IPRS, ip$iprs$tier))
add_check("G1_IP_found", if (isTRUE(ip$has_ip)) "PASS" else "FAIL",
          sprintf("IPx = %.4f, IPRS = %.3f (%s tier)", ip$ip_x, ip$iprs$IPRS, ip$iprs$tier))
add_check("G2_groups_complete", if (sum(table(dat$group)) == n) "PASS" else "FAIL",
          paste(sprintf("%s %d", levels(dat$group), as.integer(table(dat$group))), collapse = "; "))

pd <- dat[order(dat$group == "Non-SKCS-core", decreasing = TRUE), ]
grp_lab <- setNames(sprintf("%s (n=%d)", names(GROUP_COLS), as.integer(table(dat$group))), names(GROUP_COLS))

# ---------------------------------------------------------------------------
# 2. Panels
# ---------------------------------------------------------------------------
P <- list()

P$A <- ggplot(dat, aes(x, y)) +
  geom_point(colour = COL_POINT, alpha = 0.35, size = 0.6) +
  tag_title("A", "Expression space", "Composite against GOI (01a)") +
  labs(x = "GOI (log1p)", y = "Composite (log1p)") +
  theme_copula(BASE)

hx <- hist_frames(dat$x, seq(0, ceiling(max(dat$x) * 2) / 2, by = 0.25), lo = 0)
P$B <- ggplot() +
  hist_layers(hx$bars, hx$dens) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  tag_title("B", "Histogram with density curve", "GOI margin (01b)") +
  labs(x = "GOI (log1p)", y = "Cells (%)") +
  theme_copula(BASE)

P$C <- ggplot(dat, aes(u, v)) +
  geom_point(colour = COL_POINT, alpha = 0.35, size = 0.6) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  tag_title("C", "Copula space", "v against u, ranks divided by n + 1 (01c)") +
  labs(x = "u", y = "v") +
  theme_copula(BASE)

P$D <- ggplot(dat[order(abs(dat$Si)), ], aes(u, v, colour = Si)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_point(size = 0.6) +
  scale_turbo_q(dat$Si, tx("S_{i}", bold = TRUE)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  tag_title("D", "Per-cell score map", "S_{i} in copula space, colours at 5th to 95th percentile (02b)") +
  labs(x = "u", y = "v") +
  theme_copula(BASE) +
  theme(legend.key.height = unit(0.35, "cm"), legend.key.width = unit(0.25, "cm"))

gx <- seq(min(dat$x), max(dat$x), length.out = 300)
cv <- function(fit, what) {
  p <- fit_curve(fit, gx, se = TRUE)
  data.frame(what = what, x = gx, fit = as.numeric(p$fit), se = as.numeric(p$se.fit))
}
crv <- rbind(cv(f_all, "Composite"), cv(f_ctrl, "Control composite"))
pts <- rbind(data.frame(what = "Control composite", x = dat$x, y = dat$y_ctrl),
             data.frame(what = "Composite", x = dat$x, y = dat$y))
what_cols <- c("Composite" = COL_TREP, "Control composite" = COL_CTRL)
P$E <- ggplot() +
  geom_point(data = pts, aes(x, y, colour = what), alpha = 0.2, size = 0.5) +
  geom_ribbon(data = crv, aes(x = x, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se, fill = what), alpha = 0.3) +
  geom_line(data = crv, aes(x = x, y = fit, colour = what), linewidth = 0.9) +
  annotate("text", x = -Inf, y = Inf, label = sprintf("DE=%.1f%%", 100 * f_all$de),
           hjust = -0.1, vjust = 1.4, size = TXT(2.8), fontface = "bold") +
  scale_colour_manual(values = what_cols) + scale_fill_manual(values = what_cols, guide = "none") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1.8, linewidth = 0.9))) +
  tag_title("E", "GAM plot", "Composite and control composite (03a)") +
  labs(x = "GOI (log1p)", y = "Composite (log1p)", colour = NULL) +
  theme_copula(BASE) +
  legend_inside(0.99, 0.01)

g  <- seq(min(dat$x), max(dat$x), length.out = 2000)
h  <- diff(g)[1]
fg <- as.numeric(fit_curve(f_all, g))
f2 <- c(NA, (fg[-(1:2)] - 2 * fg[2:(length(g) - 1)] + fg[1:(length(g) - 2)]) / h^2, NA)
ok <- !is.na(f2[-1]) & !is.na(f2[-length(f2)])
sc2 <- which(diff(sign(f2)) != 0 & ok)
x_f2 <- g[sc2] - f2[sc2] * h / (f2[sc2 + 1] - f2[sc2])
P$F <- ggplot(data.frame(x = g, f2 = f2), aes(x, f2)) +
  geom_hline(yintercept = 0, colour = "gray40", linetype = "dotted") +
  geom_vline(xintercept = x_f2, colour = "grey30", linetype = "dotted", linewidth = 0.5) +
  geom_vline(xintercept = ip$ip_x, colour = COL_GATE, linetype = "dashed", linewidth = 0.7) +
  geom_line(colour = COL_TREP, linewidth = 0.8, na.rm = TRUE) +
  tag_title("F", "Derivatives of the fitted curve", "f''(x) with sign changes; green: IP_{x} (03d)") +
  labs(x = "GOI (log1p)", y = "f''(x)") +
  theme_copula(BASE)

# Mean TREP residual over 100 equal-width GOI bins (bars), from cepip_ip()
rd  <- ip$residual_data
br  <- seq(min(rd$GOI), max(rd$GOI), length.out = 101)
bin <- cut(rd$GOI, br, include.lowest = TRUE, labels = FALSE)
agg <- aggregate(rd$residual, list(bin = bin), mean)
bins <- data.frame(GOI = (br[agg$bin] + br[agg$bin + 1]) / 2, mean_residual = agg$x)
P$G <- ggplot() +
  annotate("rect", xmin = ip$ip_ci[1], xmax = ip$ip_ci[2], ymin = -Inf, ymax = Inf, fill = COL_GATE, alpha = 0.08) +
  geom_col(data = bins, aes(GOI, mean_residual), width = 0.8 * diff(br)[1],
           fill = COL_TREP, colour = "#3D0065", alpha = 0.40, linewidth = 0.2) +
  geom_ribbon(data = ip$grid, aes(x = x, ymin = r_lower, ymax = r_upper), fill = COL_GAM, alpha = 0.35) +
  geom_line(data = ip$grid, aes(x, r_hat), colour = COL_GAM, linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray40", linewidth = 0.4) +
  geom_vline(xintercept = ip$ip_ci, linetype = "dotdash", colour = COL_GATE, linewidth = 0.5) +
  geom_vline(xintercept = ip$ip_x, linetype = "dashed", colour = COL_GATE, linewidth = 0.8) +
  tag_title("G", "Residual curve R(x)", sprintf("IP_{x} = %.3f, IPRS = %.2f (03e)", ip$ip_x, ip$iprs$IPRS)) +
  labs(x = "GOI (log1p)", y = "R(x)") +
  theme_copula(BASE)

# Window of 2.5 CI widths on each side of IPx
win <- ip$ip_x + c(-1, 1) * max(0.15, 2.5 * diff(ip$ip_ci))
zg  <- ip$grid[ip$grid$x >= win[1] & ip$grid$x <= win[2], ]
zl  <- rbind(data.frame(x = zg$x, val = zg$r_hat,   line = "R(x)"),
             data.frame(x = zg$x, val = zg$r_upper, line = "R(x) + 1.96 SE"),
             data.frame(x = zg$x, val = zg$r_lower, line = "R(x) - 1.96 SE"))
roots <- data.frame(x = c(ip$ip_ci[1], ip$ip_x, ip$ip_ci[2]), y = 0,
                    lab = vapply(c("Lower", "IP_{x}", "Upper"), tp, character(1), bold = TRUE),
                    vj = c(1.8, -0.8, 1.8))
P$H <- ggplot(zl, aes(x, val)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray40") +
  geom_line(aes(linetype = line, linewidth = line), colour = "#D9822B") +
  geom_point(data = roots, aes(x, y), size = 1.6, colour = "black") +
  geom_text(data = roots, aes(x, y, label = lab, vjust = vj), size = TXT(2.6), parse = TRUE) +
  scale_linetype_manual(values = c("R(x)" = "solid", "R(x) + 1.96 SE" = "dashed", "R(x) - 1.96 SE" = "dotdash"),
                        guide = "none") +
  scale_linewidth_manual(values = c("R(x)" = 1.0, "R(x) + 1.96 SE" = 0.5, "R(x) - 1.96 SE" = 0.5), guide = "none") +
  tag_title("H", "Brent roots of the 95% CI", "R(x) and its band edges near IP_{x} (03f)") +
  labs(x = "GOI (log1p)", y = "R(x)") +
  theme_copula(BASE)

lab_m <- "Simulated module"; lab_c <- "Simulated control"
cmp <- rbind(data.frame(arm = lab_m, comp = c("C1", "C2", "C3", "C4", "C5", "IPRS"),
                        val = c(ip_components(ip), ip$iprs$IPRS)),
             data.frame(arm = lab_c, comp = c("C1", "C2", "C3", "C4", "C5", "IPRS"),
                        val = c(ip_components(ipc), ipc$iprs$IPRS)))
cmp$arm  <- factor(cmp$arm, c(lab_m, lab_c))
cmp$comp <- factor(cmp$comp, c("C1", "C2", "C3", "C4", "C5", "IPRS"))
P$I <- ggplot(cmp, aes(comp, val, fill = arm)) +
  geom_hline(yintercept = c(0.2, 0.4, 0.6, 0.8), colour = "gray50", linetype = "dotted") +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = setNames(c(COL_TREP, COL_ARM_CTRL), c(lab_m, lab_c))) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  tag_title("I", "IPRS components", "Module against control arm (03h)") +
  labs(x = NULL, y = "Score (0 to 1)", fill = NULL) +
  theme_copula(BASE) +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

gc <- seq(min(dat$x[core_idx]), max(dat$x[core_idx]), length.out = 300)
pc <- fit_curve(sk$fit, gc, se = TRUE)
core_curve <- data.frame(x = gc, fit = as.numeric(pc$fit),
                         lo = as.numeric(pc$fit - 1.96 * pc$se.fit), hi = as.numeric(pc$fit + 1.96 * pc$se.fit))
P$J <- ggplot(pd, aes(x, y)) +
  geom_point(aes(colour = group), alpha = 0.6, size = 0.5) +
  geom_ribbon(data = core_curve, aes(x = x, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = COL_TREP, alpha = 0.2) +
  geom_line(data = core_curve, aes(x = x, y = fit), colour = COL_TREP, linewidth = 1.1) +
  geom_hline(yintercept = ip$ip_y, colour = COL_IPY, linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = Inf, y = ip$ip_y, label = tp(sprintf("IP_{y} = %.3f ", ip$ip_y), bold = TRUE),
           parse = TRUE, hjust = 1, vjust = 1.5, size = TXT(2.8)) +
  scale_colour_manual(values = GROUP_COLS, labels = grp_lab, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1.6))) +
  tag_title("J", "CEP-IP plot", "Pre-IP and post-IP groups (04a)") +
  labs(x = "GOI (log1p)", y = "Composite (log1p)", colour = NULL) +
  theme_copula(BASE) +
  theme(legend.text = element_text(size = BASE - 3), legend.key.size = unit(0.3, "cm")) +
  legend_inside(0.01, 0.99)

P$K <- ggplot(pd, aes(Si, Ki)) +
  geom_point(aes(colour = group), alpha = 0.6, size = 0.5) +
  geom_hline(yintercept = KI_GATE, colour = COL_GATE, linewidth = 0.5) +
  geom_vline(xintercept = SI_GATE, colour = COL_GATE, linewidth = 0.5) +
  annotate("text", x = SI_GATE + 0.03, y = 0.97, label = "SKCS-core", hjust = 0, vjust = 1,
           size = TXT(2.8), fontface = "bold") +
  scale_colour_manual(values = GROUP_COLS, guide = "none") +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  tag_title("K", "Pitcher plot, score space", "K_{i} against S_{i}, gates 0.05 and 0.5 (04b)") +
  labs(x = tx("S_{i}", bold = TRUE), y = tx("K_{i}", bold = TRUE)) +
  theme_copula(BASE)

P$L <- ggplot(pd, aes(u_SK, v_SK)) +
  geom_point(aes(colour = group), alpha = 0.6, size = 0.5) +
  geom_hline(yintercept = vgate, colour = COL_GATE, linewidth = 0.5) +
  geom_vline(xintercept = ugate, colour = COL_GATE, linewidth = 0.5) +
  scale_colour_manual(values = GROUP_COLS, guide = "none") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  tag_title("L", "Pitcher plot, copula space",
            sprintf("v_{SK} against u_{SK}, u_{gate} %.3f, v_{gate} %.3f (04e)", ugate, vgate)) +
  labs(x = tx("u_{SK}", bold = TRUE), y = tx("v_{SK}", bold = TRUE)) +
  theme_copula(BASE)

add_check("G3_panel_count", if (length(P) == 12 && identical(names(P), LETTERS[1:12])) "PASS" else "FAIL",
          sprintf("%d panels: %s", length(P), paste(names(P), collapse = "")))

# ---------------------------------------------------------------------------
# 3. Assemble (grid only, no extra packages)
# ---------------------------------------------------------------------------
NC <- 4; NR <- 3
fn <- file.path(FIG_DIR, "00_overview.png")
png(fn, width = 18, height = 13.8, units = "in", res = 300)
grid::grid.newpage()
lay <- grid::grid.layout(NR + 1, NC, heights = grid::unit(c(0.5, rep(1, NR)), c("in", rep("null", NR))))
grid::pushViewport(grid::viewport(layout = lay))
grid::grid.text("Plot types in this repository, in the order scripts 01 to 04 produce them",
                x = grid::unit(0.01, "npc"), hjust = 0,
                gp = grid::gpar(fontsize = 15, fontface = "bold"),
                vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:NC))
for (i in seq_along(P)) {
  print(P[[i]], vp = grid::viewport(layout.pos.row = (i - 1) %/% NC + 2,
                                    layout.pos.col = (i - 1) %% NC + 1))
}
invisible(dev.off())
cat("Saved 00_overview.png\n")
add_check("G4_file_written", if (file.exists(fn) && file.size(fn) > 0) "PASS" else "FAIL",
          sprintf("00_overview.png, %.1f MB", file.size(fn) / 1e6))

finish_log("00_overview_figure.r")
