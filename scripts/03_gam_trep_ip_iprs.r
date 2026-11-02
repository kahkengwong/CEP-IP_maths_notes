# 03_gam_trep_ip_iprs.r
# Section 3. GAM fit, TREP cells, the inflection point (IP) and the IP
# reliability score (IPRS), run with the cepip package:
#   cepip_fit()  GAM f(x) with the PRSS search over k, donor random
#                intercept and inverse-donor weights
#   cepip_ep()   explanatory power (EP) and the TREP cells
#   cepip_skcs() SKCS-core and its own GAM fit (for the right-hand panels)
#   cepip_ip()   residual GAM R(x), IPx at its first negative-to-positive
#                zero, the Brent 95% CI, the five IPRS components and IPy
# The IPRS is also computed by hand from the parts of the cepip_ip object
# (iprs_by_hand() in 00_setup.r) and checked against the package.
# Figure 03d compares the IP with the calculus definition of an inflection
# point, f''(x) = 0 with a change of sign, on the same fitted curve.
#
# Output (figures/, tables/ and logs/ of the repository folder):
#   03a_GAM_plot.png             composite and control composite against GOI, all cells and SKCS-core
#   03b_PRSS_k_search.png        PRSS by basis dimension k (cepip plot_prss)
#   03c_CEPIP_plot.png           TREP and non-TREP cells, all cells and SKCS-core (cepip plot_fit)
#   03d_calculus_on_fit.png      f, f' and f'' of the fitted curve
#   03e_residual_curve_Rx.png    R(x), IPx and its 95% CI (cepip plot_residual)
#   03f_brent_CI.png             IPx and the two CI bounds as zeros of R(x) and its band edges
#   03g_IPx_to_IPy.png           IPx carried to the composite axis, IPy = f(IPx) (cepip plot_fit)
#   03h_IPRS_components.png      C1 to C5 and IPRS, simulated module against control arm
#   03_IP_summary.csv            k, DE, IP, CI, IPy, components, IPRS and tier for both arms
#   03_console_log.txt           console output of this run
#   03_guardrails.txt            PASS/WARN/FAIL ledger
#
# Run: setwd() to the repository folder, then source this file from scripts/.

source("scripts/00_setup.r")
start_log("03")

dat <- simulate_cells()
n <- nrow(dat)

# ---------------------------------------------------------------------------
# 1. GAM fits and TREP cells (cepip)
# ---------------------------------------------------------------------------
f_all  <- cepip_fit(dat$x, dat$y, donor = dat$donor, id = dat$cell)
ep_all <- cepip_ep(f_all)
sk     <- cepip_skcs(dat$x, dat$y, donor = dat$donor, id = dat$cell)   # SKCS-core and its fit
core   <- dat[sk$cells$is_core, ]
f_core <- sk$fit
f_ctrl_all  <- cepip_fit(dat$x, dat$y_ctrl, donor = dat$donor, id = dat$cell)
f_ctrl_core <- cepip_fit(core$x, core$y_ctrl, donor = core$donor, id = core$cell)
stopifnot(identical(ep_all$id, dat$cell))
dat$trep <- ep_all$is_TREP

de <- c(all = f_all$de, core = f_core$de, ctrl_all = f_ctrl_all$de, ctrl_core = f_ctrl_core$de)
cat(sprintf("k selected by PRSS: all cells %d, SKCS-core %d\n", f_all$k, f_core$k))
cat(sprintf("DE composite: all cells %.1f%%, SKCS-core %.1f%%\n", 100 * de["all"], 100 * de["core"]))
cat(sprintf("DE control  : all cells %.1f%%, SKCS-core %.1f%%\n", 100 * de["ctrl_all"], 100 * de["ctrl_core"]))
cat(sprintf("TREP cells  : all cells %d, SKCS-core %d\n\n", sum(dat$trep), sum(sk$ep$is_TREP)))
add_check("G1_core_raises_DE", if (de["core"] > de["all"]) "PASS" else "WARN",
          sprintf("DE %.3f (all cells) -> %.3f (SKCS-core)", de["all"], de["core"]))
de_pop <- attr(ep_all, "de_pop")
add_check("G2_trep_count", if (sum(dat$trep) == round(n * de_pop)) "PASS" else "FAIL",
          sprintf("TREP n = %d = round(n x DE) = round(%d x %.4f)", sum(dat$trep), n, de_pop))

# ---------------------------------------------------------------------------
# 2. IP and IPRS on all GOI-positive cells (cepip_ip, preprint Methods 2.4)
# ---------------------------------------------------------------------------
ip <- cepip_ip(dat$x, dat$y, donor = dat$donor, id = dat$cell, fit = f_all)
print(ip)
C_pkg <- ip_components(ip)
cat(sprintf("Pre-IP vs post-IP TREP residuals: %s, p = %.3e\n", ip$test$test, ip$test$p))
cat(sprintf("IPRS components: %s\n\n", paste(sprintf("%s=%.3f", names(C_pkg), C_pkg), collapse = ", ")))

add_check("G3_IP_found", if (isTRUE(ip$has_ip)) "PASS" else "FAIL",
          sprintf("%d negative-to-positive crossing(s); the first is the IP", ip$n_crossings))
r_ip <- R_at(ip, ip$ip_x)
add_check("G4_R_zero_at_IP", if (abs(r_ip) < 1e-3) "PASS" else "FAIL",
          sprintf("R(IPx) = %.1e (IPx is interpolated linearly between two grid points)", r_ip))
edge_lo <- R_at(ip, ip$ip_ci[1], +1.96); edge_hi <- R_at(ip, ip$ip_ci[2], -1.96)
add_check("G5_CI_band_edges",
          if (ip$ip_ci[1] < ip$ip_x && ip$ip_x < ip$ip_ci[2] &&
              abs(edge_lo) < 1e-6 && abs(edge_hi) < 1e-6) "PASS" else "FAIL",
          sprintf("CI [%.4f, %.4f] holds IPx; R+1.96SE at lower = %.1e, R-1.96SE at upper = %.1e",
                  ip$ip_ci[1], ip$ip_ci[2], edge_lo, edge_hi))
add_check("G6_components_in_0_1", if (all(C_pkg >= 0 & C_pkg <= 1)) "PASS" else "FAIL",
          paste(sprintf("%s=%.3f", names(C_pkg), C_pkg), collapse = " "))
add_check("G7_IPRS_mean", if (abs(ip$iprs$IPRS - mean(C_pkg)) < 1e-12) "PASS" else "FAIL",
          sprintf("IPRS = mean(C1..C5) = %.4f", ip$iprs$IPRS))
hand <- iprs_by_hand(ip)
d_iprs <- max(abs(hand - c(C_pkg, IPRS = ip$iprs$IPRS)))
add_check("G8_IPRS_by_hand", if (d_iprs < 1e-9) "PASS" else "FAIL",
          sprintf("C1 to C5 and IPRS by hand against cepip_ip(): max difference %.1e", d_iprs))
fy <- as.numeric(fit_curve(f_all, ip$ip_x))
add_check("G9_IPy", if (abs(ip$ip_y - fy) < 1e-9) "PASS" else "FAIL",
          sprintf("IPy = f(IPx) = %.4f (population curve, donor term excluded)", ip$ip_y))
add_check("G10_module_tier", if (ip$iprs$tier %in% c("Strong", "Very Strong")) "PASS" else "WARN",
          sprintf("IPRS %.3f, %s tier", ip$iprs$IPRS, ip$iprs$tier))

# Control arm
ctl  <- simulate_control_arm()
ipc  <- cepip_ip(ctl$x, ctl$y, donor = ctl$donor, id = ctl$cell)
C_ctl <- ip_components(ipc)
cat(sprintf("Control arm: DE %.1f%%, IPx %s, IPRS %.3f (%s tier)\n\n",
            100 * ipc$fit$de, if (isTRUE(ipc$has_ip)) sprintf("%.4f", ipc$ip_x) else "none",
            ipc$iprs$IPRS, ipc$iprs$tier))
add_check("G11_control_lower", if (ipc$iprs$IPRS < ip$iprs$IPRS) "PASS" else "WARN",
          sprintf("IPRS control %.3f < module %.3f", ipc$iprs$IPRS, ip$iprs$IPRS))

summ <- data.frame(arm = c("Simulated module", "Simulated control"),
                   n_cells = c(n, nrow(ctl)), k = c(f_all$k, ipc$fit$k), DE = c(f_all$de, ipc$fit$de),
                   n_TREP = c(sum(ip$is_top), sum(ipc$is_top)),
                   IPx = c(ip$ip_x, ipc$ip_x), CI_lower = c(ip$ip_ci[1], ipc$ip_ci[1]),
                   CI_upper = c(ip$ip_ci[2], ipc$ip_ci[2]), IPy = c(ip$ip_y, ipc$ip_y),
                   test = c(ip$test$test, ipc$test$test), p = c(ip$test$p, ipc$test$p),
                   rbind(C_pkg, C_ctl), IPRS = c(ip$iprs$IPRS, ipc$iprs$IPRS),
                   tier = c(ip$iprs$tier, ipc$iprs$tier), row.names = NULL)
write.csv(summ, file.path(TAB_DIR, "03_IP_summary.csv"), row.names = FALSE)
cat("Saved 03_IP_summary.csv\n")

# ---------------------------------------------------------------------------
# 3. Figure 03a: GAM plot (composite and housekeeping control composite)
# ---------------------------------------------------------------------------
lab_all  <- sprintf("GOI-positive cells (n=%d)", n)
lab_core <- sprintf("SKCS-core (n=%d)", nrow(core))
panel_lv <- c(lab_all, lab_core)
curve_df <- function(fit, xr, panel, what) {
  g <- seq(xr[1], xr[2], length.out = 300)
  p <- fit_curve(fit, g, se = TRUE)
  data.frame(panel = panel, what = what, x = g, fit = as.numeric(p$fit),
             lo = as.numeric(p$fit - 1.96 * p$se.fit), hi = as.numeric(p$fit + 1.96 * p$se.fit))
}
pts <- rbind(data.frame(panel = lab_all,  what = "Control composite", x = dat$x,  y = dat$y_ctrl),
             data.frame(panel = lab_core, what = "Control composite", x = core$x, y = core$y_ctrl),
             data.frame(panel = lab_all,  what = "Composite", x = dat$x,  y = dat$y),
             data.frame(panel = lab_core, what = "Composite", x = core$x, y = core$y))
crv <- rbind(curve_df(f_all, range(dat$x), lab_all, "Composite"),
             curve_df(f_core, range(core$x), lab_core, "Composite"),
             curve_df(f_ctrl_all, range(dat$x), lab_all, "Control composite"),
             curve_df(f_ctrl_core, range(core$x), lab_core, "Control composite"))
pts$panel <- factor(pts$panel, panel_lv); crv$panel <- factor(crv$panel, panel_lv)
de_lab <- data.frame(panel = factor(panel_lv, panel_lv),
                     lab = sprintf("Composite: DE=%.1f%%\nControl: DE=%.1f%%",
                                   100 * de[c("all", "core")], 100 * de[c("ctrl_all", "ctrl_core")]))
what_cols <- c("Composite" = COL_TREP, "Control composite" = COL_CTRL)
p03a <- ggplot() +
  geom_point(data = pts, aes(x, y, colour = what), alpha = 0.25, size = 0.6) +
  geom_ribbon(data = crv, aes(x = x, ymin = lo, ymax = hi, fill = what), alpha = 0.3) +
  geom_line(data = crv, aes(x = x, y = fit, colour = what), linewidth = 1.0) +
  geom_text(data = de_lab, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.2,
            size = TXT(3.0), fontface = "bold") +
  scale_colour_manual(values = what_cols) + scale_fill_manual(values = what_cols, guide = "none") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2, linewidth = 1))) +
  facet_wrap(~ panel) +
  labs(title = "GAM plot",
       subtitle = "Restricting to the SKCS-core raises DE. Curves from cepip_fit(), donor term excluded",
       x = "GOI (log1p)", y = "Composite (log1p)", colour = NULL) +
  theme_copula() +
  legend_inside(0.99, 0.01)
save_plot(p03a, "03a_GAM_plot.png", w = 9, h = 4.6)

# ---------------------------------------------------------------------------
# 4. Figure 03b: PRSS search over k (cepip plot_prss)
# ---------------------------------------------------------------------------
p03b <- restyle_cepip(plot_prss(f_all), title = "PRSS by basis dimension k, all GOI-positive cells")
save_plot(p03b, "03b_PRSS_k_search.png", w = 6, h = 4.6)

# ---------------------------------------------------------------------------
# 5. Figure 03c: CEP-IP plot, all cells and SKCS-core (cepip plot_fit)
# ---------------------------------------------------------------------------
p03c_all  <- restyle_cepip(plot_fit(f_all), title = lab_all, x = "GOI (log1p)", y = "Composite (log1p)")
p03c_core <- restyle_cepip(plot_fit(f_core), title = lab_core, x = "GOI (log1p)", y = "Composite (log1p)")
g03c <- gridExtra::arrangeGrob(p03c_all, p03c_core, ncol = 2,
                               top = grid::textGrob("CEP-IP plot: TREP cells are the top round(n x DE) cells by EP",
                                                    x = 0.01, hjust = 0,
                                                    gp = grid::gpar(fontsize = 12, fontface = "bold")))
ggsave(file.path(FIG_DIR, "03c_CEPIP_plot.png"), g03c, width = 11, height = 5.2, dpi = 300)
cat("Saved 03c_CEPIP_plot.png\n")

# ---------------------------------------------------------------------------
# 6. Figure 03d: the calculus definition on the fitted curve
# f' and f'' by central differences on a 2,000-point grid.
# Candidates: f''(x) = 0 where f'' changes sign (handwritten page A4).
# ---------------------------------------------------------------------------
g  <- seq(min(dat$x), max(dat$x), length.out = 2000)
h  <- diff(g)[1]
fg <- as.numeric(fit_curve(f_all, g))
f1 <- c(NA, (fg[-(1:2)] - fg[1:(length(g) - 2)]) / (2 * h), NA)
f2 <- c(NA, (fg[-(1:2)] - 2 * fg[2:(length(g) - 1)] + fg[1:(length(g) - 2)]) / h^2, NA)
sc2 <- which(diff(sign(f2)) != 0 & !is.na(f2[-1]) & !is.na(f2[-length(f2)]))
x_f2 <- g[sc2] - f2[sc2] * h / (f2[sc2 + 1] - f2[sc2])
cat(sprintf("f''(x) = 0 with a sign change at x = %s\n",
            if (length(x_f2)) paste(sprintf("%.3f", x_f2), collapse = ", ") else "none"))
cat(sprintf("Steepest point of f (largest f'): x = %.3f\n\n", g[which.max(f1)]))
add_check("G12_calculus_candidates", "PASS",
          sprintf("%d sign change(s) of f''; R(x) IPx = %.3f", length(x_f2), ip$ip_x))

lv3 <- c("f(x), fitted composite", "f'(x), slope", "f''(x), curvature")
calc <- rbind(data.frame(panel = lv3[1], x = g, val = fg),
              data.frame(panel = lv3[2], x = g, val = f1),
              data.frame(panel = lv3[3], x = g, val = f2))
calc$panel <- factor(calc$panel, lv3)
zero_df <- data.frame(panel = factor(lv3[3], lv3), y = 0)
p03d <- ggplot(calc, aes(x, val)) +
  geom_hline(data = zero_df, aes(yintercept = y), colour = "gray40", linetype = "dotted") +
  geom_vline(xintercept = x_f2, colour = "grey30", linetype = "dotted", linewidth = 0.5) +
  geom_vline(xintercept = ip$ip_x, colour = COL_GATE, linetype = "dashed", linewidth = 0.8) +
  geom_line(colour = COL_TREP, linewidth = 0.9, na.rm = TRUE) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  labs(title = "Calculus inflection points of the fitted curve",
       subtitle = tx(sprintf(paste0("Dotted: f''(x) = 0 with a change of sign (%d). ",
                                    "Green dashed: IP_{x} from R(x) = %.3f"), length(x_f2), ip$ip_x)),
       x = "GOI (log1p)", y = NULL) +
  theme_copula()
save_plot(p03d, "03d_calculus_on_fit.png", w = 6.5, h = 5.2)   # same shape as 03b, shown side by side

# ---------------------------------------------------------------------------
# 7. Figure 03e: residual curve R(x) (cepip plot_residual)
# ---------------------------------------------------------------------------
s <- ip$iprs
p03e <- restyle_cepip(
  plot_residual(ip),
  title = tx(sprintf("IP_{x} = %.3f (95%% CI: %.3f-%.3f; IPRS=%.3f, %s tier)",
                     ip$ip_x, ip$ip_ci[1], ip$ip_ci[2], s$IPRS, tolower(s$tier)), bold = TRUE),
  subtitle = sprintf("Residual GAM R(x) of the TREP cells. C1 %.2f, C2 %.2f, C3 %.2f, C4 %.2f, C5 %.2f",
                     s$C1, s$C2, s$C3, s$C4, s$C5),
  x = "GOI (log1p)", y = "Residual y - f(x)")
save_plot(p03e, "03e_residual_curve_Rx.png", w = 7, h = 4.8)

# ---------------------------------------------------------------------------
# 8. Figure 03f: IPx and its 95% CI as three roots (CI bounds by Brent's method)
# IPx solves R(x) = 0. The CI bounds solve R(x) + 1.96 SE = 0 (lower) and
# R(x) - 1.96 SE = 0 (upper). Values from cepip_ip(); curves from its grid.
# ---------------------------------------------------------------------------
# Window of 2.5 CI widths on each side of IPx
win <- ip$ip_x + c(-1, 1) * max(0.15, 2.5 * diff(ip$ip_ci))
zg  <- ip$grid[ip$grid$x >= win[1] & ip$grid$x <= win[2], ]
zl  <- rbind(data.frame(x = zg$x, val = zg$r_hat,   line = "R(x)"),
             data.frame(x = zg$x, val = zg$r_upper, line = "R(x) + 1.96 SE"),
             data.frame(x = zg$x, val = zg$r_lower, line = "R(x) - 1.96 SE"))
roots <- data.frame(x = c(ip$ip_ci[1], ip$ip_x, ip$ip_ci[2]), y = 0,
                    lab = vapply(sprintf(c("CI lower\n%.3f", "IP_{x}\n%.3f", "CI upper\n%.3f"),
                                         c(ip$ip_ci[1], ip$ip_x, ip$ip_ci[2])), tp, character(1), bold = TRUE),
                    vj = c(1.3, -0.5, 1.3))   # IPx label above the line, CI labels below
p03f <- ggplot(zl, aes(x, val)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray40") +
  geom_line(aes(linetype = line, linewidth = line), colour = "#D9822B") +
  geom_point(data = roots, aes(x, y), size = 2.2, colour = "black") +
  geom_text(data = roots, aes(x, y, label = lab, vjust = vj), size = TXT(3.0), parse = TRUE) +
  scale_linetype_manual(values = c("R(x)" = "solid", "R(x) + 1.96 SE" = "dashed",
                                   "R(x) - 1.96 SE" = "dotdash")) +
  scale_linewidth_manual(values = c("R(x)" = 1.2, "R(x) + 1.96 SE" = 0.6, "R(x) - 1.96 SE" = 0.6)) +
  labs(title = tx("IP_{x} and its 95% CI as three roots", bold = TRUE),
       subtitle = "CI bounds by Brent's method (uniroot) from the band sign change nearest the IP",
       x = "GOI (log1p)", y = "R(x)", linetype = NULL, linewidth = NULL) +
  theme_copula() +
  legend_inside(0.99, 0.01)
save_plot(p03f, "03f_brent_CI.png", w = 6.5, h = 4.6)

# ---------------------------------------------------------------------------
# 9. Figure 03g: IPx to IPy (cepip plot_fit with the IP)
# ---------------------------------------------------------------------------
p03g <- restyle_cepip(
  plot_fit(ip),
  title = tx("From IP_{x} to IP_{y}", bold = TRUE),
  subtitle = tx(sprintf("IP_{y} = f(IP_{x}) = %.3f is the pre-IP and post-IP threshold. Shaded: 95%% CI of IP_{x}",
                        ip$ip_y)),
  x = "GOI (log1p)", y = "Composite (log1p)")
save_plot(p03g, "03g_IPx_to_IPy.png", w = 6.5, h = 5.6)

# ---------------------------------------------------------------------------
# 10. Figure 03h: IPRS components, module against control arm
# ---------------------------------------------------------------------------
lab_m <- sprintf("Simulated module (%s)", tolower(ip$iprs$tier))
lab_c <- sprintf("Simulated control (%s)", tolower(ipc$iprs$tier))
cmp <- rbind(data.frame(arm = lab_m, comp = c(names(C_pkg), "IPRS"), val = c(C_pkg, ip$iprs$IPRS)),
             data.frame(arm = lab_c, comp = c(names(C_ctl), "IPRS"), val = c(C_ctl, ipc$iprs$IPRS)))
cmp$arm  <- factor(cmp$arm, c(lab_m, lab_c))
cmp$comp <- factor(cmp$comp, c("C1", "C2", "C3", "C4", "C5", "IPRS"))
tiers <- data.frame(y = c(0.2, 0.4, 0.6, 0.8), lab = c("Weak", "Moderate", "Strong", "Very strong"))
p03h <- ggplot(cmp, aes(comp, val, fill = arm)) +
  geom_hline(data = tiers, aes(yintercept = y), colour = "gray50", linetype = "dotted") +
  annotate("text", x = Inf, y = tiers$y, label = tiers$lab,
           hjust = 1.05, vjust = -0.4, size = TXT(2.6), colour = "gray30") +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", val)), position = position_dodge(width = 0.8),
            vjust = -0.4, size = TXT(2.6), fontface = "bold") +
  scale_fill_manual(values = setNames(c(COL_TREP, COL_ARM_CTRL), c(lab_m, lab_c))) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete(expand = expansion(add = c(0.6, 1.3))) +
  labs(title = "IPRS components (cepip_ip)",
       subtitle = "C1 precision, C2 pre/post test, C3 abruptness, C4 band width, C5 smoothness and sign",
       x = NULL, y = "Score (0 to 1)", fill = NULL) +
  theme_copula() +
  theme(legend.position = "bottom")
save_plot(p03h, "03h_IPRS_components.png", w = 7, h = 4.8)

finish_log("03_gam_trep_ip_iprs.r")
