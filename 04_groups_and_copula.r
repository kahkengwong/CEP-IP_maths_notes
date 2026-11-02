# 04_groups_and_copula.r
# Section 4. Pre-IP and post-IP groups, and the IP and the gates in copula
# space (the simulated versions of preprint Figures 3 and 4).
# SKCS-core cells below IPy are pre-IP and the others post-IP. Each group is
# split into TREP and non-TREP cells by the GAM fit of the SKCS-core.
# The IP and the Si and Ki gates are carried into copula space by counting
# cells at or below each value and dividing by n + 1:
#   IPu = #{GOI <= IPx}/(n+1),  IPv = #{composite <= IPy}/(n+1)
#   ugate = #{Si <= 0.05}/(n+1), vgate = #{Ki <= 0.5}/(n+1)
# The scores, the SKCS-core and its TREP cells come from cepip_skcs(), the IP
# from cepip_ip(). Figure 04b is plot_pitcher() from cepip; the copula-space
# figures are drawn here because the package draws them only after CMRF.
#
# Output (written to OUT_DIR):
#   04a_CEPIP_groups.png          composite against GOI by group, with IPy
#   04b_pitcher_groups.png        Ki against Si by group, the two spines (cepip plot_pitcher)
#   04c_expression_to_copula.png  expression space and copula space, with IPy and IPv
#   04d_scores_to_copula.png      score space and its copula space, with the gates
#   04e_pitcher_copula_space.png  pitcher plot in copula space (vSK against uSK) by group
#   04_cell_groups.csv            group, scores and copula coordinates per cell
#   04_console_log.txt            console output of this run
#   04_guardrails.txt             PASS/WARN/FAIL ledger
#
# Run: setwd() to the folder that holds 00_setup.r, then source this file.

source("00_setup.r")
start_log("04")

dat <- simulate_cells()
n <- nrow(dat)
# cepip: SKCS scores, SKCS-core and its GAM fit (TREP cells within the core)
sk <- cepip_skcs(dat$x, dat$y, donor = dat$donor, id = dat$cell)
stopifnot(identical(sk$cells$id, dat$cell))
dat$Si <- sk$cells$Si_norm; dat$Ki <- sk$cells$Ki
dat$core <- sk$cells$is_core
core_idx <- which(dat$core)
stopifnot(identical(sk$ep$id, dat$cell[core_idx]))
dat$trep_core <- FALSE
dat$trep_core[core_idx] <- sk$ep$is_TREP

# cepip: IP on all GOI-positive cells, as in script 03
ip  <- cepip_ip(dat$x, dat$y, donor = dat$donor, id = dat$cell)
IPx <- ip$ip_x; IPy <- ip$ip_y
cat(sprintf("IPx = %.4f, IPy = %.4f, IPRS = %.3f (%s tier)\n", IPx, IPy, ip$iprs$IPRS, ip$iprs$tier))

# ---------------------------------------------------------------------------
# 1. Groups
# ---------------------------------------------------------------------------
dat$group <- "Non-SKCS-core"
side <- ifelse(dat$y < IPy, "Pre-IP", "Post-IP")
dat$group[dat$core] <- paste(side[dat$core], ifelse(dat$trep_core[dat$core], "TREP", "non-TREP"))
grp_lv <- names(GROUP_COLS)
dat$group <- factor(dat$group, grp_lv)
tab <- table(dat$group)
print(tab)
add_check("G1_groups_complete", if (sum(tab) == n && all(!is.na(dat$group))) "PASS" else "FAIL",
          sprintf("%s", paste(sprintf("%s %d", names(tab), as.integer(tab)), collapse = "; ")))
pre_ok <- all(dat$y[grepl("^Pre-IP", dat$group)] < IPy) &&
          all(dat$y[grepl("^Post-IP", dat$group)] >= IPy)
add_check("G2_split_by_IPy", if (pre_ok) "PASS" else "FAIL",
          sprintf("pre-IP below IPy = %.4f, post-IP at or above it", IPy))
add_check("G2b_matches_cepip_labels", if (identical(side, ip$labels)) "PASS" else "FAIL",
          "Pre-IP and post-IP sides equal the labels of cepip_ip() (ip_axis = \"y\")")
grp_lab <- setNames(sprintf("%s (n=%d)", grp_lv, as.integer(tab[grp_lv])), grp_lv)

# ---------------------------------------------------------------------------
# 2. Copula coordinates
# ---------------------------------------------------------------------------
dat$u    <- to_u(dat$x);  dat$v    <- to_u(dat$y)
dat$u_SK <- to_u(dat$Si); dat$v_SK <- to_u(dat$Ki)
IPu   <- sum(dat$x <= IPx) / (n + 1)
IPv   <- sum(dat$y <= IPy) / (n + 1)
ugate <- sum(dat$Si <= SI_GATE) / (n + 1)
vgate <- sum(dat$Ki <= KI_GATE) / (n + 1)
cat(sprintf("IPu = %.4f, IPv = %.4f, ugate = %.4f, vgate = %.4f\n\n", IPu, IPv, ugate, vgate))

core_cop <- dat$u_SK > ugate & dat$v_SK > vgate
add_check("G3_gates_in_copula", if (identical(core_cop, dat$core)) "PASS" else "WARN",
          sprintf("SKCS-core by uSK > ugate and vSK > vgate: %d cells, by Si and Ki: %d cells (differ only with ties)",
                  sum(core_cop), sum(dat$core)))
pre_cop <- dat$v[dat$core] <= IPv
add_check("G4_IPv_split", if (identical(pre_cop, dat$y[dat$core] < IPy)) "PASS" else "WARN",
          "Pre-IP in copula space (v <= IPv) matches pre-IP in expression space (y < IPy)")

write.csv(dat[, c("cell", "x", "y", "u", "v", "Si", "Ki", "u_SK", "v_SK", "core", "trep_core", "group")],
          file.path(OUT_DIR, "04_cell_groups.csv"), row.names = FALSE)
cat("Saved 04_cell_groups.csv\n")

# Worked example cells: a post-IP TREP cell near the middle of its group,
# and the same cell in score space
pick <- which(dat$group == "Post-IP TREP")
ex   <- pick[which.min(abs(dat$u[pick] - 0.75) + abs(dat$v[pick] - 0.75))]
cat(sprintf("Worked cell %s: GOI %.3f -> u %.4f, composite %.3f -> v %.4f, Si %.4f -> uSK %.4f, Ki %.4f -> vSK %.4f\n",
            dat$cell[ex], dat$x[ex], dat$u[ex], dat$y[ex], dat$v[ex],
            dat$Si[ex], dat$u_SK[ex], dat$Ki[ex], dat$v_SK[ex]))

plot_order <- order(dat$group == "Non-SKCS-core", decreasing = TRUE)
pd <- dat[plot_order, ]
group_scale <- scale_colour_manual(values = GROUP_COLS, labels = grp_lab, drop = FALSE)
group_guide <- guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.2)))

# ---------------------------------------------------------------------------
# 3. Figure 04a: CEP-IP plot by group (preprint Figure 3, left)
# ---------------------------------------------------------------------------
gc <- seq(min(dat$x[core_idx]), max(dat$x[core_idx]), length.out = 300)
pc <- fit_curve(sk$fit, gc, se = TRUE)   # SKCS-core GAM from cepip, donor term excluded
core_curve <- data.frame(x = gc, fit = as.numeric(pc$fit),
                         lo = as.numeric(pc$fit - 1.96 * pc$se.fit), hi = as.numeric(pc$fit + 1.96 * pc$se.fit))
p04a <- ggplot(pd, aes(x, y)) +
  geom_point(aes(colour = group), alpha = 0.6, size = 0.8) +
  geom_ribbon(data = core_curve, aes(x = x, ymin = lo, ymax = hi), inherit.aes = FALSE,
              fill = COL_TREP, alpha = 0.2) +
  geom_line(data = core_curve, aes(x = x, y = fit), colour = COL_TREP, linewidth = 1.5) +
  geom_hline(yintercept = IPy, colour = COL_IPY, linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = Inf, y = IPy, label = tp(sprintf("IP_{y} = %.3f ", IPy), bold = TRUE), parse = TRUE,
           hjust = 1, vjust = 1.5, size = TXT(3.6)) +
  group_scale + group_guide +
  labs(title = "CEP-IP plot",
       subtitle = "Purple curve: GAM fit of the SKCS-core, which defines TREP and non-TREP cells",
       x = "GOI (log1p)", y = "Composite (log1p)", colour = NULL) +
  theme_copula() +
  legend_inside(0.01, 0.99)
save_plot(p04a, "04a_CEPIP_groups.png", w = 6, h = 5.2)

# ---------------------------------------------------------------------------
# 4. Figure 04b: pitcher plot by group (preprint Figure 3, right)
# ---------------------------------------------------------------------------
p04b <- restyle_cepip(
  plot_pitcher(sk, ip = ip),
  text_map = list("Si > 0.05" = tp("S_{i} > 0.05", bold = TRUE),
                  "Ki > 0.50" = tp("K_{i} > 0.50", bold = TRUE)),
  subtitle = "The two spines of the SKCS-core are the pre-IP and post-IP groups",
  x = tx("Spearman rank-product contribution (S_{i})", bold = TRUE),
  y = tx("Kendall concordance contribution (K_{i})", bold = TRUE))
save_plot(p04b, "04b_pitcher_groups.png", w = 6, h = 5.6)

# ---------------------------------------------------------------------------
# 5. Figures 04c and 04d: expression and score space into copula space
# (preprint Figure 4)
# ---------------------------------------------------------------------------
# Each panel shows one worked cell (circled) with a label joined to it by a
# short line. place = "upleft" puts the label above and to the left of the
# cell, "upright" above and to the right. Lines are drawn before the label so
# that no gate or IP line crosses the label text.
two_spaces <- function(a1, b1, a2, b2, lv, hl, vl, cell_txt, place, xlab, ylab, title, subtitle) {
  pts <- rbind(data.frame(panel = lv[1], a = pd[[a1]], b = pd[[b1]], group = pd$group),
               data.frame(panel = lv[2], a = pd[[a2]], b = pd[[b2]], group = pd$group))
  pts$panel <- factor(pts$panel, lv)
  ra <- c(diff(range(pd[[a1]])), diff(range(pd[[a2]])))
  rb <- c(diff(range(pd[[b1]])), diff(range(pd[[b2]])))
  ca <- c(dat[[a1]][ex], dat[[a2]][ex]); cb <- c(dat[[b1]][ex], dat[[b2]][ex])
  right <- place == "upright"
  exd <- data.frame(panel = factor(lv, lv), a = ca, b = cb,
                    la = ca + ifelse(right, 0.07, -0.05) * ra, lb = cb + 0.09 * rb,
                    hj = ifelse(right, 0, 1), lab = cell_txt)
  room <- data.frame(panel = factor(lv, lv), a = ca, b = exd$lb + 0.24 * rb)   # space for the label
  hl$panel <- factor(hl$panel, lv)
  p <- ggplot(pts, aes(a, b)) +
    geom_point(aes(colour = group), alpha = 0.55, size = 0.7) +
    geom_hline(data = hl, aes(yintercept = y, linetype = kind), colour = hl$col, linewidth = 0.7)
  if (!is.null(vl)) {
    vl$panel <- factor(vl$panel, lv)
    p <- p + geom_vline(data = vl, aes(xintercept = x), colour = vl$col, linewidth = 0.7)
  }
  p +
    geom_blank(data = room, aes(a, b)) +
    geom_point(data = exd, aes(a, b), shape = 21, size = 4.5, stroke = 1.1, colour = "black", fill = NA) +
    geom_segment(data = exd, aes(x = a, y = b, xend = la, yend = lb), colour = "black", linewidth = 0.4) +
    geom_label(data = exd, aes(la, lb, label = lab, hjust = hj), vjust = 0, size = TXT(2.8),
               parse = TRUE, label.size = 0.25, fill = "white") +
    scale_linetype_manual(values = c(dashed = "dashed", solid = "solid"), guide = "none") +
    facet_wrap(~ panel, scales = "free") +
    group_scale + group_guide +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab, colour = NULL) +
    theme_copula() +
    theme(legend.position = "bottom")
}

lv_c <- c("Expression space", "Copula space")
hl_c <- data.frame(panel = lv_c, y = c(IPy, IPv), kind = "dashed", col = COL_IPY)
txt_c <- c(tp(sprintf("Cell A\nGOI %.3f\ncomposite %.3f", dat$x[ex], dat$y[ex]), bold = TRUE),
           tp(sprintf("Cell A\nu %.4f\nv %.4f", dat$u[ex], dat$v[ex]), bold = TRUE))
p04c <- two_spaces("x", "y", "u", "v", lv_c, hl_c, NULL, txt_c, c("upleft", "upleft"),
                   "GOI (left) or u (right)", "Composite (left) or v (right)",
                   "Expression space and copula space",
                   tx(sprintf("Dashed: IP_{y} = %.3f, carried to IP_{v} = #{composite <= IP_{y}}/(n+1) = %.3f", IPy, IPv)))
save_plot(p04c, "04c_expression_to_copula.png", w = 9.5, h = 5.4)

lv_s <- c("Pitcher plot, score space", "Pitcher plot, copula space")
hl_s <- data.frame(panel = lv_s, y = c(KI_GATE, vgate), kind = "solid", col = COL_GATE)
vl_s <- data.frame(panel = lv_s, x = c(SI_GATE, ugate), col = COL_GATE)
txt_s <- c(tp(sprintf("Cell A\nS_{i} %.4f\nK_{i} %.4f", dat$Si[ex], dat$Ki[ex]), bold = TRUE),
           tp(sprintf("Cell A\nu_{SK} %.4f\nv_{SK} %.4f", dat$u_SK[ex], dat$v_SK[ex]), bold = TRUE))
p04d <- two_spaces("Si", "Ki", "u_SK", "v_SK", lv_s, hl_s, vl_s, txt_s, c("upright", "upleft"),
                   tx("S_{i} (left) or u_{SK} (right)", bold = TRUE),
                   tx("K_{i} (left) or v_{SK} (right)", bold = TRUE),
                   "Pitcher plot in score space and in copula space",
                   tx(sprintf("Gates S_{i} = %.2f and K_{i} = %.1f become u_{gate} = %.3f and v_{gate} = %.3f",
                              SI_GATE, KI_GATE, ugate, vgate)))
save_plot(p04d, "04d_scores_to_copula.png", w = 9.5, h = 5.4)

# ---------------------------------------------------------------------------
# 6. Figure 04e: pitcher plot in copula space (preprint Figure 4C, middle)
# uSK = rank(Si)/(n+1), vSK = rank(Ki)/(n+1). The SKCS-core is the quadrant
# uSK > ugate and vSK > vgate.
# ---------------------------------------------------------------------------
p04e <- ggplot(pd, aes(u_SK, v_SK)) +
  geom_point(aes(colour = group), alpha = 0.6, size = 0.8) +
  geom_hline(yintercept = vgate, colour = COL_GATE, linewidth = 0.6) +
  geom_vline(xintercept = ugate, colour = COL_GATE, linewidth = 0.6) +
  annotate("text", x = 0.01, y = vgate + 0.015, label = tp(sprintf("v_{SK}>%.3f", vgate), bold = TRUE),
           parse = TRUE, hjust = 0, vjust = 0, size = TXT(3.4)) +
  annotate("text", x = ugate + 0.015, y = 0.01, label = tp(sprintf("u_{SK}>%.3f", ugate), bold = TRUE),
           parse = TRUE, hjust = 0, vjust = 0, size = TXT(3.4)) +
  annotate("text", x = ugate + 0.015, y = 0.99, label = "SKCS-core", hjust = 0, vjust = 1,
           size = TXT(3.6), fontface = "bold") +
  group_scale +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Pitcher plot in copula space",
       subtitle = tx(sprintf("Gates S_{i} = %.2f and K_{i} = %.1f become u_{gate} = %.3f and v_{gate} = %.3f",
                             SI_GATE, KI_GATE, ugate, vgate)),
       x = tx("u_{SK} = rank(S_{i})/(n+1)", bold = TRUE), y = tx("v_{SK} = rank(K_{i})/(n+1)", bold = TRUE),
       colour = NULL) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(alpha = 1, size = 2.2))) +
  theme_copula() +
  theme(legend.position = "bottom")
save_plot(p04e, "04e_pitcher_copula_space.png", w = 6.2, h = 7)

finish_log("04_groups_and_copula.r")
