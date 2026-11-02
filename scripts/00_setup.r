# 00_setup.r
# Shared settings, simulated data and helper functions for scripts 01 to 04
# and the overview figure. Each script sources this file first, so all of
# them work on the same simulated cells.
#
# The CEP-IP steps are run with the R package cepip (CRAN, install with
# install.packages("cepip")): cepip_composite(), cepip_fit(), cepip_ep(),
# cepip_ip() and cepip_skcs(), and its plots plot_prss(), plot_fit(),
# plot_residual() and plot_pitcher(). Plots that the package does not
# provide are drawn here with ggplot2. The formulas are also written out by
# hand below (cell_scores(), iprs_by_hand()), and the scripts check that the
# hand calculation gives the same numbers as the package.
#
# The simulated cells imitate one arm of the CEP-IP pipeline: a gene of
# interest (GOI, x) and a composite of partner genes (y), both on a log1p
# scale, in GOI-positive cells. Two groups of cells follow a monotone
# relationship with a gradual step near x = 2.55 (the pre-IP and post-IP
# groups), gentle below the step and steeper above it, and a third group
# between them follows the composite only loosely (background cells).
# A housekeeping composite that does not depend on the GOI serves as the
# control. A second, separate data set with almost no relationship serves as
# the control arm for the IP reliability score (IPRS).
#
# Each cell is also given a donor (20 donors, about equal numbers of cells).
# The donors have no effect on expression in the simulation, but cepip fits
# its donor random intercept and inverse-donor weights as in the published
# pipeline. The composite is the mean of eight simulated partner genes, made
# with cepip_composite().
#
# Requires: cepip, ggplot2 and mgcv.

if (!requireNamespace("cepip", quietly = TRUE))
  stop("Package cepip is not installed. Run install.packages(\"cepip\") first.", call. = FALSE)
suppressPackageStartupMessages({
  library(cepip)
  library(ggplot2)
  library(mgcv)
})
N_DONORS <- 20L

SEED    <- 20261102
# Output folders, relative to the repository folder (the working directory).
# The scripts are in scripts/ and are run from the repository folder, for
# example source("scripts/01_ranks_and_PIT.r").
FIG_DIR <- "figures"
TAB_DIR <- "tables"
LOG_DIR <- "logs"
for (d in c(FIG_DIR, TAB_DIR, LOG_DIR)) dir.create(d, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Colours
# Basic plots (scripts 01 and 02): the colours of the first version.
# Per-cell score maps (02b, 02c): the viridis turbo palette (blue low, red
# high), scaled to the 5th to 95th percentile of the values shown.
# Histograms: fill and border of the Part 12 donor histograms (PVM_3 colour
# #F2A950 at 45% on white, border at 70% brightness), with a kernel density
# curve in the border colour.
# Pipeline plots (scripts 02 to 04): the colours of the preprint figures.
# ---------------------------------------------------------------------------
COL_POINT    <- "black"
COL_TIE      <- "#D55E00"
COL_OTHER    <- "grey30"
COL_LINE_RED <- "#B2182B"

HIST_FILL    <- "#F9D8AF"
HIST_BORDER  <- "#A97638"
HIST_BAR_LW  <- 0.32   # lwd 0.9 in base graphics
HIST_CURVE_LW <- 0.70  # lwd 2 in base graphics

COL_TREP     <- "#4B0082"   # TREP cells and the SKCS-core
COL_NONTREP  <- "gray80"    # non-TREP and non-SKCS-core cells
COL_GAM      <- "#FFCC99"   # GAM fit and R(x)
COL_CTRL     <- "#4F75DE"   # housekeeping control composite
COL_GATE     <- "#009E73"   # Si and Ki gates, IPx
COL_IPY      <- "grey10"    # IPy line
COL_ARM_CTRL <- "#B8B8B8"   # control arm
GROUP_COLS <- c("Pre-IP TREP"      = "#E69F00",
                "Pre-IP non-TREP"  = "#1A237E",
                "Post-IP TREP"     = "#009E73",
                "Post-IP non-TREP" = "#D94F4F",
                "Non-SKCS-core"    = "gray80")

theme_copula <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(plot.title       = element_text(face = "bold", size = base_size + 1),
          plot.subtitle    = element_text(size = base_size - 1, color = "gray30"),
          axis.title       = element_text(face = "bold", size = base_size + 0.5),
          axis.text        = element_text(size = base_size - 2),
          strip.text       = element_text(face = "bold", size = base_size - 1),
          legend.title     = element_text(face = "bold", size = base_size - 1),
          legend.text      = element_text(size = base_size - 1),
          panel.grid.major = element_line(color = "gray90"),
          panel.grid.minor = element_line(color = "gray95"),
          panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.35),
          plot.background  = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
}

# Turbo colour scale over the 5th to 95th percentile of v. Values outside
# that range take the end colours (blue below, red above).
scale_turbo_q <- function(v, name, probs = c(0.05, 0.95)) {
  lim <- as.numeric(stats::quantile(v, probs, na.rm = TRUE))
  scale_colour_viridis_c(option = "turbo", direction = 1, limits = lim,
                         oob = scales::squish, name = name)
}

# ---------------------------------------------------------------------------
# Plot text with subscripts, superscripts and Greek letters.
# Markup: X_{sub} subscript, X^{sup} superscript, \rho \tau etc. Greek.
# Example: "Gates S_{i} = 0.05 and K_{i} = 0.5" or "\\tau_{b} = mean(K_{i})".
# tp() returns a plotmath string for layers with parse = TRUE.
# tx() returns an expression for labs(), legend titles and scale labels.
# Lines are separated by "\n"; several lines are stacked with atop().
.pm_q <- function(s) paste0('"', gsub('"', '\\\\"', gsub("\\\\", "\\\\\\\\", s)), '"')
.pm_line <- function(s, bold) {
  pat <- "(\\\\[a-zA-Z]+|[A-Za-z0-9]+)([_^])\\{([^}]*)\\}|\\\\([a-zA-Z]+)"
  out <- character(0); rest <- s
  wrap <- function(e) if (bold) sprintf("bold(%s)", e) else e
  while (nchar(rest) > 0) {
    m <- regexpr(pat, rest, perl = TRUE)
    if (m == -1) { out <- c(out, wrap(.pm_q(rest))); break }
    if (m > 1) out <- c(out, wrap(.pm_q(substr(rest, 1, m - 1))))
    tok <- regmatches(rest, m)
    g <- regmatches(tok, regexec(pat, tok, perl = TRUE))[[1]]
    if (nzchar(g[2])) {
      base <- if (startsWith(g[2], "\\")) substring(g[2], 2) else .pm_q(g[2])
      op <- if (g[3] == "_") "[%s]" else "^%s"
      out <- c(out, wrap(paste0(base, sprintf(op, .pm_q(g[4])))))
    } else {
      out <- c(out, wrap(g[5]))
    }
    rest <- substring(rest, m + attr(m, "match.length"))
  }
  paste(out, collapse = "*")
}
tp <- function(s, bold = FALSE) {
  lines <- strsplit(s, "\n", fixed = TRUE)[[1]]
  pl <- vapply(lines, .pm_line, character(1), bold = bold)
  pl <- unname(pl)
  if (length(pl) == 1) return(pl)
  acc <- sprintf("displaystyle(%s)", pl[length(pl)])
  for (k in rev(seq_len(length(pl) - 1))) acc <- sprintf("atop(displaystyle(%s), %s)", pl[k], acc)
  acc
}
tx <- function(s, bold = FALSE) parse(text = tp(s, bold))[[1]]

# Size of text drawn inside panels (mm), scaled by 1.15 from the earlier version
TXT <- function(mm) round(mm * 1.15, 2)

# Legend inside the panel. ggplot2 >= 3.5.0 uses legend.position.inside.
legend_inside <- function(x, y) {
  box <- theme(legend.background = element_rect(fill = "white", color = "gray70", linewidth = 0.3),
               legend.justification = c(x, y))
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    box + theme(legend.position = "inside", legend.position.inside = c(x, y))
  } else {
    box + theme(legend.position = c(x, y))
  }
}

# ---------------------------------------------------------------------------
# Logging and guardrail ledger
# ---------------------------------------------------------------------------
start_log <- function(prefix) {
  LOG_PREFIX <<- prefix
  sink(file.path(LOG_DIR, paste0(prefix, "_console_log.txt")), split = TRUE)
  cat("Run started:", format(Sys.time()), "\n")
  cat("R version:", R.version.string,
      " ggplot2:", as.character(packageVersion("ggplot2")),
      " mgcv:", as.character(packageVersion("mgcv")), "\n\n")
  ledger <<- data.frame(id = character(), status = character(), detail = character())
}

add_check <- function(id, status, detail) {
  ledger[nrow(ledger) + 1, ] <<- c(id, status, detail)
  cat(sprintf("[%s] %s: %s\n", status, id, detail))
}

finish_log <- function(script_name) {
  writeLines(c(sprintf("%s guardrails, %s", script_name, format(Sys.time())),
               sprintf("%-24s %-5s %s", ledger$id, ledger$status, ledger$detail),
               sprintf("Summary: %d PASS, %d WARN, %d FAIL",
                       sum(ledger$status == "PASS"), sum(ledger$status == "WARN"),
                       sum(ledger$status == "FAIL"))),
             file.path(LOG_DIR, paste0(LOG_PREFIX, "_guardrails.txt")))
  cat("\nGuardrail summary:", sum(ledger$status == "PASS"), "PASS,",
      sum(ledger$status == "WARN"), "WARN,", sum(ledger$status == "FAIL"), "FAIL\n")
  cat("Run finished:", format(Sys.time()), "\n")
  sink()
}

save_plot <- function(p, name, w = 6, h = 5) {
  ggsave(file.path(FIG_DIR, name), p, width = w, height = h, dpi = 300)
  cat("Saved", name, "\n")
}

# ---------------------------------------------------------------------------
# Histogram with a kernel density curve, in percent of cells.
# The density uses the bw.nrd0 bandwidth and is reflected at finite bounds
# (lo, hi) so that it does not drop to zero at a hard edge. It is scaled to
# the bar axis: percent per bin = 100 x bin width x density.
# ---------------------------------------------------------------------------
hist_frames <- function(v, breaks, lo = -Inf, hi = Inf, panel = "") {
  h  <- hist(v, breaks = breaks, plot = FALSE)
  bw <- diff(breaks)[1]
  bars <- data.frame(panel = panel, xmin = head(breaks, -1), xmax = tail(breaks, -1),
                     pct = 100 * h$counts / length(v))
  aug <- v
  if (is.finite(lo)) aug <- c(aug, 2 * lo - v)
  if (is.finite(hi)) aug <- c(aug, 2 * hi - v)
  k  <- length(aug) / length(v)
  xs <- seq(min(breaks), max(breaks), length.out = 512)
  d  <- stats::density(aug, bw = stats::bw.nrd0(v), from = min(xs), to = max(xs), n = 512)
  dens <- data.frame(panel = panel, x = xs, pct = 100 * bw * k * stats::approx(d$x, d$y, xs)$y)
  list(bars = bars, dens = dens)
}

hist_layers <- function(bars, dens) {
  list(geom_rect(data = bars, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = pct),
                 fill = HIST_FILL, colour = HIST_BORDER, linewidth = HIST_BAR_LW,
                 inherit.aes = FALSE),
       geom_line(data = dens, aes(x = x, y = pct), colour = HIST_BORDER,
                 linewidth = HIST_CURVE_LW, inherit.aes = FALSE))
}

# ---------------------------------------------------------------------------
# Simulated cells
# rtnorm: normal draws truncated below at lo, by inversion (no ties).
# pos: smooth positive map, close to the identity above 0.5, so the
# composite stays positive without clipping (clipping would create ties).
# ---------------------------------------------------------------------------
rtnorm <- function(n, mean, sd, lo) stats::qnorm(stats::runif(n, stats::pnorm(lo, mean, sd), 1), mean, sd)
pos    <- function(z) log1p(exp(4 * z)) / 4

simulate_cells <- function(seed = SEED) {
  set.seed(seed)
  n_pre <- 650; n_post <- 750; n_bg <- 550
  x_pre  <- rtnorm(n_pre,  1.8, 0.6, 0.1)
  x_post <- rtnorm(n_post, 3.4, 0.9, 0.1)
  x_mod  <- c(x_pre, x_post)
  # A gentle slope below the step and a steeper one above it. The two groups
  # overlap in GOI and the step is gradual (logistic scale 0.35), which keeps
  # the edges of the pitcher plot smooth in score space and in copula space.
  f_true <- function(x) 0.5 + 0.2 * x + 1.3 * stats::plogis((x - 2.55) / 0.35) +
    0.45 * pmax(0, x - 2.55)
  y_mod  <- pos(f_true(x_mod) + stats::rnorm(length(x_mod), 0, rep(c(0.30, 0.25), c(n_pre, n_post))))
  # Background cells sit between the two groups and follow the composite
  # only loosely
  x_bg   <- rtnorm(n_bg, 2.6, 0.7, 0.1)
  y_bg   <- pos(1.0 + 0.25 * x_bg + stats::rnorm(n_bg, 0, 0.70))
  x <- c(x_mod, x_bg)
  y_true    <- c(y_mod, y_bg)
  ctrl_true <- pos(1.0 + stats::rnorm(length(x), 0, 0.35))
  n <- length(x)
  # Eight partner genes per cell whose mean equals the composite exactly,
  # and eight housekeeping genes whose mean equals the control composite.
  partners <- function(target, prefix) {
    u <- matrix(stats::runif(n * 8, -0.5, 0.5), n, 8)
    w <- 1 + 0.6 * (u - rowMeans(u))          # row means of w are exactly 1
    m <- target * w
    colnames(m) <- paste0(prefix, 1:8)
    m
  }
  P <- partners(y_true, "P"); H <- partners(ctrl_true, "H")
  donor <- sprintf("D%02d", sample(rep(seq_len(N_DONORS), length.out = n)))
  d <- data.frame(cell = sprintf("cell_%04d", seq_len(n)), donor = donor, x = x,
                  P, H, source = rep(c("pre group", "post group", "background"),
                                     c(n_pre, n_post, n_bg)),
                  stringsAsFactors = FALSE)
  d$y      <- cepip_composite(d, paste0("P", 1:8))   # composite of the partner genes
  d$y_ctrl <- cepip_composite(d, paste0("H", 1:8))   # housekeeping control composite
  attr(d, "y_true") <- y_true
  d
}

# Control arm: the composite barely depends on the GOI.
simulate_control_arm <- function(seed = SEED + 7) {
  set.seed(seed)
  x <- c(rtnorm(1000, 1.6, 0.6, 0.1), rtnorm(1000, 3.4, 0.9, 0.1))
  data.frame(cell = sprintf("ctrl_%04d", seq_along(x)),
             donor = sprintf("D%02d", sample(rep(seq_len(N_DONORS), length.out = length(x)))),
             x = x, y = pos(1.4 + 0.08 * x + stats::rnorm(length(x), 0, 0.5)))
}

to_u <- function(z) rank(z, ties.method = "average") / (length(z) + 1)

# ---------------------------------------------------------------------------
# Per-cell scores (preprint Methods 2.3)
#   Si = 4 (Rx_i - Rbar)(Ry_i - Rbar) / (n - 1)^2
#   Ki = (C_i - D_i) / (n - 1), with C_i and D_i the concordant and
#        discordant pairs that involve cell i
# ---------------------------------------------------------------------------
cell_scores <- function(x, y) {
  n  <- length(x)
  Rx <- rank(x, ties.method = "average")
  Ry <- rank(y, ties.method = "average")
  Rb <- (n + 1) / 2
  Si <- 4 * (Rx - Rb) * (Ry - Rb) / (n - 1)^2
  Ki <- vapply(seq_len(n), function(i)
    sum(sign(x[-i] - x[i]) * sign(y[-i] - y[i])), numeric(1)) / (n - 1)
  data.frame(Rx = Rx, Ry = Ry, Si = Si, Ki = Ki)
}

SI_GATE <- 0.05
KI_GATE <- 0.5

# ---------------------------------------------------------------------------
# Helpers around cepip objects
# pop_predict(): population-level curve of a cepip GAM (donor term excluded),
# as cepip uses internally for f(x) and R(x).
# ---------------------------------------------------------------------------
pop_predict <- function(model, x, donor_levels, se = FALSE) {
  nd <- data.frame(GOI = x)
  has_re <- any(vapply(model$smooth, function(sm) sm$label, character(1)) == "s(donor)")
  if (has_re) {
    nd$donor <- factor(donor_levels[1], levels = donor_levels)
    return(stats::predict(model, newdata = nd, se.fit = se, exclude = "s(donor)"))
  }
  stats::predict(model, newdata = nd, se.fit = se)
}
fit_curve <- function(fit, x, se = FALSE) {
  pop_predict(fit$model, x, levels(factor(fit$data$donor)), se = se)
}
R_at <- function(ip, x, mult = 0) {
  p <- pop_predict(ip$resid_model, x, levels(factor(ip$residual_data$donor)), se = TRUE)
  as.numeric(p$fit) + mult * as.numeric(p$se.fit)
}

# ---------------------------------------------------------------------------
# IPRS by hand (Supplementary Methods S5), from the parts of a cepip_ip
# object: the R(x) grid, the 95% CI, the test p-value and the TREP residuals.
#   C1 precision    1 - CI width / x range
#   C2 statistical  tanh(s x -log10 p), s set so that p = 0.05 gives 0.5
#   C3 abruptness   tanh(|R'(IPx)|), slope by central difference on the grid
#   C4 ribbon       exp(-mean band width / amplitude of R(x))
#   C5 smoothness   exp(-max(0, TV / amplitude - 1)) x (1 - post-IP negative share)^2
#   IPRS = (C1 + C2 + C3 + C4 + C5) / 5
# ---------------------------------------------------------------------------
iprs_by_hand <- function(ip) {
  g <- ip$grid; xs <- g$x; r <- g$r_hat
  C1 <- max(0, 1 - diff(ip$ip_ci) / diff(range(xs)))
  C2 <- tanh(atanh(0.5) / (-log10(0.05)) * (-log10(max(ip$test$p, 1e-300))))
  i  <- min(max(which.min(abs(xs - ip$ip_x)), 2L), length(xs) - 1L)
  C3 <- tanh(abs((r[i + 1] - r[i - 1]) / (xs[i + 1] - xs[i - 1])))
  amp <- diff(range(r))
  C4 <- exp(-mean(g$r_upper - g$r_lower) / (amp + 1e-6))
  rd  <- ip$residual_data
  neg <- mean(rd$residual[rd$GOI >= ip$ip_x] < 0)
  C5 <- exp(-max(0, sum(abs(diff(r))) / (amp + 1e-6) - 1)) * (1 - neg)^2
  C <- c(C1 = C1, C2 = C2, C3 = C3, C4 = C4, C5 = C5)
  c(C, IPRS = mean(C))
}
ip_components <- function(ip) unlist(ip$iprs[c("C1", "C2", "C3", "C4", "C5")])

# ---------------------------------------------------------------------------
# restyle_cepip(): applies the conventions of this repository to a plot made
# by cepip. Text drawn inside the panel is enlarged by the same factor as in
# the other figures, plain labels listed in text_map are replaced by
# subscripted versions, and titles and axis titles can be replaced.
# Layers of a ggplot are edited in place; if that fails (for example with a
# future ggplot2), the plot is returned with only the titles changed.
# ---------------------------------------------------------------------------
restyle_cepip <- function(p, text_map = list(), title = NULL, subtitle = NULL, x = NULL, y = NULL) {
  tryCatch({
    for (l in p$layers) {
      if (inherits(l$geom, "GeomText") || inherits(l$geom, "GeomLabel")) {
        lab <- l$aes_params$label
        if (is.character(lab) && length(lab) == 1 && lab %in% names(text_map)) {
          l$aes_params$label <- text_map[[lab]]
          l$geom_params$parse <- TRUE
        }
        if (is.numeric(l$aes_params$size)) l$aes_params$size <- TXT(l$aes_params$size)
      }
    }
  }, error = function(e) message("restyle_cepip: layers left unchanged (", conditionMessage(e), ")"))
  # One labs() call per title: do.call() would evaluate the plotmath calls
  if (!is.null(title))    p <- p + labs(title = title)
  if (!is.null(subtitle)) p <- p + labs(subtitle = subtitle)
  if (!is.null(x))        p <- p + labs(x = x)
  if (!is.null(y))        p <- p + labs(y = y)
  p + theme(axis.title = element_text(face = "bold", size = 10.5))
}
