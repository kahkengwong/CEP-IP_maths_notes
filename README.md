# Rank correlation, copulas and the inflection point: notes behind CEP-IP

These notes record the mathematics behind the CEP-IP framework and the CMRF filter. Each part pairs handwritten working from my self-study of calculus and statistics (2024 to 2026) with R scripts that reproduce the main identities on simulated data.

The handwritten pages are photographs of the original notebooks, shown two per row. Click a page to open it at full size. Corrections found later are added as typed notes below each page; the pages themselves are left unedited.

## Overview of the plots

The figure below shows one example of each plot type, in the order scripts 01 to 04 produce them. It is drawn by `scripts/00_overview_figure.r`.

![](figures/00_overview.png)

## Contents of this repository

| Part | Files | Content |
|---|---|---|
| A | `handwritten/hw_A*.jpg` | Calculus: limits, derivatives, smoothness, concavity, inflection points, integrals, small changes |
| B | `handwritten/hw_B*.jpg` | Statistics: sums of squares, penalised regression, Pearson, Spearman and Kendall |
| 0 | `scripts/00_setup.r` | Shared settings, simulated cells and functions, sourced by all other scripts |
| 0 | `scripts/00_overview_figure.r`, `figures/00_overview.png` | One panel for each plot type in scripts 01 to 04 |
| 1 | `scripts/01_ranks_and_PIT.r`, `figures/01*.png` | Ranks and the probability integral transform |
| 2 | `scripts/02_spearman_kendall_scores.r`, `figures/02*.png` | Spearman's rho and Kendall's tau, the per-cell scores S<sub>i</sub> and K<sub>i</sub>, and the pitcher plot |
| 3 | `scripts/03_gam_trep_ip_iprs.r`, `figures/03*.png` | GAM fit with the PRSS search over k, TREP cells, the inflection point and the IP reliability score (IPRS) |
| 4 | `scripts/04_groups_and_copula.r`, `figures/04*.png` | Pre-IP and post-IP groups, and the IP and gates in copula space |

Folders:

- `scripts/`. The R scripts.
- `figures/`. All figures.
- `tables/`. Per-cell values and the IP summary (CSV).
- `logs/`. The console log (`*_console_log.txt`) and the guardrail ledger (`*_guardrails.txt`) of each script.
- `handwritten/`. Photographs of the notebook pages.

The scripts need R, ggplot2, mgcv and the R package [cepip](https://github.com/kahkengwong/cepip), available on CRAN: The scripts need R, ggplot2, mgcv and the R package [cepip](https://github.com/kahkengwong/cepip), available on CRAN:

```r
install.packages("cepip")
```

- The CEP-IP steps are run with cepip: `cepip_composite()`, `cepip_fit()`, `cepip_ep()`, `cepip_skcs()` and `cepip_ip()`.
- The figures that cepip already draws use its plot functions: `plot_prss()`, `plot_fit()`, `plot_residual()` and `plot_pitcher()`.
- The other figures are drawn with ggplot2.
- The formulas for S<sub>i</sub>, K<sub>i</sub> and the five IPRS components are also written out by hand in `scripts/00_setup.r`. The guardrails check that the hand calculation gives the same numbers as the package.

The simulated cells come from 20 donors, so cepip fits its donor random intercept and inverse-donor weights as in the published pipeline. The donors have no effect on expression in the simulation.

## Running the scripts

Set the working directory to the repository folder, then run the scripts in order. Each script sources `scripts/00_setup.r` and rebuilds the same simulated cells, so the scripts can also be run one at a time.

```r
setwd("path/to/this/repository")
source("scripts/00_overview_figure.r")
source("scripts/01_ranks_and_PIT.r")
source("scripts/02_spearman_kendall_scores.r")
source("scripts/03_gam_trep_ip_iprs.r")
source("scripts/04_groups_and_copula.r")
```

The figures are written to `figures/`, the CSV tables to `tables/`, and the console logs and guardrail ledgers to `logs/`.

---

## Part A. Calculus foundations

### A1. Limits and the derivative from first principles

The derivative is the limit of the gradient of a secant line as the second point moves onto the first.

1. A limit describes the value a function approaches near a point, whether or not the function is defined at that point (`hw_A1a`, `hw_A1b`).
2. A limit exists only when the left and right limits agree. For |x|/x at 0 the right limit is 1 and the left limit is −1, so the limit does not exist (`hw_A1a`).
3. The gradient between two points on a curve is δy/δx. As x₂ → x₁, δx → 0, and the secant becomes the tangent (`hw_A1c`, `hw_A1d`).
4. The gradient function is dy/dx = lim (δx → 0) δy/δx, written in standard form as f'(x) = lim (h → 0) [f(x + h) − f(x)] / h (`hw_A1d`).

<img src="handwritten/hw_A1a_limit_relation_to_differentiation.jpg" width="49%"> <img src="handwritten/hw_A1b_limits_concept.jpg" width="49%">
<img src="handwritten/hw_A1c_first_principles_gradient_function.jpg" width="49%"> <img src="handwritten/hw_A1d_secant_to_tangent.jpg" width="49%">

### A2. Continuity, and when a derivative exists

A derivative can only be read from a curve that is continuous and smooth at that point.

1. Continuity at x = a needs three conditions: f(a) is defined, the limit as x → a exists, and the limit equals f(a). Each condition is checked in turn (`hw_A2a`, `hw_A2b`).
2. A derivative fails to exist at a jump, at a sharp corner or cusp (for example |x| at 0), and at a vertical tangent (for example x^(1/3) at 0) (`hw_A2c`).
3. d²y/dx² is the operator d/dx applied twice to y. Functions are graded by how many continuous derivatives they have: C⁰ is continuous, C¹ has one continuous derivative, C² has two, and C∞ has infinitely many (`hw_A2d`).

Relevance to CEP-IP: the fitted curve is a penalised regression spline, which is C² over the observed range. Its first and second derivatives therefore exist at every point used in the analysis.

<img src="handwritten/hw_A2a_continuity_definition_step1.jpg" width="49%"> <img src="handwritten/hw_A2b_continuity_steps2_3.jpg" width="49%">
<img src="handwritten/hw_A2c_when_a_derivative_exists.jpg" width="49%"> <img src="handwritten/hw_A2d_second_derivative_notation_and_Ck_smoothness.jpg" width="49%">

### A3. The product rule from first principles

Expanding (u + δu)(v + δv) and dividing by δx gives three terms. The term δu·δv/δx goes to zero as δx → 0, leaving d(uv)/dx = u dv/dx + v du/dx. The page notes that this step needs u and v to be continuous.

<p align="center"><img src="handwritten/hw_A3_product_rule_from_first_principles.jpg" width="49%"></p>

### A4. Second derivative, concavity and inflection points

An inflection point is where the curve changes from concave up to concave down, or the reverse.

1. f''(x) > 0 means the gradient is increasing (concave up). f''(x) < 0 means it is decreasing (concave down) (`hw_A4a`).
2. f''(x) = 0 marks a candidate only. The point is an inflection point if f'' changes sign across it (`hw_A4a`).
3. Worked example, f(x) = x / (x² − 4): f''(x) = 2x(x² + 12) / (x² − 4)³, so f''(0) = 0, and f'' changes sign at x = 0. The point (0, 0) is an inflection point (`hw_A4b`).
4. Curve sketching follows a fixed order: domain, key points, edges, asymptotes, extrema and inflection points. A sign chart of f, f' and f'' over each interval gives the shape of every branch (`hw_A4c`, `hw_A4d`).

<img src="handwritten/hw_A4a_second_derivative_concavity.jpg" width="49%"> <img src="handwritten/hw_A4b_inflection_point_worked_example.jpg" width="49%">
<img src="handwritten/hw_A4c_curve_sketching_steps.jpg" width="49%"> <img src="handwritten/hw_A4d_curve_sketching_sign_chart.jpg" width="49%">

> Note added later. On `hw_A4b` the line "for x < 0, f''(x) > 0" holds on −2 < x < 0 only. For x < −2, f''(x) < 0. The sign chart on `hw_A4d` has the correct intervals.

> Note added later. Point 6 on `hw_A4d` reads "if f''(x) = 0 has solutions, concavity changes". The sign change is also required. Counterexample: y = x⁴ has f''(0) = 0 but stays concave up.

### A5. From a Riemann sum to an integral

A volume of revolution is built from thin cylinders of volume πy²δx. The page first sums five cylinders numerically, then lets n → ∞ and δx → 0 to obtain V = ∫ πy² dx.

Relevance: the curvature penalty ∫ f''(x)² dx in Part B1 is an integral of this kind. It sets how smooth the GAM fit in script 03 is.

<p align="center"><img src="handwritten/hw_A5_riemann_sum_to_integral.jpg" width="49%"></p>

### A6. Two questions that needed correction

- **A cubic displacement with no maximum** (`hw_A6a`). The question asked for the maximum displacement of s(t) = 3t³ + (23/4)t² − 3t − 3. The only turning point for t > 0 is a local minimum (s'' > 0 there), and s → ∞ as t → ∞, so no maximum exists. The page rewrites the question as "the displacement at the critical point after passing X".
- **Disc against cylinder method** (`hw_A6b`). For y = x² + 2 rotated about the y-axis, the y-values from 0 to 2 are not bounded by the curve, so ∫ πx² dy over that range does not give the volume. The cylinder must be computed separately and the remainder subtracted.

<img src="handwritten/hw_A6a_textbook_question_corrected.jpg" width="49%"> <img src="handwritten/hw_A6b_disc_vs_cylinder_method.jpg" width="49%">

### A7. Small changes and approximations

For a small change δx, the change in y is approximately δy ≈ (dy/dx) δx (`hw_A7`). Example 2 on the page uses this to estimate ∛63.98 ≈ 4 − 1/2400 ≈ 3.9996.

Relevance to the scripts: the same ratio of small changes, read the other way, estimates a derivative. Script 03 computes f'(x) and f''(x) of the fitted curve, and the slope R'(IPx) used in IPRS component C3, from differences between neighbouring points of a 2,000-point grid. IPx is located between the two neighbouring grid points where R(x) changes sign from negative to positive. The bounds of its 95% CI are bracketed the same way and then refined with Brent's method.

<p align="center"><img src="handwritten/hw_A7_small_changes_approximation.jpg" width="49%"></p>

> Note added later. The heading defines δx = x_old − x_new, with x_old as the target value (63.98 in Example 2). Question 17(a) writes δx = x_new − x_old = 1.05 − 1, with x_new as the target. Both answers are correct, because in both cases δx is the target value minus the point where dy/dx is evaluated. In 17(b), dx/dy = 4/y² is evaluated at y = 7.96. The usual choice is the known point y = 8, which gives δx ≈ (1/16)(−0.04) = −0.0025. Both give −0.0025 to two significant figures.

---

## Part B. Statistics foundations

### B1. Sums of squares and the penalised residual sum of squares

Every sum of squares has the form Σ(observed − reference)². The reference changes with what is being decomposed: SS_total = SS_model + SS_residual.

The penalised residual sum of squares is

PRSS = Σ (yᵢ − f(xᵢ))² + λ ∫ f''(x)² dx.

The first term measures the distance of the data from the fitted curve. The second term is the integrated squared second derivative, which penalises curvature. λ sets the balance between fit and smoothness.

Relevance to CEP-IP: this is the criterion minimised to fit the CEP-IP curve, with λ chosen by REML. The f'' in the penalty is the same second derivative as in Part A4.

The page also has the standard deviation, the 68-95-99.7 rule for the normal distribution and the Shapiro-Wilk statistic W. Script 03 uses the Shapiro-Wilk test to choose the test for IPRS component C2: Welch's t-test if the pre-IP and post-IP TREP residuals both pass, and the Mann-Whitney U test otherwise.

<p align="center"><img src="handwritten/hw_B1_sums_of_squares_and_PRSS.jpg" width="49%"></p>

### B2. Pearson's r and its denominator

r = SS_xy / √(SS_xx · SS_yy). The denominator is the geometric mean of the two sums of squares, which returns the covariance to a unitless scale between −1 and +1.

<p align="center"><img src="handwritten/hw_B2_pearson_denominator.jpg" width="49%"></p>

### B3. Spearman's rho, ranks and ties

1. Spearman's rho is Pearson's r computed on the ranks.
2. With no ties, ρ = 1 − 6Σdᵢ² / (n(n² − 1)), where dᵢ is the rank difference.
3. Tied values receive the average of the ranks they span (worked example with two tied scores of 0).
4. The mean rank is (n + 1)/2, from 1 + 2 + … + n = n(n + 1)/2 (derived on the page by writing the sum forwards and backwards).

Script 01 shows the effect of ties in copula space (`figures/01e_ties_zero_GOI.png`). Script 02 checks that ρ equals 3 × mean(S<sub>i</sub>) × (n − 1)/(n + 1) and that τ<sub>b</sub> equals mean(K<sub>i</sub>), both to machine precision.

<p align="center"><img src="handwritten/hw_B3_spearman_ranks_and_ties.jpg" width="49%"></p>

### B4. Spearman's rho as a sum of per-cell terms

Writing the numerator of rho as Σ (R_x,i − R̄)(R_y,i − R̄) shows that rho is a sum of one term per cell. Each cell's term is the product of its two rank deviations. This is the starting point for the per-cell Spearman score used in CEP-IP.

<p align="center"><img src="handwritten/hw_B4_spearman_per_cell_terms.jpg" width="49%"></p>

> Note added later. The boxed line ρ = ΣSᵢ / (σ_Rx · σ_Ry) holds when Sᵢ includes the factor 1/(n − 1). With Sᵢ defined as the rank product (R_x,i − R̄)(R_y,i − R̄), as on this page, ρ = ΣSᵢ / ((n − 1) σ_Rx σ_Ry) = ΣSᵢ / √(SS_xx · SS_yy). The three equivalent forms at the bottom of the page are correct.

### B5. Kendall's tau as a sum of per-cell terms

Kendall's tau splits into per-cell terms in the same way as Spearman's rho in B4.

1. τ<sub>b</sub> = (C − D) / √((N<sub>p</sub> − T<sub>x</sub>)(N<sub>p</sub> − T<sub>y</sub>)), where C and D are the numbers of concordant and discordant pairs, N<sub>p</sub> = n(n − 1)/2 is the number of pairs, and T<sub>x</sub> and T<sub>y</sub> are the pairs tied on the GOI and on the composite.
2. Each cell is part of n − 1 pairs. Its Kendall score is K<sub>i</sub> = (C<sub>i</sub> − D<sub>i</sub>) / (n − 1), where C<sub>i</sub> and D<sub>i</sub> are the concordant and discordant pairs involving cell i.
3. Without ties, τ<sub>b</sub> is the mean of K<sub>i</sub>. Script 02 checks this to machine precision.
4. Without ties, K<sub>i</sub> > 0.5 means that more than three quarters of the pairs of cell i are concordant. This is the K<sub>i</sub> gate of the SKCS-core.

---

## Part 1. Ranks and the probability integral transform (`scripts/01_ranks_and_PIT.r`)

The simulation has 1,950 GOI-positive cells. The composite is the mean of eight partner genes (`cepip_composite()`).

**01a (left).** Composite against GOI in expression space. **01c (right).** Copula space: u = rank(GOI)/(n + 1) and v = rank(composite)/(n + 1). Both margins are uniform on (0, 1).

<img src="figures/01a_expression_space.png" width="49%"> <img src="figures/01c_copula_space.png" width="49%">

**01b.** Margins of the GOI and the composite, with a kernel density curve reflected at 0. Both are right-skewed.

<img src="figures/01b_expression_margins.png" width="80%">

**01d (left).** Taking the values back from the log1p scale to the CP10K scale changes expression space but leaves copula space identical point for point, because a strictly increasing transformation keeps the ranks. **01e (right).** When 30% of cells are set to GOI = 0, all tied zeros receive the same average rank and form a vertical stripe at u = ((n<sub>0</sub> + 1)/2)/(n + 1). This is why the inflection point is located on GOI-positive cells only.

<img src="figures/01d_invariance.png" width="49%"> <img src="figures/01e_ties_zero_GOI.png" width="49%">

## Part 2. Spearman's rho, Kendall's tau and the per-cell scores (`scripts/02_spearman_kendall_scores.r`)

**02a (left).** Spearman's ρ is Pearson's r of u and v. **02b (right).** S<sub>i</sub> of each cell in copula space (`cepip_skcs()`). S<sub>i</sub> is positive in the lower-left and upper-right quadrants and near 0 at the middle GOI rank. Colours span the 5th to 95th percentile.

<img src="figures/02a_rho_as_pearson_of_ranks.png" width="49%"> <img src="figures/02b_Si_in_copula_space.png" width="49%">

**02c (left).** K<sub>i</sub> of each cell in copula space, the share of concordant minus discordant pairs of the cell. **02e (right).** Pitcher plot (`plot_pitcher()`). The SKCS-core holds 1,287 of 1,950 cells (66.0%).

<img src="figures/02c_Ki_in_copula_space.png" width="49%"> <img src="figures/02e_pitcher_plot.png" width="49%">

**02d.** Distributions of S<sub>i</sub> and K<sub>i</sub>, with the SKCS-core gates S<sub>i</sub> = 0.05 and K<sub>i</sub> = 0.5.

<img src="figures/02d_Si_Ki_distributions.png" width="80%">

## Part 3. GAM fit, TREP cells, the inflection point and IPRS (`scripts/03_gam_trep_ip_iprs.r`)

**03a.** GAM plot of the composite and of a housekeeping control composite, for all GOI-positive cells and for the SKCS-core. Restricting to the SKCS-core raises the deviance explained (DE) from 73.6% to 94.1%. The control composite has a DE close to 0.

![](figures/03a_GAM_plot.png)

**03b (left).** The basis dimension k is chosen by the lowest penalized residual sum of squares (PRSS) (`cepip_fit()`, `plot_prss()`). **03d (right).** The calculus definition applied to the fitted curve: f(x), f'(x) and f''(x). Here f''(x) = 0 with a change of sign at x = 3.303, while the CEP-IP inflection point from R(x) is IP<sub>x</sub> = 2.739.

<img src="figures/03b_PRSS_k_search.png" width="49%"> <img src="figures/03d_calculus_on_fit.png" width="49%">

**03c.** CEP-IP plot (`plot_fit()`). TREP cells are the top round(n × DE) cells by explanatory power (`cepip_ep()`).

![](figures/03c_CEPIP_plot.png)

**03e (left).** Residual curve R(x) of the TREP cells (`cepip_ip()`, `plot_residual()`). IP<sub>x</sub> is the first negative-to-positive zero crossing of R(x): 2.739 (95% CI 2.680 to 2.795), IPRS = 0.705 (strong tier). **03f (right).** IP<sub>x</sub> and its 95% CI as three roots: R(x) = 0, R(x) + 1.96 SE = 0 (lower bound) and R(x) − 1.96 SE = 0 (upper bound). The CI bounds are found with Brent's method.

<img src="figures/03e_residual_curve_Rx.png" width="49%"> <img src="figures/03f_brent_CI.png" width="49%">

**03g (left).** IP<sub>y</sub> = f(IP<sub>x</sub>) = 1.841 is the threshold on the composite axis that splits pre-IP from post-IP cells (`plot_fit()`). **03h (right).** The five IPRS components for the simulated module (IPRS 0.705, strong tier) and for a simulated control arm with almost no relationship (IPRS 0.252, weak tier). The guardrails check that the components computed by hand equal those of `cepip_ip()`.

<img src="figures/03g_IPx_to_IPy.png" width="49%"> <img src="figures/03h_IPRS_components.png" width="49%">

## Part 4. Pre-IP and post-IP groups in expression, score and copula space (`scripts/04_groups_and_copula.r`)

**04a (left).** CEP-IP plot of the four SKCS-core groups. SKCS-core cells below IP<sub>y</sub> are pre-IP and the others post-IP. Each group is split into TREP and non-TREP cells by the GAM fit of the SKCS-core. **04b (right).** Pitcher plot by group (`plot_pitcher()`). The two spines of the SKCS-core are the pre-IP and post-IP groups.

<img src="figures/04a_CEPIP_groups.png" width="49%"> <img src="figures/04b_pitcher_groups.png" width="49%">

**04c.** One cell (cell A) and IP<sub>y</sub> carried from expression space to copula space. IP<sub>v</sub> = #{composite ≤ IP<sub>y</sub>}/(n + 1) = 0.554.

![](figures/04c_expression_to_copula.png)

**04d.** The same cell and the SKCS-core gates carried from score space to copula space. S<sub>i</sub> = 0.05 and K<sub>i</sub> = 0.5 become u<sub>gate</sub> = 0.299 and v<sub>gate</sub> = 0.146.

![](figures/04d_scores_to_copula.png)

**04e.** Pitcher plot in copula space (v<sub>SK</sub> against u<sub>SK</sub>). The SKCS-core is the quadrant u<sub>SK</sub> > u<sub>gate</sub> and v<sub>SK</sub> > v<sub>gate</sub>.

<p align="center"><img src="figures/04e_pitcher_copula_space.png" width="55%"></p>

---

## From the calculus inflection point to the CEP-IP inflection point

- The calculus definition in Part A4 was the starting point.
- In CEP-IP, the inflection point IP<sub>x</sub> is the first negative-to-positive zero crossing of the residual curve R(x). The bounds of its 95% CI are found by root-finding with Brent's method.
- The f''(x) = 0 definition on the fitted mean curve was tested as an alternative and was not used. On the SEA-AD microglia data, the population mean curve had no interior inflection, and its steepest segment lay in the upper tail.
- Two parts of Part A carried over: the sign analysis around a root (A4), and the check that the curve is continuous and differentiable where a point is read off (A2).
- Figure 03d shows the two definitions on the simulated cells. f''(x) changes sign at x = 3.303, and R(x) crosses zero at IP<sub>x</sub> = 2.739.
