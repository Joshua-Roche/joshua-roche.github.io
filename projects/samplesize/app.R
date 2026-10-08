# Sample size when scans are read with error
#
# Shiny app accompanying "Sample size for trials whose endpoints are read from
# scans: correcting for false and missed progression calls in objective
# response, disease control and progression-free survival".
#
# Runs with ordinary Shiny, or in the browser with Shinylive (see README.md).
# Uses only shiny, bslib and base R so that it loads quickly under webR.

library(shiny)
library(bslib)

source("misclass_core.R")

# Edit these two lines to point readers to the paper and the code.
PAPER_CITATION <- paste(
  "[Authors]. Sample size for trials whose endpoints are read from scans:",
  "correcting for false and missed progression calls in objective response,",
  "disease control and progression-free survival. Submitted to Clinical Trials."
)
CODE_URL <- ""   # e.g. "https://github.com/<user>/<repo>"

# ---------------------------------------------------------------------------
# Colours (validated categorical slots 1-4, neutral reference, sequential ramp)
# ---------------------------------------------------------------------------

COL <- list(
  corrected = "#2a78d6", uncorrected = "#eb6834", misspec = "#1baf7a",
  smaller = "#eda100", perfect = "#2b2a28",
  text1 = "#0b0b0b", text2 = "#52514e", axis = "#9a9893", grid = "#e7e6e2"
)
SEQ <- c("#104281", "#256abf", "#5598e7", "#86b6ef")   # Se = 1, 0.9, 0.8, 0.7
LTY <- c(perfect = 2, uncorrected = 1, corrected = 1, misspec = 4, smaller = 3)
LWD <- c(perfect = 1.6, uncorrected = 2.4, corrected = 3.2, misspec = 2.4, smaller = 2.2)

EXAMPLES <- list(
  M = list(m0 = 7.2, hr = 0.70, pi_d = 10, mu_med = 3.4, int1 = 6, switch_wk = 52,
           int2 = 12, tau = 24, alpha = 0.05, sided = "2", power = 0.80, ratio = 1,
           se = 1, sp = 0.90, se_r = 1, fr = 0.05, orr0 = 0.30, orr1 = 0.50,
           hr_orr = 1, hr_small = 0.85, orr1_small = 0.42, dcr_week = "24"),
  B = list(m0 = 9, hr = 0.69, pi_d = 5, mu_med = 24, int1 = 12, switch_wk = 156,
           int2 = 12, tau = 36, alpha = 0.025, sided = "1", power = 0.90, ratio = 2,
           se = 1, sp = 0.98, se_r = 1, fr = 0.025, orr0 = 0.30, orr1 = 0.50,
           hr_orr = 1, hr_small = 0.80, orr1_small = 0.42, dcr_week = "24")
)

fmt_n <- function(x) formatC(round(x), format = "d", big.mark = ",")
fmt_p <- function(x, d = 0) paste0(formatC(100 * x, format = "f", digits = d), "%")
fmt_x <- function(x, d = 2) formatC(x, format = "f", digits = d)
fmt_acc <- function(x) format(round(x, 3), nsmall = 2)
strip_tags <- function(x) gsub("<[^>]+>", "", x)
SER <- "Se<sub>R</sub>"

# ---------------------------------------------------------------------------
# Planning for each endpoint
# ---------------------------------------------------------------------------

design_common <- function(s) {
  list(a = s$alpha, sd = as.numeric(s$sided), pw = s$power, w = s$ratio / (1 + s$ratio),
       z2 = (qnorm(1 - s$alpha / as.numeric(s$sided)) + qnorm(s$power))^2)
}

plan_pfs <- function(s, sch, G) {
  d <- design_common(s)
  arms <- arms_from_median(s$m0, s$hr, s$pi_d / 100)
  el_p <- expected_logrank(arms[[1]], arms[[2]], G, 1, 1, d$w)
  el_r <- expected_logrank(arms[[1]], arms[[2]], G, s$se, s$sp, d$w)
  nr_r <- d$z2 * el_r$v / el_r$u^2
  nr_p <- d$z2 * el_p$v / el_p$u^2
  naive <- n_pfs_naive(arms[[1]], arms[[2]], s$tau, d$a, d$pw, d$sd, d$w)
  calc <- n_pfs_calculator(arms[[1]], arms[[2]], s$tau, s$sp, sch$visits[1],
                           d$a, d$pw, d$sd, d$w)
  n <- c(uncorrected = unname(naive["n"]), calculator = unname(calc["n"]),
         corrected = ceiling(nr_r), perfect = ceiling(nr_p))
  ev <- c(unname(naive["events"]), unname(calc["events"]), nr_r * el_r$ev, nr_p * el_p$ev)
  pow <- vapply(n, function(k) power_logrank(k, el_r, d$a, d$sd), 0)
  out <- list(n = n, events = ev, power = pow, d = d, arms = arms,
              f = fp_hazard(s$sp, sch$visits[1]),
              hr_rec = unname(calc["hr_obs"]))
  if (isTRUE(s$misspec)) {
    el_t <- expected_logrank(arms[[1]], arms[[2]], G, s$se_t, s$sp_t, d$w)
    out$power_true <- vapply(n, function(k) power_logrank(k, el_t, d$a, d$sd), 0)
  }
  if (isTRUE(s$smaller)) {
    a2 <- arms_from_median(s$m0, s$hr_small, s$pi_d / 100)
    el_s <- expected_logrank(a2[[1]], a2[[2]], G, 1, 1, d$w)
    out$n_small <- ceiling(d$z2 * el_s$v / el_s$u^2)
  }
  out
}

dcr_kappa <- function(s, sch) {
  k <- match(as.numeric(s$dcr_week), round(sch$weeks, 6))
  if (is.na(k)) k <- which.min(abs(sch$weeks - 24))
  k
}

plan_dcr <- function(s, sch) {
  d <- design_common(s)
  kp <- dcr_kappa(s, sch)
  arms <- arms_from_median(s$m0, s$hr, s$pi_d / 100)
  th <- c(dcr_true(arms[[1]], sch, kp), dcr_true(arms[[2]], sch, kp))
  pr <- c(dcr_observed(arms[[1]], sch, kp, s$se, s$sp), dcr_observed(arms[[2]], sch, kp, s$se, s$sp))
  n <- c(uncorrected = n_two_prop(th[1], th[2], d$a, d$pw, d$sd, d$w),
         corrected = n_two_prop(pr[1], pr[2], d$a, d$pw, d$sd, d$w))
  n["perfect"] <- n[["uncorrected"]]
  pow <- vapply(n, function(k) power_two_prop(k, pr[1], pr[2], d$a, d$sd, d$w), 0)
  out <- list(n = n, power = pow, d = d, arms = arms, kappa = kp, true = th, rec = pr)
  if (isTRUE(s$misspec)) {
    pt <- c(dcr_observed(arms[[1]], sch, kp, s$se_t, s$sp_t), dcr_observed(arms[[2]], sch, kp, s$se_t, s$sp_t))
    out$power_true <- vapply(n, function(k) power_two_prop(k, pt[1], pt[2], d$a, d$sd, d$w), 0)
  }
  if (isTRUE(s$smaller)) {
    a2 <- arms_from_median(s$m0, s$hr_small, s$pi_d / 100)
    out$n_small <- n_two_prop(th[1], dcr_true(a2[[2]], sch, kp), d$a, d$pw, d$sd, d$w)
  }
  out
}

orr_arm <- function(theta, lam, pi_d, sch) {
  rho <- theta / exp(-lam * sch$visits[1])
  make_arm((1 - pi_d) * lam, pi_d * lam, rho)
}

plan_orr <- function(s, sch) {
  d <- design_common(s)
  lam0 <- log(2) / s$m0; lam1 <- lam0 * s$hr_orr; pd <- s$pi_d / 100
  a0 <- orr_arm(s$orr0, lam0, pd, sch); a1 <- orr_arm(s$orr1, lam1, pd, sch)
  M <- reader_matrix(1, 1, s$se_r, s$fr)
  pr <- c(orr_observed(a0, sch, M, s$confirmed), orr_observed(a1, sch, M, s$confirmed))
  n <- c(uncorrected = n_two_prop(s$orr0, s$orr1, d$a, d$pw, d$sd, d$w),
         corrected = n_two_prop(pr[1], pr[2], d$a, d$pw, d$sd, d$w))
  n["perfect"] <- n[["uncorrected"]]
  pow <- vapply(n, function(k) power_two_prop(k, pr[1], pr[2], d$a, d$sd, d$w), 0)
  out <- list(n = n, power = pow, d = d, a0 = a0, lam1 = lam1, M = M,
              true = c(s$orr0, s$orr1), rec = pr,
              rho = c(a0$rho, a1$rho))
  if (isTRUE(s$misspec)) {
    Mt <- reader_matrix(1, 1, s$se_r_t, s$fr_t)
    pt <- c(orr_observed(a0, sch, Mt, s$confirmed), orr_observed(a1, sch, Mt, s$confirmed))
    out$power_true <- vapply(n, function(k) power_two_prop(k, pt[1], pt[2], d$a, d$sd, d$w), 0)
    out$Mt <- Mt
  }
  if (isTRUE(s$smaller))
    out$n_small <- n_two_prop(s$orr0, s$orr1_small, d$a, d$pw, d$sd, d$w)
  out
}

# ---------------------------------------------------------------------------
# Power curves: power against the true effect for each design
# ---------------------------------------------------------------------------

curves_pfs <- function(s, P, G) {
  lower <- max(0.3, min(0.5, s$hr - 0.2))
  x <- round(seq(1, lower, by = -0.01), 2)
  d <- P$d; a0 <- P$arms[[1]]; pd <- s$pi_d / 100
  readers <- list(perfect = c(1, 1), reader = c(s$se, s$sp))
  if (isTRUE(s$misspec)) readers$true <- c(s$se_t, s$sp_t)
  els <- lapply(readers, function(r) {
    S0 <- pfs_survival_grid(a0, G, r[1], r[2])
    lapply(x, function(h) expected_logrank(a0, arms_from_median(s$m0, h, pd)[[2]], G, r[1], r[2], d$w, S0))
  })
  pw <- function(n, who) vapply(els[[who]], function(el) power_logrank(n, el, d$a, d$sd), 0)
  y <- list(perfect = pw(P$n[["perfect"]], "perfect"),
            uncorrected = pw(P$n[["uncorrected"]], "reader"),
            corrected = pw(P$n[["corrected"]], "reader"))
  if (isTRUE(s$misspec)) y$misspec <- pw(P$n[["corrected"]], "true")
  if (isTRUE(s$smaller)) y$smaller <- pw(P$n_small, "perfect")
  list(x = x, y = y, xlim = c(1, lower), design_x = s$hr, xlab = "True hazard ratio",
       xname = "True HR", xfmt = function(v) fmt_x(v, 2))
}

curves_dcr <- function(s, P, sch) {
  lower <- max(0.3, min(0.5, s$hr - 0.2))
  x <- round(seq(1, lower, by = -0.01), 2)
  d <- P$d; kp <- P$kappa; a0 <- P$arms[[1]]; pd <- s$pi_d / 100
  rec <- function(se, sp) {
    p0 <- dcr_observed(a0, sch, kp, se, sp)
    list(p0 = p0, p1 = vapply(x, function(h) dcr_observed(arms_from_median(s$m0, h, pd)[[2]], sch, kp, se, sp), 0))
  }
  rp <- rec(1, 1); rr <- rec(s$se, s$sp)
  pw <- function(n, r) power_two_prop(n, r$p0, r$p1, d$a, d$sd, d$w)
  y <- list(perfect = pw(P$n[["perfect"]], rp),
            uncorrected = pw(P$n[["uncorrected"]], rr),
            corrected = pw(P$n[["corrected"]], rr))
  if (isTRUE(s$misspec)) y$misspec <- pw(P$n[["corrected"]], rec(s$se_t, s$sp_t))
  if (isTRUE(s$smaller)) y$smaller <- pw(P$n_small, rp)
  list(x = x, y = y, xlim = c(1, lower), design_x = s$hr,
       xlab = "True hazard ratio (sets the true DCR)",
       xname = "True HR", xfmt = function(v) fmt_x(v, 2))
}

curves_orr <- function(s, P, sch) {
  upper <- min(0.95, s$orr0 + 1.6 * (s$orr1 - s$orr0))
  x <- round(seq(s$orr0, upper, by = 0.01), 2)
  pd <- s$pi_d / 100
  rho1 <- x / exp(-P$lam1 * sch$visits[1])
  x <- x[rho1 <= 1]
  d <- P$d
  # The recorded ORR is linear in the latent responder fraction rho, so the
  # recursion is run once per reader for responders and non-responders.
  rho1 <- x / exp(-P$lam1 * sch$visits[1])
  rec <- function(M) {
    a1 <- orr_arm(0, P$lam1, pd, sch)
    pR <- orr_observed(modifyList(a1, list(rho = 1)), sch, M, s$confirmed)
    pS <- orr_observed(modifyList(a1, list(rho = 0)), sch, M, s$confirmed)
    list(p0 = orr_observed(P$a0, sch, M, s$confirmed), p1 = rho1 * pR + (1 - rho1) * pS)
  }
  rp <- rec(reader_matrix()); rr <- rec(P$M)
  pw <- function(n, r) power_two_prop(n, r$p0, r$p1, d$a, d$sd, d$w)
  y <- list(perfect = pw(P$n[["perfect"]], rp),
            uncorrected = pw(P$n[["uncorrected"]], rr),
            corrected = pw(P$n[["corrected"]], rr))
  if (isTRUE(s$misspec)) y$misspec <- pw(P$n[["corrected"]], rec(P$Mt))
  if (isTRUE(s$smaller)) y$smaller <- pw(P$n_small, rp)
  list(x = x, y = y, xlim = range(x), design_x = s$orr1,
       xlab = "True ORR in the experimental arm",
       xname = "True ORR", xfmt = function(v) fmt_p(v))
}

series_labels <- function(s, P) {
  acc <- if (s$endpoint == "ORR") sprintf("%s %s, FR %s", SER, fmt_acc(s$se_r), fmt_acc(s$fr))
         else sprintf("Se %s, Sp %s", fmt_acc(s$se), fmt_acc(s$sp))
  acc_t <- if (s$endpoint == "ORR") sprintf("%s %s, FR %s", SER, fmt_acc(s$se_r_t), fmt_acc(s$fr_t))
           else sprintf("Se %s, Sp %s", fmt_acc(s$se_t), fmt_acc(s$sp_t))
  small <- if (s$endpoint == "ORR") sprintf("ORR %s", fmt_p(s$orr1_small)) else sprintf("HR %s", fmt_x(s$hr_small))
  list(
    perfect = sprintf("Perfect reading (n = %s)", fmt_n(P$n[["perfect"]])),
    uncorrected = sprintf("Uncorrected size, read at %s (n = %s)", acc, fmt_n(P$n[["uncorrected"]])),
    corrected = sprintf("Corrected size, read at %s (n = %s)", acc, fmt_n(P$n[["corrected"]])),
    misspec = sprintf("Corrected size, but true accuracy %s", acc_t),
    smaller = if (!is.null(P$n_small)) sprintf("Perfect reading, sized for %s (n = %s)", small, fmt_n(P$n_small)) else ""
  )
}

ORDER <- c("smaller", "misspec", "uncorrected", "corrected", "perfect")

draw_power <- function(cv, target, alpha_line) {
  par(mar = c(4.2, 4.4, 0.6, 1.4), mgp = c(2.7, 0.55, 0), tcl = -0.25, las = 1,
      col.axis = COL$text2, col.lab = COL$text1, cex.axis = 0.9, cex.lab = 0.95,
      bg = "white", family = "sans")
  plot(NA, xlim = cv$xlim, ylim = c(0, 1), xaxs = "i", yaxs = "i", axes = FALSE,
       xlab = cv$xlab, ylab = "Probability of rejecting the null")
  abline(h = seq(0.2, 1, 0.2), col = COL$grid, lwd = 1)
  xt <- pretty(cv$xlim, n = 7); xt <- xt[xt >= min(cv$xlim) - 1e-9 & xt <= max(cv$xlim) + 1e-9]
  axis(1, at = xt, labels = cv$xfmt(xt), col = COL$axis, lwd = 1, lwd.ticks = 1)
  axis(2, at = seq(0, 1, 0.2), labels = paste0(seq(0, 100, 20), "%"), lwd = 0, lwd.ticks = 0)
  segments(cv$design_x, 0, cv$design_x, 1, lty = 3, col = COL$axis)
  abline(h = target, lty = 2, col = COL$axis)
  text(cv$xlim[1], target, sprintf("target %s", fmt_p(target)), adj = c(-0.08, -0.5),
       cex = 0.78, col = COL$text2)
  text(cv$design_x, 0.035, "design", adj = c(-0.12, 0), cex = 0.78, col = COL$text2)
  for (k in ORDER) if (!is.null(cv$y[[k]]))
    lines(cv$x, cv$y[[k]], col = COL[[k]], lwd = LWD[[k]], lty = LTY[[k]], lend = 1)
}

legend_ui <- function(keys, labels) {
  sw <- function(k) {
    dash <- switch(as.character(LTY[[k]]), "1" = "", "2" = "6 4", "3" = "2 3", "4" = "8 3 2 3")
    HTML(sprintf('<svg width="30" height="10" aria-hidden="true"><line x1="1" y1="5" x2="29" y2="5" stroke="%s" stroke-width="%s" stroke-dasharray="%s"/></svg>',
                 COL[[k]], max(2, LWD[[k]]), dash))
  }
  div(class = "legend-row", lapply(keys, function(k) span(class = "legend-item", sw(k), HTML(labels[[k]]))))
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

num <- function(id, label, value, step = NULL, min = NA, max = NA)
  numericInput(id, label, value, min = min, max = max, step = step)

# The default bslib theme ships precompiled, so it loads fast under webR;
# fonts and a few colours are adjusted in CSS instead of through bs_theme().

css <- "
:root, [data-bs-theme=light] {
  --bs-body-font-family: system-ui, -apple-system, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
  --bs-body-font-size: .9375rem;
}
body { background: #f6f6f4; }
.bslib-sidebar-layout > .sidebar { background: #fbfbfa; }
.sidebar .accordion-button { font-weight: 600; font-size: .9rem; padding: .6rem .9rem; }
.sidebar .form-label, .sidebar label { font-size: .84rem; color: #3a3936; }
.sidebar .form-group { margin-bottom: .55rem; }
.sidebar .help-block, .sidebar .help { font-size: .78rem; color: #6b6a66; margin-top: -.2rem; }
.examples { font-size: .84rem; color: #52514e; margin: -.25rem 0 .4rem; }
.lede { color: #52514e; font-size: .92rem; max-width: 62rem; margin: 0 0 .25rem; }
.stat { background: #fff; border: 1px solid #e3e2de; border-radius: .6rem; padding: .8rem 1rem;
        height: 100%; border-top-width: 3px; }
.stat.unc { border-top-color: #eb6834; } .stat.cor { border-top-color: #2a78d6; }
.stat.rec { border-top-color: #9a9893; }
.stat .k { font-size: .76rem; letter-spacing: .04em; text-transform: uppercase; color: #52514e; }
.stat .v { font-size: 1.85rem; font-weight: 650; line-height: 1.15; color: #0b0b0b;
           font-variant-numeric: tabular-nums; margin: .15rem 0; }
.stat .v small { font-size: .95rem; font-weight: 500; color: #52514e; }
.stat .s { font-size: .82rem; color: #52514e; }
.legend-row { display: flex; flex-wrap: wrap; gap: .3rem 1.2rem; font-size: .84rem;
              color: #3a3936; margin: .2rem 0 .1rem .4rem; }
.legend-item { display: inline-flex; align-items: center; gap: .4rem; }
.readout { font-size: .84rem; color: #52514e; min-height: 1.6rem; margin: .35rem 0 0 .4rem;
           font-variant-numeric: tabular-nums; }
.readout b { color: #0b0b0b; font-weight: 600; }
.readout .chip { display: inline-flex; align-items: center; gap: .3rem; margin-right: 1rem; }
.readout .dot { width: .6rem; height: .6rem; border-radius: 50%; display: inline-block; }
.note { font-size: .82rem; color: #52514e; max-width: 62rem; }
.shiny-table { font-variant-numeric: tabular-nums; font-size: .88rem; }
.shiny-table td, .shiny-table th { padding: .35rem .7rem !important; }
.about { max-width: 52rem; font-size: .92rem; }
.about h5 { margin-top: 1.1rem; font-size: 1rem; }
.shiny-output-error-validation { color: #8a3b12; font-size: .9rem; }
"

sidebar_ui <- sidebar(
  width = 350,
  radioButtons("endpoint", "Endpoint",
               c("PFS" = "PFS", "Disease control (DCR)" = "DCR", "Objective response (ORR)" = "ORR"),
               selected = "PFS"),
  div(class = "examples", "Load an example: ",
      actionLink("ex_M", "mesothelioma"), " · ", actionLink("ex_B", "PALOMA-2")),
  accordion(
    multiple = TRUE, open = c("Trial design", "Reader accuracy, per scan"),
    accordion_panel(
      "Trial design",
      num("m0", "Control median PFS (months)", 7.2, 0.1, 0.5),
      conditionalPanel("input.endpoint != 'ORR'",
                       num("hr", "True hazard ratio to detect", 0.70, 0.01, 0.05, 0.99)),
      conditionalPanel("input.endpoint == 'DCR'",
                       selectInput("dcr_week", "Disease control assessed at the scan in week",
                                   choices = c(6, 12, 18, 24, 30, 36, 42, 48, 60, 72, 84, 96),
                                   selected = 24)),
      conditionalPanel(
        "input.endpoint == 'ORR'",
        layout_columns(col_widths = c(6, 6), gap = ".5rem",
                       num("orr0", "Control ORR", 0.30, 0.01, 0.01, 0.99),
                       num("orr1", "Experimental ORR", 0.50, 0.01, 0.01, 0.99)),
        num("hr_orr", "Hazard ratio for progression or death", 1, 0.01, 0.05, 2),
        div(class = "help", "Below 1, patients on the new drug attend more scans, and so have more chances of a false response."),
        checkboxInput("confirmed", "Response must be confirmed (two consecutive PR reads)", FALSE)
      ),
      layout_columns(col_widths = c(6, 6), gap = ".5rem",
                     num("alpha", "Significance level α", 0.05, 0.005, 0.001, 0.5),
                     selectInput("sided", "Test", c("Two-sided" = "2", "One-sided" = "1"))),
      layout_columns(col_widths = c(6, 6), gap = ".5rem",
                     num("power", "Power", 0.80, 0.05, 0.5, 0.99),
                     num("ratio", "Ratio, new : control", 1, 0.5, 0.2, 5))
    ),
    accordion_panel(
      "Reader accuracy, per scan",
      conditionalPanel(
        "input.endpoint != 'ORR'",
        sliderInput("sp", "Specificity (Sp): a progression-free scan is not called PD",
                    min = 0.80, max = 1, value = 0.90, step = 0.005, ticks = FALSE),
        sliderInput("se", "Sensitivity (Se): a scan after progression is called PD",
                    min = 0.50, max = 1, value = 1, step = 0.01, ticks = FALSE)
      ),
      conditionalPanel(
        "input.endpoint == 'ORR'",
        sliderInput("fr", "False response (FR): stable disease read as PR",
                    min = 0, max = 0.20, value = 0.05, step = 0.005, ticks = FALSE),
        sliderInput("se_r", HTML("Response sensitivity (Se<sub>R</sub>): a partial response read as PR"),
                    min = 0.50, max = 1, value = 1, step = 0.01, ticks = FALSE),
        div(class = "help", "Progression is assumed to be read correctly for ORR.")
      ),
      checkboxInput("misspec", "Show the corrected trial if the true accuracy is different", FALSE),
      conditionalPanel(
        "input.misspec",
        conditionalPanel(
          "input.endpoint != 'ORR'",
          sliderInput("sp_t", "True Sp", min = 0.80, max = 1, value = 1, step = 0.005, ticks = FALSE),
          sliderInput("se_t", "True Se", min = 0.50, max = 1, value = 1, step = 0.01, ticks = FALSE)
        ),
        conditionalPanel(
          "input.endpoint == 'ORR'",
          sliderInput("fr_t", "True FR", min = 0, max = 0.20, value = 0, step = 0.005, ticks = FALSE),
          sliderInput("se_r_t", HTML("True Se<sub>R</sub>"), min = 0.50, max = 1, value = 1, step = 0.01, ticks = FALSE)
        )
      )
    ),
    accordion_panel(
      "Scan schedule and follow-up",
      layout_columns(col_widths = c(6, 6), gap = ".5rem",
                     num("int1", "Scan every (weeks)", 6, 1, 1, 52),
                     num("switch_wk", "up to week", 52, 1, 1, 520)),
      layout_columns(col_widths = c(6, 6), gap = ".5rem",
                     num("int2", "then every (weeks)", 12, 1, 0, 104),
                     num("tau", "Follow-up (months)", 24, 1, 1, 240)),
      uiOutput("sched_text", class = "help"),
      div(class = "help", "Every patient is followed for the same time; there is no staggered entry.")
    ),
    accordion_panel(
      "Disease course",
      num("pi_d", "Deaths before progression (% of PFS events)", 10, 1, 0, 99),
      num("mu_med", "Median survival after an undetected progression (months)", 3.4, 0.1, 0.1)
    ),
    accordion_panel(
      "Comparison",
      checkboxInput("smaller", "Add a perfectly read trial sized for a smaller effect", TRUE),
      conditionalPanel(
        "input.smaller",
        conditionalPanel("input.endpoint != 'ORR'",
                         num("hr_small", "Smaller effect: hazard ratio", 0.85, 0.01, 0.05, 0.99)),
        conditionalPanel("input.endpoint == 'ORR'",
                         num("orr1_small", "Smaller effect: experimental ORR", 0.42, 0.01, 0.01, 0.99))
      )
    )
  )
)

about_ui <- div(
  class = "about",
  h5("What this calculates"),
  p("Objective response, disease control and progression-free survival are recorded from scans, and readers make mistakes. ",
    "A false progression call ends a patient's PFS early; a missed progression delays it; a false or missed partial response changes the recorded ORR. ",
    "Standard sample-size formulas assume every scan is read correctly, so a trial sized that way is under-powered for the effect it was designed to detect."),
  p("The corrected sample size is the standard one, worked out for the endpoint as it will actually be recorded. ",
    "The planned analysis does not change: a log-rank test or a comparison of two proportions on the recorded data."),
  h5("Methods"),
  tags$ul(
    tags$li(tags$b("Uncorrected:"), " Schoenfeld's formula (PFS) or the two-proportion formula on the true rates, as if every scan were read correctly."),
    tags$li(tags$b("Calculator (PFS):"), " false calls add a hazard f = −ln(Sp)/Δ to both arms, so the recorded hazard ratio is (HR + f/λ₀)/(1 + f/λ₀). Δ is the first scan interval. Missed progressions are ignored."),
    tags$li(tags$b("Exact (recommended):"), " the recorded PFS curves, DCR or ORR are worked out exactly under the model, and the expected log-rank statistic or the two-proportion formula gives the sample size. Power curves use the same calculation.")
  ),
  h5("Model assumptions"),
  tags$ul(
    tags$li("Constant hazards of progression and of death before progression; a patient whose progression is missed dies at a constant rate."),
    tags$li("Each scan is read independently given the true state, with the same accuracy in both arms. Differential reading, for example unblinded review in an open-label trial, is not covered."),
    tags$li("A PD read ends assessment; deaths are recorded exactly; every patient keeps the same scan schedule and is followed for the same time."),
    tags$li(HTML("For ORR, progression is read correctly, a true PR is read as PR with probability Se<sub>R</sub> and otherwise as SD, and stable disease is read as PR with probability FR."))
  ),
  h5("Reading the power curve"),
  p("When the assumed accuracy is right, the corrected trial has almost the same power as a perfectly read trial at every true effect, so the correction does not make small effects easier to detect. ",
    "A trial sized for a smaller effect does. If the true accuracy is better than assumed, the corrected trial is larger than it needed to be; if worse, it is still under-powered. ",
    "Reader accuracy can be estimated during the trial, blind to treatment, and the sample size re-estimated."),
  h5("Reference"),
  p(PAPER_CITATION),
  if (nzchar(CODE_URL)) p("Code: ", a(href = CODE_URL, CODE_URL, target = "_blank")),
  p(class = "note", "This tool supports planning; it is not a substitute for a statistician's review of a specific trial.")
)

ui <- page_sidebar(
  title = "Sample size when scans are read with error",
  window_title = "Sample size with misread scans",
  sidebar = sidebar_ui,
  tags$head(tags$style(HTML(css))),
  p(class = "lede",
    "How many patients a trial needs when progression and response are read from scans with error, ",
    "and how its power changes when the sample size is, or is not, corrected."),
  layout_columns(
    col_widths = c(4, 4, 4), fill = FALSE, gap = ".75rem",
    uiOutput("vb_unc"), uiOutput("vb_cor"), uiOutput("vb_rec")
  ),
  navset_card_underline(
    id = "tabs", full_screen = FALSE,
    nav_panel(
      "Power curve",
      plotOutput("pc", height = "420px",
                 hover = hoverOpts("pc_hover", delay = 60, delayType = "throttle", nullOutside = TRUE)),
      uiOutput("pc_legend"),
      uiOutput("pc_readout"),
      accordion(open = FALSE,
                accordion_panel("Values in a table", tableOutput("pc_table")))
    ),
    nav_panel(
      "Sample sizes",
      tableOutput("ss_table"),
      uiOutput("ss_notes")
    ),
    nav_panel(
      "Reader accuracy",
      plotOutput("sens", height = "400px"),
      uiOutput("sens_legend"),
      uiOutput("sens_note")
    ),
    nav_panel("About", about_ui)
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  load_example <- function(key) {
    e <- EXAMPLES[[key]]
    for (id in c("m0", "hr", "pi_d", "mu_med", "int1", "switch_wk", "int2", "tau",
                 "alpha", "power", "ratio", "orr0", "orr1", "hr_orr", "hr_small", "orr1_small"))
      updateNumericInput(session, id, value = e[[id]])
    for (id in c("se", "sp", "se_r", "fr")) updateSliderInput(session, id, value = e[[id]])
    if (identical(input$endpoint, "DCR")) {           # the paper's DCR design: HR 0.6 at week 24
      updateNumericInput(session, "hr", value = 0.6)
      updateNumericInput(session, "hr_small", value = 0.75)
      updateSelectInput(session, "dcr_week", selected = "24")
    }
    updateSelectInput(session, "sided", selected = e$sided)
  }
  observeEvent(input$ex_M, load_example("M"))
  observeEvent(input$ex_B, load_example("B"))

  raw <- reactive(list(
    endpoint = input$endpoint, m0 = input$m0, hr = input$hr, dcr_week = input$dcr_week,
    orr0 = input$orr0, orr1 = input$orr1, hr_orr = input$hr_orr, confirmed = isTRUE(input$confirmed),
    alpha = input$alpha, sided = input$sided, power = input$power, ratio = input$ratio,
    se = input$se, sp = input$sp, se_r = input$se_r, fr = input$fr,
    misspec = isTRUE(input$misspec), se_t = input$se_t, sp_t = input$sp_t,
    se_r_t = input$se_r_t, fr_t = input$fr_t,
    int1 = input$int1, switch_wk = input$switch_wk, int2 = input$int2, tau = input$tau,
    pi_d = input$pi_d, mu_med = input$mu_med,
    smaller = isTRUE(input$smaller), hr_small = input$hr_small, orr1_small = input$orr1_small
  ))
  par_in <- debounce(raw, 300)

  ok <- function(x, lo, hi) is.numeric(x) && length(x) == 1 && !is.na(x) && x >= lo && x <= hi

  sch <- reactive({
    s <- par_in()
    validate(
      need(ok(s$int1, 1, 104), "Scan interval must be between 1 and 104 weeks."),
      need(ok(s$switch_wk, s$int1, 1e4), "The first scan interval must end at or after the first scan."),
      need(ok(s$int2, 0, 104), "The second scan interval must be between 0 (no further scans) and 104 weeks."),
      need(ok(s$tau, 0.5, 240), "Follow-up must be between 0.5 and 240 months."),
      need(s$int1 * WEEK < s$tau, "Follow-up must be longer than the first scan interval."),
      need(ok(s$mu_med, 0.05, 600), "Median survival after an undetected progression must be positive.")
    )
    make_schedule(s$int1, s$switch_wk, s$int2, s$tau, log(2) / s$mu_med)
  })
  grid <- reactive(pfs_grid(sch()))

  output$sched_text <- renderUI({
    w <- sch()$weeks
    shown <- if (length(w) > 12) c(head(w, 6), NA, tail(w, 3)) else w
    txt <- paste(ifelse(is.na(shown), "…", formatC(shown, format = "fg")), collapse = ", ")
    sprintf("%d scans, at weeks %s.", length(w), txt)
  })

  observeEvent(sch(), {
    w <- formatC(sch()$weeks, format = "fg")
    cur <- isolate(input$dcr_week)
    sel <- if (!is.null(cur) && cur %in% w) cur else w[which.min(abs(sch()$weeks - 24))]
    updateSelectInput(session, "dcr_week", choices = w, selected = sel)
  })

  plan <- reactive({
    s <- par_in(); sc <- sch()
    validate(
      need(ok(s$m0, 0.1, 600), "Control median PFS must be positive."),
      need(ok(s$alpha, 1e-4, 0.5), "α must be between 0.0001 and 0.5."),
      need(ok(s$power, 0.5, 0.999), "Power must be between 0.5 and 0.999."),
      need(ok(s$ratio, 0.1, 10), "Allocation ratio must be between 0.1 and 10."),
      need(ok(s$pi_d, 0, 99), "Deaths before progression must be between 0% and 99% of PFS events.")
    )
    if (s$endpoint != "ORR") {
      validate(need(ok(s$hr, 0.05, 0.995), "The hazard ratio to detect must be between 0.05 and 0.995."))
      if (s$smaller) validate(need(ok(s$hr_small, s$hr + 1e-6, 0.995),
                                   "The smaller effect needs a hazard ratio between the target and 1."))
    } else {
      validate(
        need(ok(s$orr0, 0.001, 0.99) && ok(s$orr1, s$orr0 + 0.005, 0.999),
             "ORRs must lie between 0 and 1, with the experimental ORR higher."),
        need(ok(s$hr_orr, 0.05, 3), "The hazard ratio for progression or death must be between 0.05 and 3.")
      )
      if (s$smaller) validate(need(ok(s$orr1_small, s$orr0 + 0.005, s$orr1 - 1e-6),
                                   "The smaller effect needs an experimental ORR between the two ORRs."))
      lam0 <- log(2) / s$m0; v1 <- sc$visits[1]
      validate(need(s$orr0 <= exp(-lam0 * v1) && s$orr1 <= exp(-lam0 * s$hr_orr * v1),
                    sprintf("Only %s of control patients are alive and progression-free at the first scan, so these ORRs are not possible.",
                            fmt_p(exp(-lam0 * v1)))))
    }
    validate(need(length(sc$visits) >= 1, "The schedule has no scans within follow-up."))
    switch(s$endpoint,
           PFS = plan_pfs(s, sc, grid()),
           DCR = plan_dcr(s, sc),
           ORR = plan_orr(s, sc))
  })

  curves <- reactive({
    s <- par_in(); P <- plan()
    switch(s$endpoint,
           PFS = curves_pfs(s, P, grid()),
           DCR = curves_dcr(s, P, sch()),
           ORR = curves_orr(s, P, sch()))
  })

  # --- headline numbers -----------------------------------------------------

  stat <- function(cls, k, v, s) div(class = paste("stat", cls), div(class = "k", k), div(class = "v", v), div(class = "s", s))

  output$vb_unc <- renderUI({
    P <- plan(); s <- isolate(par_in())
    ev <- if (s$endpoint == "PFS") HTML(sprintf(" <small>patients · %s events</small>", fmt_n(P$events[1]))) else HTML(" <small>patients</small>")
    stat("unc", "Uncorrected size", tagList(fmt_n(P$n[["uncorrected"]]), ev),
         sprintf("Assumes every scan is read correctly. With the stated reader its power is %s, not %s.",
                 fmt_p(P$power[["uncorrected"]]), fmt_p(s$power)))
  })

  output$vb_cor <- renderUI({
    P <- plan(); s <- isolate(par_in())
    k <- if (s$endpoint == "PFS") 3 else 2
    ev <- if (s$endpoint == "PFS") HTML(sprintf(" <small>patients · %s events</small>", fmt_n(P$events[k]))) else HTML(" <small>patients</small>")
    stat("cor", "Corrected size (exact)", tagList(fmt_n(P$n[["corrected"]]), ev),
         sprintf("%s× the uncorrected size; restores %s power for the same true effect.",
                 formatC(P$n[["corrected"]] / P$n[["uncorrected"]], format = "f", digits = 2), fmt_p(s$power)))
  })

  output$vb_rec <- renderUI({
    P <- plan(); s <- isolate(par_in())
    if (s$endpoint == "PFS") {
      stat("rec", "Recorded hazard ratio", HTML(sprintf("≈ %s <small>true %s</small>", fmt_x(P$hr_rec), fmt_x(s$hr))),
           sprintf("False calls add a hazard of %s per month to both arms (calculator approximation, first scan interval).",
                   formatC(P$f, format = "f", digits = 4)))
    } else {
      lab <- if (s$endpoint == "DCR") "Recorded DCR, control vs new" else "Recorded ORR, control vs new"
      stat("rec", lab, HTML(sprintf("%s vs %s", fmt_p(P$rec[1]), fmt_p(P$rec[2]))),
           sprintf("True %s vs %s%s.", fmt_p(P$true[1]), fmt_p(P$true[2]),
                   if (s$endpoint == "DCR") sprintf(", at the week-%s scan", formatC(sch()$weeks[P$kappa], format = "fg")) else ""))
    }
  })

  # --- power curve ------------------------------------------------------------

  output$pc <- renderPlot({
    cv <- curves(); s <- isolate(par_in())
    draw_power(cv, s$power, s$alpha)
  }, res = 96)

  shown_keys <- reactive({
    cv <- curves()
    intersect(c("perfect", "corrected", "uncorrected", "misspec", "smaller"), names(cv$y))
  })

  output$pc_legend <- renderUI({
    s <- par_in(); P <- plan()
    legend_ui(shown_keys(), series_labels(s, P))
  })

  output$pc_readout <- renderUI({
    h <- input$pc_hover; cv <- curves()
    if (is.null(h) || is.null(h$x))
      return(div(class = "readout", "Point at the plot to read off the power at any true effect."))
    i <- which.min(abs(cv$x - h$x))
    nm <- c(perfect = "Perfect reading", corrected = "Corrected", uncorrected = "Uncorrected",
            misspec = "Corrected, true accuracy", smaller = "Smaller-effect design")
    div(class = "readout",
        tags$b(sprintf("%s %s:", cv$xname, cv$xfmt(cv$x[i]))), " ",
        lapply(shown_keys(), function(k)
          span(class = "chip", span(class = "dot", style = sprintf("background:%s", COL[[k]])),
               sprintf("%s %s", nm[[k]], fmt_p(cv$y[[k]][i], 1)))))
  })

  output$pc_table <- renderTable({
    cv <- curves(); s <- par_in(); P <- plan()
    keep <- which(abs(cv$x * 20 - round(cv$x * 20)) < 1e-6 | abs(cv$x - cv$design_x) < 1e-9)
    lab <- series_labels(s, P)
    df <- data.frame(cv$xfmt(cv$x[keep]), check.names = FALSE)
    names(df) <- cv$xname
    for (k in shown_keys()) df[[strip_tags(lab[[k]])]] <- fmt_p(cv$y[[k]][keep], 1)
    df
  }, striped = TRUE, spacing = "s", align = "l")

  # --- sample-size table ------------------------------------------------------

  output$ss_table <- renderTable({
    P <- plan(); s <- par_in()
    if (s$endpoint == "PFS") {
      df <- data.frame(
        Method = c("Uncorrected (Schoenfeld; assumes perfect reading)",
                   "Calculator (false calls as an extra hazard)",
                   "Exact (recommended)",
                   "Size if every scan were read correctly (exact)"),
        Patients = fmt_n(P$n), Events = fmt_n(P$events),
        `Power with the stated reader` = fmt_p(P$power, 1), check.names = FALSE)
    } else {
      df <- data.frame(
        Method = c("Uncorrected (two-proportion formula on the true rates)",
                   "Exact (recommended)"),
        Patients = fmt_n(P$n[c("uncorrected", "corrected")]),
        `Power with the stated reader` = fmt_p(P$power[c("uncorrected", "corrected")], 1),
        check.names = FALSE)
    }
    if (!is.null(P$power_true)) {
      pt <- if (s$endpoint == "PFS") P$power_true else P$power_true[c("uncorrected", "corrected")]
      df[["Power if the true accuracy differs (as set)"]] <- fmt_p(pt, 1)
    }
    df
  }, striped = TRUE, spacing = "s", align = "l")

  output$ss_notes <- renderUI({
    P <- plan(); s <- par_in()
    d <- P$d
    common <- sprintf("%s test at α = %s, %s power, %s:1 allocation (new treatment : control).",
                      if (d$sd == 2) "Two-sided" else "One-sided", s$alpha, fmt_p(s$power),
                      formatC(s$ratio, format = "fg", digits = 3))
    extra <- switch(s$endpoint,
      PFS = sprintf("Events for the exact method are expected recorded events at %s months' follow-up. The calculator ignores missed progressions and uses the first scan interval; on a regular schedule with Se = 1 it is exact, and it needs at least %s× the uncorrected events.",
                    fmt_x(s$tau, 0), formatC((1 + P$f / (log(2) / s$m0))^2, format = "f", digits = 2)),
      DCR = sprintf("Disease control is recorded at the scan in week %s (scan %d) if the patient is alive and no scan up to then was read as PD. True DCR %s vs %s; recorded %s vs %s.",
                    formatC(sch()$weeks[P$kappa], format = "fg"), P$kappa, fmt_p(P$true[1], 1), fmt_p(P$true[2], 1), fmt_p(P$rec[1], 1), fmt_p(P$rec[2], 1)),
      ORR = sprintf("Response is recorded if any scan before the first PD read is read as PR%s. True ORR %s vs %s; recorded %s vs %s.",
                    if (s$confirmed) " on two consecutive scans" else "", fmt_p(P$true[1], 1), fmt_p(P$true[2], 1), fmt_p(P$rec[1], 1), fmt_p(P$rec[2], 1)))
    div(class = "note", p(common), p(extra))
  })

  # --- sample size against reader accuracy ------------------------------------

  sens <- reactive({
    s <- par_in(); P <- plan(); sc <- sch(); d <- P$d
    lev <- c(1, 0.9, 0.8, 0.7)
    if (s$endpoint == "ORR") {
      x <- seq(0, max(0.15, ceiling(s$fr * 100) / 100), by = 0.01)
      a1 <- orr_arm(s$orr1, P$lam1, s$pi_d / 100, sc)
      y <- lapply(lev, function(se_r) vapply(x, function(fr) {
        M <- reader_matrix(1, 1, se_r, fr)
        n_two_prop(orr_observed(P$a0, sc, M, s$confirmed), orr_observed(a1, sc, M, s$confirmed), d$a, d$pw, d$sd, d$w)
      }, 0))
      list(x = x, y = y, lev = lev, cur = c(s$fr, s$se_r), xlab = "Per-scan false-response probability, FR",
           leg = SER, xfmt = function(v) fmt_x(v, 2))
    } else {
      x <- seq(min(0.85, floor(s$sp * 100) / 100), 1, by = 0.01)
      a <- P$arms
      if (s$endpoint == "PFS") {
        G <- grid()
        f <- function(se, sp) {
          el <- expected_logrank(a[[1]], a[[2]], G, se, sp, d$w)
          ceiling(d$z2 * el$v / el$u^2)
        }
      } else {
        kp <- P$kappa
        f <- function(se, sp) n_two_prop(dcr_observed(a[[1]], sc, kp, se, sp), dcr_observed(a[[2]], sc, kp, se, sp), d$a, d$pw, d$sd, d$w)
      }
      y <- lapply(lev, function(se) vapply(x, function(sp) f(se, sp), 0))
      list(x = x, y = y, lev = lev, cur = c(s$sp, s$se), xlab = "Per-scan specificity, Sp",
           leg = "Se", xfmt = function(v) fmt_x(v, 2))
    }
  })

  output$sens <- renderPlot({
    S <- sens(); P <- plan()
    yl <- c(0, max(pretty(c(0, unlist(S$y), P$n[["corrected"]]))))
    par(mar = c(4.2, 4.8, 0.6, 1.4), mgp = c(3.1, 0.55, 0), tcl = -0.25, las = 1,
        col.axis = COL$text2, col.lab = COL$text1, cex.axis = 0.9, cex.lab = 0.95, bg = "white")
    plot(NA, xlim = range(S$x), ylim = yl, xaxs = "i", yaxs = "i", axes = FALSE,
         xlab = S$xlab, ylab = "Patients needed (exact method)")
    yt <- pretty(yl, 6); yt <- yt[yt <= yl[2]]
    abline(h = yt[yt > 0], col = COL$grid)
    axis(1, at = pretty(S$x, 6), labels = S$xfmt(pretty(S$x, 6)), col = COL$axis)
    axis(2, at = yt, labels = fmt_n(yt), lwd = 0)
    abline(h = P$n[["uncorrected"]], lty = 2, col = COL$uncorrected, lwd = 1.6)
    text(S$x[1], P$n[["uncorrected"]],
         sprintf("uncorrected, %s", fmt_n(P$n[["uncorrected"]])), adj = c(-0.05, -0.5), cex = 0.78, col = COL$text2)
    for (i in rev(seq_along(S$lev))) lines(S$x, S$y[[i]], col = SEQ[i], lwd = 2.4)
    points(S$cur[1], P$n[["corrected"]], pch = 21, bg = "white", col = COL$text1, lwd = 2, cex = 1.4)
  }, res = 96)

  output$sens_legend <- renderUI({
    S <- sens()
    items <- lapply(seq_along(S$lev), function(i)
      span(class = "legend-item",
           HTML(sprintf('<svg width="30" height="10" aria-hidden="true"><line x1="1" y1="5" x2="29" y2="5" stroke="%s" stroke-width="2.4"/></svg>', SEQ[i])),
           HTML(sprintf("%s = %s", S$leg, fmt_x(S$lev[i], 1)))))
    items[[length(items) + 1]] <- span(class = "legend-item",
      HTML('<svg width="16" height="16" aria-hidden="true"><circle cx="8" cy="8" r="5.5" fill="white" stroke="#0b0b0b" stroke-width="2"/></svg>'),
      "current setting")
    div(class = "legend-row", items)
  })

  output$sens_note <- renderUI({
    s <- par_in()
    msg <- if (s$endpoint == "ORR")
      "For ORR, false responses drive the correction: every scan a patient attends is another chance of one."
    else
      "False progression calls drive the correction; missed progressions (lower Se) matter much less, because they delay recorded events in both arms alike."
    div(class = "note", p(msg))
  })
}

shinyApp(ui, server)
