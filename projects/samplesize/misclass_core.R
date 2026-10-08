# misclass_core.R -- planning functions for the sample-size app.
#
# Sample size for trials whose endpoints are read from scans: correcting for
# false and missed progression calls in objective response, disease control
# and progression-free survival.
#
# Same model and formulas as misclass_ss.R in the paper's code supplement.
# The recorded PFS survival function is evaluated in vectorised form so the
# power curves draw quickly when R runs in the browser (webR); it agrees with
# misclass_ss.R to machine precision on the same time grid.
#
# Model: constant hazards. In arm a, progression hazard lamP and death without
# progression lamD. An undetected progressor dies at hazard mu. Each scan is
# read independently given the true state, identically in both arms, with
# per-scan sensitivity Se and specificity Sp for progression. A PD read ends
# assessment; deaths are observed exactly. Times are in months.

WEEK <- 7 / 30.4375

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

make_arm <- function(lamP, lamD, rho = 0) list(lamP = lamP, lamD = lamD, rho = rho)

lam_tot <- function(arm) arm$lamP + arm$lamD

arms_from_median <- function(m0, hr, pi_death, rho0 = 0, rho1 = 0) {
  lam0 <- log(2) / m0
  list(make_arm((1 - pi_death) * lam0, pi_death * lam0, rho0),
       make_arm((1 - pi_death) * lam0 * hr, pi_death * lam0 * hr, rho1))
}

# Scans every int1 weeks up to week switch_wk, then every int2 weeks, until
# the end of follow-up tau (months).
make_schedule <- function(int1, switch_wk, int2, tau, mu) {
  until <- tau / WEEK
  v <- numeric(0); t <- 0
  while (t + int1 <= switch_wk + 1e-9 && t + int1 <= until + 1e-9) {
    t <- t + int1; v <- c(v, t)
  }
  if (!is.na(int2) && int2 > 0)
    while (t + int2 <= until + 1e-9) { t <- t + int2; v <- c(v, t) }
  list(visits = v * WEEK, weeks = v, tau = tau, mu = mu)
}

# Confusion matrix: rows latent R, S, P; columns read PR, SD, PD.
reader_matrix <- function(se = 1, sp = 1, se_r = 1, fr = 0) {
  fp <- 1 - sp
  M <- rbind(c(se_r, 1 - se_r - fp, fp),
             c(fr, 1 - fr - fp, fp),
             c(0, 1 - se, se))
  if (any(M < -1e-12)) stop("infeasible reader")
  pmin(pmax(M, 0), 1)
}

# P(progress in (a, b], alive and undetected at t >= b); vectorised in t
prog_alive_undetected <- function(arm, mu, a, b, t) {
  d <- lam_tot(arm) - mu
  if (abs(d) < 1e-12) return(arm$lamP * exp(-mu * t) * (b - a))
  arm$lamP * exp(-mu * t) * (exp(-d * a) - exp(-d * b)) / d
}

# ---------------------------------------------------------------------------
# Binary endpoints
# ---------------------------------------------------------------------------

n_two_prop <- function(p0, p1, alpha = 0.05, power = 0.80, sided = 2, w = 0.5) {
  za <- qnorm(1 - alpha / sided); zb <- qnorm(power)
  pbar <- (1 - w) * p0 + w * p1
  num <- (za * sqrt(pbar * (1 - pbar) * (1 / (1 - w) + 1 / w)) +
          zb * sqrt(p0 * (1 - p0) / (1 - w) + p1 * (1 - p1) / w))^2
  ceiling(num / (p1 - p0)^2)
}

power_two_prop <- function(n, p0, p1, alpha = 0.05, sided = 2, w = 0.5) {
  za <- qnorm(1 - alpha / sided)
  pbar <- (1 - w) * p0 + w * p1
  s0 <- sqrt(pbar * (1 - pbar) * (1 / (1 - w) + 1 / w) / n)
  s1 <- sqrt((p0 * (1 - p0) / (1 - w) + p1 * (1 - p1) / w) / n)
  d <- p1 - p0
  out <- pnorm((d - za * s0) / s1)
  if (sided == 2) out <- out + pnorm((-d - za * s0) / s1)
  out
}

dcr_true <- function(arm, sch, kappa) exp(-lam_tot(arm) * sch$visits[kappa])

# Exact recorded DCR at scan kappa
dcr_observed <- function(arm, sch, kappa, se, sp) {
  v <- c(0, sch$visits); t <- v[kappa + 1]
  out <- sp^kappa * dcr_true(arm, sch, kappa)
  for (j in seq_len(kappa))
    out <- out + sp^(j - 1) * (1 - se)^(kappa - j + 1) *
      prog_alive_undetected(arm, sch$mu, v[j], v[j + 1], t)
  out
}

orr_true <- function(arm, sch) arm$rho * exp(-lam_tot(arm) * sch$visits[1])

# Exact recorded ORR (unconfirmed or confirmed) by recursion over scans
orr_observed <- function(arm, sch, M, confirmed = FALSE) {
  v <- c(0, sch$visits); lam <- lam_tot(arm); mu <- sch$mu
  total <- 0
  for (cls in 1:2) {                     # 1 = latent R, 2 = latent S
    wgt <- if (cls == 1) arm$rho else 1 - arm$rho
    x0 <- exp(-lam * v[2]); x1 <- 0
    y0 <- prog_alive_undetected(arm, mu, 0, v[2], v[2]); y1 <- 0
    pr <- 0
    K <- length(v) - 1
    for (k in seq_len(K)) {
      pr <- pr + if (confirmed) x1 * M[cls, 1] + y1 * M[3, 1]
                 else (x0 + x1) * M[cls, 1] + (y0 + y1) * M[3, 1]
      if (k == K) break
      d <- v[k + 2] - v[k + 1]
      eL <- exp(-lam * d); eM <- exp(-mu * d)
      q <- prog_alive_undetected(arm, mu, 0, d, d)
      if (confirmed) {
        nx0 <- eL * (x0 + x1) * M[cls, 2]; nx1 <- eL * x0 * M[cls, 1]
        ny0 <- eM * (y0 + y1) * M[3, 2] + q * (x0 + x1) * M[cls, 2]
        ny1 <- eM * y0 * M[3, 1] + q * x0 * M[cls, 1]
      } else {
        nx0 <- eL * x0 * M[cls, 2]; nx1 <- 0
        ny0 <- eM * y0 * M[3, 2] + q * x0 * M[cls, 2]; ny1 <- 0
      }
      x0 <- nx0; x1 <- nx1; y0 <- ny0; y1 <- ny1
    }
    total <- total + wgt * pr
  }
  total
}

# ---------------------------------------------------------------------------
# PFS
# ---------------------------------------------------------------------------

# Time grid used by the exact method: steps of h between scans, and every scan
# time twice (before and after its read). N = number of reads by each time.
pfs_grid <- function(sch, h = 0.02) {
  v <- c(0, sch$visits); K <- length(sch$visits)
  stepseq <- function(from, to) if (from > to) numeric(0) else seq(from, to, by = h)
  tt <- vector("list", 3 * K + 1); nn <- tt; i <- 0
  for (k in seq_len(K) + 1) {
    s <- stepseq(v[k - 1] + h, v[k] - 1e-9); s <- s[s < v[k] - 1e-9]
    tt[[i + 1]] <- s;    nn[[i + 1]] <- rep(k - 2, length(s))
    tt[[i + 2]] <- v[k]; nn[[i + 2]] <- k - 2
    tt[[i + 3]] <- v[k]; nn[[i + 3]] <- k - 1
    i <- i + 3
  }
  s <- stepseq(v[K + 1] + h, sch$tau + 1e-9); s <- s[s <= sch$tau + 1e-9]
  tt[[i + 1]] <- s; nn[[i + 1]] <- rep(K, length(s))
  list(t = unlist(tt), N = unlist(nn), v = v, tau = sch$tau, mu = sch$mu)
}

# Recorded PFS survival on the grid (vectorised). With g(a, b) the integral of
# lamP e^{-lam x} e^{mu x} over (a, b], and C_N the contribution of
# progressions in completed intervals that every later read has missed,
#   S(t) = e^{-lam t} Sp^N + lamP e^{-mu t} [C_N + Sp^N g(v_N, t)].
pfs_survival_grid <- function(arm, G, se, sp) {
  lam <- lam_tot(arm); d <- lam - G$mu; v <- G$v; K <- length(v) - 1
  g <- if (abs(d) < 1e-12) function(a, b) b - a else
    function(a, b) (exp(-d * a) - exp(-d * b)) / d
  C <- numeric(K + 1)
  for (n in seq_len(K))
    C[n + 1] <- (1 - se) * C[n] + sp^(n - 1) * (1 - se) * g(v[n], v[n + 1])
  N <- G$N; t <- G$t
  c(1, exp(-lam * t) * sp^N + arm$lamP * exp(-G$mu * t) * (C[N + 1] + sp^N * g(v[N + 1], t)))
}

# Per-patient drift u and variance v of the log-rank statistic, and the
# probability of a recorded event by tau; w = fraction on the experimental arm.
expected_logrank <- function(arm0, arm1, G, se, sp, w = 0.5, S0 = NULL) {
  if (is.null(S0)) S0 <- pfs_survival_grid(arm0, G, se, sp)
  S1 <- pfs_survival_grid(arm1, G, se, sp)
  m <- length(S0)
  y0 <- (1 - w) * S0[-m]; y1 <- w * S1[-m]
  dL0 <- 1 - S0[-1] / S0[-m]; dL1 <- 1 - S1[-1] / S1[-m]
  Y <- y0 + y1; r <- y1 / Y
  dbar <- (y0 * dL0 + y1 * dL1) / Y
  list(u = sum(y0 * y1 / Y * (dL1 - dL0)),
       v = sum(r * (1 - r) * (y0 * dL0 + y1 * dL1) * (1 - dbar)),
       ev = 1 - ((1 - w) * S0[m] + w * S1[m]))
}

n_pfs_exact <- function(arm0, arm1, G, se, sp, alpha = 0.05, power = 0.80,
                        sided = 2, w = 0.5) {
  el <- expected_logrank(arm0, arm1, G, se, sp, w)
  n <- (qnorm(1 - alpha / sided) + qnorm(power))^2 * el$v / el$u^2
  c(n = ceiling(n), events = n * el$ev)
}

# Power of the log-rank test with n patients; el from expected_logrank.
# Benefit (HR < 1) gives u < 0. Two-sided power counts both directions.
power_logrank <- function(n, el, alpha = 0.05, sided = 2) {
  za <- qnorm(1 - alpha / sided)
  d <- sqrt(n) * (-el$u) / sqrt(el$v)
  out <- pnorm(d - za)
  if (sided == 2) out <- out + pnorm(-d - za)
  out
}

schoenfeld_events <- function(hr, alpha = 0.05, power = 0.80, sided = 2, w = 0.5)
  (qnorm(1 - alpha / sided) + qnorm(power))^2 / (w * (1 - w) * log(hr)^2)

fp_hazard <- function(sp, Delta) -log(sp) / Delta

hr_recorded <- function(hr, lam0, f) (hr + f / lam0) / (1 + f / lam0)

n_pfs_calculator <- function(arm0, arm1, tau, sp, Delta, alpha = 0.05, power = 0.80,
                             sided = 2, w = 0.5) {
  f <- fp_hazard(sp, Delta)
  r0 <- lam_tot(arm0) + f; r1 <- lam_tot(arm1) + f
  D <- schoenfeld_events(r1 / r0, alpha, power, sided, w)
  ev <- 1 - ((1 - w) * exp(-r0 * tau) + w * exp(-r1 * tau))
  c(n = ceiling(D / ev), events = D, hr_obs = r1 / r0)
}

n_pfs_naive <- function(arm0, arm1, tau, alpha = 0.05, power = 0.80, sided = 2, w = 0.5) {
  hr <- lam_tot(arm1) / lam_tot(arm0)
  D <- schoenfeld_events(hr, alpha, power, sided, w)
  ev <- 1 - ((1 - w) * exp(-lam_tot(arm0) * tau) + w * exp(-lam_tot(arm1) * tau))
  c(n = ceiling(D / ev), events = D)
}
