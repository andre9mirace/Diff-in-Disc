# Diff-in-Disc estimators compared in the Monte Carlo.
# Requires Bandwidth/SGD/BW_Functions.R (Imbens_BW, rdd_k_weights, penalized_MSE) and rdrobust.
#
# Every local linear fit uses the Grembi et al. (2016) specification. Its standard errors are
# clustered by property, since each property is observed in every year.

# Local linear Diff-in-Disc within |x - c| < h: returns the bandwidth, tau and its clustered SE
grembi_fit <- function(data, h, kernel = "triangular", c = 0, t0 = 3){
  s <- data[abs(data$Distance_to_Border - c) < h, ]
  x <- s$Distance_to_Border - c
  group <- as.numeric(x >= 0)
  post <- as.numeric(s$Year >= t0)
  did <- group * post
  X <- cbind(1, x, group, post, did, x*group, x*post, x*did)
  w <- rdd_k_weights(x = x, c = 0, h = h, kernel = kernel)

  # Weighted least squares
  XtW <- t(X * w)
  bread <- solve(XtW %*% X)
  beta <- bread %*% (XtW %*% s$Price_m2)
  e <- as.vector(s$Price_m2 - X %*% beta)

  # Cluster-robust variance (CR1) by property
  scores <- rowsum(X * (w * e), s$PropertyID)
  G <- nrow(scores)
  n <- nrow(X)
  V <- bread %*% crossprod(scores) %*% bread * (G / (G - 1)) * ((n - 1) / (n - ncol(X)))
  c(h = h, tau = beta[5], se = sqrt(V[5, 5]))
}

## Bandwidth selectors

# Imbens & Kalyanaraman (2012) plug-in on the pooled sample, as in the thesis
bw_ik <- function(data, kernel){
  Imbens_BW(y = data$Price_m2, x = data$Distance_to_Border, c = 0, kernel = kernel)
}

# Grembi et al. (2016): average of the rdrobust MSE-optimal bandwidths before and after t0
bw_avg <- function(data, kernel, t0 = 3){
  k <- if (kernel == "triangular") "tri" else "uni"
  pre <- data$Year < t0
  # Each property's x repeats across years, so rdrobust warns about mass points
  h_pre <- suppressWarnings(
    rdrobust(data$Price_m2[pre], data$Distance_to_Border[pre], kernel = k)$bws[1])
  h_post <- suppressWarnings(
    rdrobust(data$Price_m2[!pre], data$Distance_to_Border[!pre], kernel = k)$bws[1])
  (h_pre + h_post) / 2
}

# Minimizer of the adaptive SGD objective J(h) = penalty(h) * MSE(h). DiRD_BW minimizes the
# same J; this is a faster minimizer for many replications (run_montecarlo.R compares both).
# J can have narrow minima (e.g. just below a break in the data) and several local minima, so
# it is evaluated on a dense grid of quantiles of |x - c|, and the best few local minima of
# the grid are refined by golden-section search.
bw_sgd_criterion <- function(data, kernel, t0 = 3, n_grid = 400, n_refine = 3){
  d <- data.frame(ID = data$PropertyID, y = data$Price_m2, x = data$Distance_to_Border,
                  time_var = data$Year)
  N <- nrow(d)
  J <- function(h) penalized_MSE(d, h, c = 0, t0 = t0, N = N, kernel = kernel)
  grid <- unique(quantile(abs(d$x), seq(0.005, 0.995, length.out = n_grid), names = FALSE))
  J_grid <- sapply(grid, J)
  J_grid[is.na(J_grid)] <- Inf

  # Local minima of the grid, best first
  n <- length(grid)
  is_min <- J_grid <= c(Inf, J_grid[-n]) & J_grid <= c(J_grid[-1], Inf)
  candidates <- head(order(ifelse(is_min, J_grid, Inf)), n_refine)
  candidates <- candidates[is.finite(J_grid[candidates])]

  best_h <- grid[candidates[1]]
  best_J <- J_grid[candidates[1]]
  for (i in candidates) {
    opt <- optimize(J, c(grid[max(i - 1, 1)], grid[min(i + 1, n)]))
    if (!is.na(opt$objective) && opt$objective < best_J) {
      best_h <- opt$minimum
      best_J <- opt$objective
    }
  }
  best_h
}

## Estimators

# rdrobust on the property-level change in price (mean after t0 minus mean before t0).
# With a panel and a time-invariant running variable this is the Diff-in-Disc estimand, and
# rdrobust's MSE-optimal bandwidth and robust bias-corrected CI target tau directly.
fd_rdrobust <- function(data, kernel, t0 = 3){
  k <- if (kernel == "triangular") "tri" else "uni"
  before <- data$Year < t0
  y_pre <- tapply(data$Price_m2[before], data$PropertyID[before], mean)
  y_post <- tapply(data$Price_m2[!before], data$PropertyID[!before], mean)
  x <- tapply(data$Distance_to_Border, data$PropertyID, `[`, 1)
  r <- rdrobust(y_post - y_pre, x, kernel = k)
  c(h = r$bws[1, 1], tau = r$coef[1], se = r$se[1], ci_low = r$ci[3, 1], ci_high = r$ci[3, 2])
}

# Every method and kernel on one dataset, plus Grembi fits at fixed bandwidths (used to find
# the infeasible RMSE-optimal bandwidth across replications)
estimate_all <- function(data, fixed_h, t0 = 3){
  rows <- list()
  add <- function(method, kernel, est){
    ci_low <- if ("ci_low" %in% names(est)) est[["ci_low"]] else est[["tau"]] - 1.96 * est[["se"]]
    ci_high <- if ("ci_high" %in% names(est)) est[["ci_high"]] else est[["tau"]] + 1.96 * est[["se"]]
    rows[[length(rows) + 1]] <<- data.frame(method = method, kernel = kernel, h = est[["h"]],
                                            tau = est[["tau"]], se = est[["se"]],
                                            ci_low = ci_low, ci_high = ci_high)
  }
  for (k in c("triangular", "uniform")) {
    add("IK plug-in", k, grembi_fit(data, bw_ik(data, k), k, t0 = t0))
    add("AVG rdrobust", k, grembi_fit(data, bw_avg(data, k, t0), k, t0 = t0))
    add("SGD criterion", k, grembi_fit(data, bw_sgd_criterion(data, k, t0), k, t0 = t0))
    add("FD rdrobust", k, fd_rdrobust(data, k, t0))
    for (h in fixed_h) add("Fixed h", k, grembi_fit(data, h, k, t0 = t0))
  }
  do.call(rbind, rows)
}
