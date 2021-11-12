# Define global variables to avoid
# notes from R CMD check
# P <- NULL
# N <- NULL
# ptail <- NULL
# sim_id <- NULL
# a1 <- b1 <- a2 <- b2 <- NULL
# parm <- stage <- variable <- value <- grp <- NULL


#' Calculate density of beta-binomial distribution
#'
#' @param x The value at which to evaluate the density
#' @param n The sample size
#' @param a First parameter
#' @param b Second parameter
#' @return Value of beta-binomial(n,a,b) evaluated at x
#' @examples
#' dbetabinom(5, 10, 2, 3)
#' @export
dbetabinom <- function(x, n, a = 1, b = 1) {
  stopifnot("x must be >= 0" = all(x >= 0))
  stopifnot("n must be > 0" = all(n >= 1))
  stopifnot("a and b must be > 0" = all(c(a, b) > 0))
  num <- lgamma(a + b) + lgamma(n + 1) + lgamma(x + a) + lgamma(n - x + b)
  den <- lgamma(a) + lgamma(b) + lgamma(x + 1) + lgamma(n - x + 1) + lgamma(n + a + b)
  prob <- exp(num - den)
  return(prob)
}


#' Compare Normal approximation to Beta
#'
#' @param a First Beta parameter
#' @param b Second Beta parameter
#' @param ... Other arguments to `curve()`
#' @return A plot of Beta density and Normal approximation
#' @export
plot_beta_norm <- function(a, b, ...) {
  graphics::curve(stats::dbeta(x, a, b), ...)
  graphics::curve(stats::dnorm(x, a / (a + b), sqrt(a * b / ((a + b)^2 * (a + b + 1)))), add = TRUE, col = "red", ...)
}


#' Draw random variates from beta-binomial distribution
#'
#' @param n The number of random values to sample
#' @param m The sample size
#' @param a First parameter
#' @param b Second parameter
#' @examples
#' rbetabinom(2, 10, 2, 3)
#' @export
rbetabinom <- function(n, m, a = 1, b = 1) {
  stopifnot("n must be > 0" = all(n > 0))
  stopifnot("a and b must be > 9" = all(c(a, b) > 0))
  return(stats::rbinom(n, m, stats::rbeta(n, a, b)))
}


#' Calculate Pr(X > Y + delta) where X and Y are independent Beta random variables
#' using numerical integration
#'
#' @param a Parameter one of beta density for X
#' @param b Parameter two of beta density for X
#' @param c Parameter one of beta density for Y
#' @param d Parameter two of beta density for Y
#' @param delta The difference we wish to assess (i.e. X - Y > delta)
#' @param ... other arguments passed to integrate/quadgk function
#' @return The value of the integral
#' @examples
#' beta_ineq(5, 5, 3, 7)
#' @export
#' @importFrom pracma quadgk
beta_ineq <- function(a, b, c, d, delta = 0, ...) {
  stopifnot("a, b, c, d must be > 0" = all(c(a, b, c, d) > 0))
  integrand <- function(x) {
    stats::dbeta(x, a, b) * stats::pbeta(x - delta, c, d)
  }
  tryCatch(
    pracma::quadgk(integrand, delta, 1, ...),
    error = function(err) NA
  )
}


#' Calculate Pr(X > Y + delta) where X and Y are independent Beta random variables
#' using Normal approximation.
#'
#' @param a Parameter one of beta density for X
#' @param b Parameter two of beta density for X
#' @param c Parameter one of beta density for Y
#' @param d Parameter two of beta density for Y
#' @param delta The difference we wish to assess (i.e. X - Y > delta)
#' @return The value of the integral
#' @examples
#' beta_ineq_approx(5, 5, 3, 7)
#' @export
beta_ineq_approx <- function(a, b, c, d, delta = 0) {
  stopifnot("a, b, c, d must be > 0" = all(c(a, b, c, d) > 0))
  m1 <- a / (a + b)
  v1 <- a * b / ((a + b)^2 * (a + b + 1))
  m2 <- c / (c + d)
  v2 <- c * d / ((c + d)^2 * (c + d + 1))
  z <- (m1 - m2 - delta) / sqrt(v1 + v2)
  return(stats::pnorm(z))
}


#' Calculate Pr(X > Y + delta) where X and Y are independent Beta random variables
#' using Monte Carlo method.
#'
#' @param a Parameter one of beta density for X
#' @param b Parameter two of beta density for X
#' @param c Parameter one of beta density for Y
#' @param d Parameter two of beta density for Y
#' @param delta The difference we wish to assess (i.e. X - Y > delta)
#' @param sims The number of Monte Carlo variates to generate for estimation
#' @return The value of the integral
#' @examples
#' beta_ineq_sim(5, 5, 3, 7)
#' @export
beta_ineq_sim <- function(a, b, c, d, delta = 0, sims = 10000) {
  stopifnot("a, b, c, d must be > 0" = all(c(a, b, c, d) > 0))
  lens <- unlist(lapply(list(a, b, c, d), length))
  stopifnot("a, b, c, d must be same length" = all(max(lens) - min(lens) == 0))
  X <- lapply(1:length(a), function(x) stats::rbeta(sims, a[x], b[x]))
  Y <- lapply(1:length(a), function(x) stats::rbeta(sims, c[x], d[x]))
  means <- lapply(1:length(a), function(x) mean(X[[x]] > Y[[x]] + delta))
  return(unlist(means))
}


#' Calculate the predicted probability of success
#'
#' @import data.table
#' @param a First parameter of first beta random variable
#' @param b Second parameter of first beta random variable
#' @param c First paramter of second beta random variable
#' @param d Second parameter of second beta random variable
#' @param m1 Sample size to predict for first beta random variable
#' @param m2 Sample size to predict for second beta random variable
#' @param k_ppos The posterior probability cut-point to be assessed
#' @param post_method The method to use for calculating posterior probabilities,
#' one of "exact" (numerical), "approx", "sim".
#' @param post_sim Number of posterior simulations if post_method = "sim".
#' @return The predicted probability of success
#' @export
calc_ppos <- function(a, b, c, d, m1, m2, k_ppos,
                      post_method = "exact", post_sim = 1e4) {
  stopifnot("a, b, c, d, m1, and m2 must be > 0" = all(c(a, b, c, d, m1, m2) > 0))
  stopifnot("k_ppos must be in [0, 1]" = (k_ppos >= 0 & k_ppos <= 1))
  calc_post <- switch(
    post_method,
    "exact" = beta_ineq,
    "approx" = beta_ineq_approx,
    "sim" = beta_ineq_sim)
  y1pred <- rbetabinom(post_sim, m1, a, b)
  y2pred <- rbetabinom(post_sim, m2, c, d)
  ypred <- data.table(y1pred = y1pred, y2pred = y2pred)[, .N, keyby = list(y1pred, y2pred)]
  ypred[, `:=`(P = Vectorize(calc_post)(a + y1pred,
                                        b + m1 - y1pred,
                                        c + y2pred,
                                        d + m2 - y2pred))]
  return (ypred[, c(sum(N * (P > k_ppos)) / sum(N))])
}


#' Simulate trial data using Poisson process for accr
#'
#' @param p1tru True response rate under control
#' @param p2tru True response rate under treatment
#' @param nmax The maximum total sample size
#' @param enro_rate The baseline accrual rate
#' @param enro_intensity Function for changing accrual intensity
#' @param resp_delay Function returning response time value
#' @param simple_rand Use complete randomisation or assume equal allocation at each interim
#' @return A data.table
#' @export
#' @importFrom randomizr complete_ra
sim_trial_dat <- function(
  p1tru  = 0.1,
  p2tru  = 0.1,
  nmax = 100,
  enro_rate = 1,
  enro_intensity = function(t) 1,
  resp_delay = function(n) stats::runif(n),
  simple_rand = TRUE
) {
  if(!(all(c(p1tru, p2tru) > 0) & all(c(p1tru, p2tru) < 1)))
    stop("p1tru and p2tru must be in (0, 1)")

  enro_t <- poisson::nhpp.sim(rate = enro_rate, num.events = nmax, enro_intensity, prepend.t0 = F)
  resp_t <- enro_t + resp_delay(nmax)
  if (simple_rand) {
    x <- as.numeric(randomizr::complete_ra(nmax, num_arms = 2))
  } else {
    x <- rep(c(1, 2), times = nmax/2)
  }
  d <- data.table(p1tru = p1tru,
                  p2tru = p2tru,
                  pid = 1:nmax,
                  enro_t = enro_t,
                  resp_t = resp_t,
                  x = x)
  d[order(resp_t), resp_o := 1:.N]
  d[, y := stats::rbinom(nmax, 1, prob = ifelse(x == 1, p1tru, p2tru))]
  return(d)
}


#' Aggregate trial data over treatment
#'
#' @param d The raw trial data.
#' @param stage_n The number of responses required to trigger each interim.
#' @param min_rem The minimum number of subjects still to be enrolled for an interim to trigger.
#' @return An aggregated trial dataset.
#' @export
agg_trial_dat <- function(d, stage_n, min_rem = 10) {
  resp_cut <- sort(d[, resp_t])[stage_n]

  # for all cutpoints, by group
  # n participants with resp_time <= cut point
  # y participants with events with resp_time <= cut point
  # m participants enrolled but without resp_time available at cut point
  # w participants with events that are enrolled but without resp_time available at cut point
  # l participants with enrolment occurring after cut point
  # z participants with events with enrolment occurring after cut point
  dd <- d[,
          {
            list(resp = stage_n,
                 n = sapply(resp_cut, function(a) sum(resp_t <= a)),
                 y = sapply(resp_cut, function(a) sum((resp_t <= a)*y)),
                 m = sapply(resp_cut, function(a) sum(resp_t > a & enro_t <= a)),
                 w = sapply(resp_cut, function(a) sum((resp_t > a & enro_t <= a)*y)),
                 l = sapply(resp_cut, function(a) sum(enro_t > a)),
                 z = sapply(resp_cut, function(a) sum((enro_t > a)*y))
            )
          }, by = .(p1tru, p2tru, x)]

  # spread to wide.
  dcast(dd, resp + p1tru + p2tru ~ x,
        value.var = c("n", "y", "m", "w", "l", "z"), sep = "")[(l1 + l2) > min_rem | resp == max(stage_n)]
}


#' Calculate trial probabilities (posterior and predictive)
#'
#' @param d The Trial data
#' @param ppos_q The vector of cut-points to consider
#' @param ppos_sim Number of simulations from posterior predictive to use in PPoS(q) calculation
#' @param post_method Method used to calculate P(X > Y + delta), one of `exact` (integration), `approx` (normal), `sim` (monte carlo)
#' @param a1 Prior shape for treatment 1
#' @param b1 Prior scale for treatment 1
#' @param a2 Prior shape for treatment 2
#' @param b2 Prior scale for treatment 2
#' @return Updates `d` inplace but also returns the updated `d`.
#' `post` - current posterior probability,
#' `post_int` - posterior probability when follow-up enrolled individuals,
#' `post_fin` - posterior probability when follow-up to maximum sample size,
#' `ppos_intq` - predictive probability of success given q if follow-up enrolled individuals,
#' `ppos_finq` - predictive probability of success given q if follow-up to maximum sample size,
#'
#' @export
est_trial_prob <- function(
  d,
  ppos_q = 0.95,
  ppos_sim = 1000,
  post_method = "approx",
  a1 = 1, b1 = 1, a2 = 1, b2 = 1
) {
  if(!(all(c(a1, b1, a2, b2) > 0)))
    stop("a1, b1, a2, b2 must be > 0")
  if(!post_method %in% c("exact", "approx", "sim"))
    stop("post_method must be 'exact', 'approx', or 'sim'.")

  # What method used to estimate posterior probability?
  calc_post <- switch(post_method,
                      "exact" = beta_ineq,
                      "approx" = beta_ineq_approx,
                      "sim" = beta_ineq_sim)
  # Setup storage for prediction
  ypred <- data.table(y1 = rep(0, ppos_sim),
                      y2 = rep(0, ppos_sim))
  # Determine posterior parameters
  d[, `:=`(a1 = a1 + y1,
           b1 = b1 + n1 - y1,
           a2 = a2 + y2,
           b2 = a2 + n2 - y2)]
  # Calculate posterior probabilities given current data and assuming follow-up of currently enrolled
  d[, `:=`(post = Vectorize(calc_post)(a1, b1, a2, b2),
           post_int = Vectorize(calc_post)(a1 + w1, b1 + m1 - w1, a2 + w2, b2 + m2 - w2))]
  d[, post_fin := post[.N]]

  # Calculate PPoS at each interim for currently enrolled and complete enrolment
  if(!is.null(ppos_q)) {
    for(i in 1:(nrow(d) - 1)) {
      # Do interim PPoS calculation
      ypred[, `:=`(y1 = rbetabinom(ppos_sim, d[i, m1], d[i, a1], d[i, b1]),
                   y2 = rbetabinom(ppos_sim, d[i, m2], d[i, a2], d[i, b2]))]
      ypred_agg <- ypred[, .N, by = .(y1, y2)]
      ypred_agg[, P := Vectorize(calc_post)(d[i, a1] + y1,
                                            d[i, b1 + m1] - y1,
                                            d[i, a2] + y2,
                                            d[i, b2 + m2] - y2)]
      d[i, paste0('ppos_int', ppos_q) := ypred_agg[, lapply(ppos_q, function(q) sum(N * (P > q)) / sum(N))]]

      # Do final PPoS calculation
      ypred[, `:=`(y1 = rbetabinom(ppos_sim, d[i, m1 + l1], d[i, a1], d[i, b1]),
                   y2 = rbetabinom(ppos_sim, d[i, m2 + l2], d[i, a2], d[i, b2]))]
      ypred_agg <- ypred[, .N, by = .(y1, y2)]
      ypred_agg[, P := Vectorize(calc_post)(d[i, a1] + y1,
                                            d[i, b1 + m1 + l1] - y1,
                                            d[i, a2] + y2,
                                            d[i, b2 + m2 + l2] - y2)]
      d[i, paste0('ppos_fin', ppos_q) := ypred_agg[, lapply(ppos_q, function(q) sum(N * (P > q)) / sum(N))]]
    }
  }
  return(d)
}


#' Apply decision rule to trial data
#'
#' @param trial Aggregated trial data with estimated probabilities (after call to est_trial_prob).
#' @param fut_k Futility boundary vector. Must has as many elements as the maximum number of interims in `trial`.
#' @param suc_k Success boundary vector.
#' @param inf_k Inferiority cut-off.
#' @param sup_k Superiority boundary vector.
#' @return A new `data.table` given the decision made and boundaries used
#' @export
dec_trial <- function(
  trial,
  fut_k = 0.1,
  suc_k = 0.9,
  inf_k = 0.05,
  sup_k = 0.95) {

  # Use the ppos column which matches sup_k
  ppos_cols <- grep(paste0(sup_k, "$"), names(trial), value = T)
  if(length(ppos_cols) == 0) stop("sup_k did not match any columns in trial.")
  max_stage <- trial[, .N, by = sim_id][, max(N)]

  if (length(fut_k) == 1)
    fut_k <- rep(fut_k, max_stage - 1)
  if (length(suc_k) == 1)
    suc_k <- rep(suc_k, max_stage - 1)
  if (length(fut_k) < (max_stage - 1))
    stop("fut_k has too few elements")
  if (length(suc_k) < (max_stage - 1))
    stop("suc_k has too few elements")

  trial[,
        {
          fut <- match(TRUE, get(ppos_cols[2])[-.N] < fut_k[1:(.N-1)])
          suc <- match(TRUE, get(ppos_cols[1])[-.N] > suc_k[1:(.N-1)])
          if(any(!is.na(c(fut, suc)))) {
            res <- which.min(c(fut, suc))
            int <- min(c(fut, suc), na.rm = TRUE)
            list(p1tru = p1tru[1],
                 p2tru = p2tru[1],
                 res = switch(res, "futile", "expect success"),
                 fin = ifelse(post_int[int] > sup_k, "success", "failure"),
                 stage = int,
                 resp = n1[int] + n2[int],
                 enro = n1[int] + n2[int] + m1[int] + m2[int],
                 post = post[int],
                 post_int = post_int[int],
                 post_fin = post_fin[int],
                 ppos_int = get(ppos_cols[1])[int],
                 ppos_fin = get(ppos_cols[2])[int],
                 fut_k = paste(unique(fut_k), collapse = ","),
                 suc_k = paste(unique(suc_k), collapse = ","),
                 inf_k = inf_k,
                 sup_k = sup_k)
          } else {
            inf <- utils::tail(post, 1) < utils::tail(inf_k, 1)
            sup <- utils::tail(post, 1) > utils::tail(sup_k, 1)
            list(p1tru = p1tru[1],
                 p2tru = p2tru[1],
                 res = ifelse(inf, "inferior", ifelse(sup, "superior", "inconclusive")),
                 fin = ifelse(post[.N] > sup_k, "success", "failure"),
                 stage = .N,
                 resp = n1[.N] + n2[.N],
                 enro = n1[.N] + n2[.N] + m1[.N] + m2[.N],
                 post = post[.N],
                 post_int = post_int[.N],
                 post_fin = post_fin[.N],
                 ppos_int = get(ppos_cols[1])[.N],
                 ppos_fin = get(ppos_cols[2])[.N],
                 fut_k = paste(unique(fut_k), collapse = ","),
                 suc_k = paste(unique(suc_k), collapse = ","),
                 inf_k = inf_k,
                 sup_k = sup_k)
          }
        }, by = sim_id]
}
