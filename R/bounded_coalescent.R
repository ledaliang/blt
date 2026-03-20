#' Bounded-coalescent MLE of effective population size
#'
#' Computes a piecewise-constant maximum likelihood estimate of effective
#' population size under the bounded coalescent.
#'
#' @param data A \code{phylo} object or a list/data.frame with
#'   \code{coal_times}, \code{samp_times}, and \code{n_sampled}.
#' @param bound Numeric upper bound on the TMRCA.
#' @param lengthout Number of grid points.
#' @param eps Convergence tolerance for gradient ascent.
#' @param eta Step size for gradient ascent.
#'
#' @return A list with components \code{grid}, \code{x}, and \code{effpop}.
#' @export
mle_bounded_coal <- function(data, bound, lengthout = 5, eps = .02, eta = .01) {
  if (inherits(data, "phylo")) {
    phy <- phylodyn::summarize_phylo(data)
    n <- data$Nnode + 1
  } else if (all(c("coal_times", "samp_times", "n_sampled") %in%
                 names(data))) {
    phy <- with(data, list(
      samp_times = samp_times, coal_times = coal_times,
      n_sampled = n_sampled
    ))
  }
  grid_bds <- range(0, bound + .0001)
  grid <- seq(grid_bds[1], grid_bds[2], length.out = lengthout)
  lik_init <- coal_lik_init(
    samp_times = phy$samp_times, n_sampled = phy$n_sampled,
    coal_times = phy$coal_times, grid = grid
  )
  f_init <- rep(0, lik_init$ng)
  eta <- 0.01
  eps <- 0.02
  gradResult <- n_e_gradient_ascent(f_init, lik_init, bound, eps, eta)
  return(list(grid = grid, x = grid, effpop = exp(c(gradResult[[1]][1], gradResult[[1]]))))
}


#' Bounded coalescent log-likelihood
#'
#' Computes the bounded coalescent log-likelihood and its gradient for a
#' piecewise-constant effective population size trajectory.
#'
#' @param init Initialization object returned by
#'   \code{coal_lik_init()}.
#' @param f Vector of log effective population size values.
#'
#' @return A list with components \code{ll} (log-likelihood) and
#'   \code{dll} (gradient).
coal_loglik_bounded <- function(init, f) {
  if (init$ng != length(f)) {
    stop(paste("Incorrect length for f; should be", init$ng))
  }
  f <- rep(f, init$gridrep)

  ntip <- sum(init$ns)

  r_func <- function(k, j) {
    if (j == 1) {
      return(1)
    }
    prod <- 1
    for (m in 1:(j - 1)) {
      prod <- prod * ((2 * m + 1) / (2 * m - 1)) * ((k - m) / (k + m))
    }
    (-1)^(j - 1) * prod
  }

  r_ntip <- sapply(seq_len(ntip), function(i) r_func(ntip, i))
  com_vec <- choose(seq_len(ntip), 2)

  llnocoal <- init$D * init$C * exp(-f)
  sllnocoal <- init$D * exp(-f)

  lambda <- sum(sllnocoal)
  bound_prob <- sum(r_ntip * exp(-com_vec * lambda))
  ll_vec <- -init$y * f - llnocoal
  ll <- sum(ll_vec[!is.nan(ll_vec)]) - log(bound_prob)

  grad_bound <- sum(r_ntip * com_vec * exp(-com_vec * lambda))

  dll <- apply(init$rep_idx, 1, function(idx) {
    sum(-init$y[idx[1]:idx[2]] + llnocoal[idx[1]:idx[2]])
  }) - (grad_bound / bound_prob) * apply(init$rep_idx, 1, function(idx) {
    sum(sllnocoal[idx[1]:idx[2]])
  })
  list(ll = ll, dll = dll)
}


#' Gradient ascent for bounded coalescent MLE
#'
#' Runs gradient ascent on the bounded coalescent log-likelihood.
#'
#' @param f_init Initial vector of log effective population size values.
#' @param lik_init Initialization object returned by
#'   \code{coal_lik_init()}.
#' @param bound Numeric upper bound on the TMRCA.
#' @param eps Convergence tolerance.
#' @param eta Step size for gradient ascent.
#'
#' @return A list containing the final parameter vector and the sequence of
#'   log-likelihood values.
n_e_gradient_ascent <- function(f_init, lik_init, bound, eps, eta) {
  diff <- 1
  num_steps <- 1

  curr_f <- f_init
  ll <- c()
  while (diff > eps) {
    result <- coal_loglik_bounded(lik_init, curr_f)
    ll <- c(ll, result$ll)
    new_f <- curr_f + eta * result$dll

    diff <- sum(abs(new_f - curr_f))
    if (num_steps %% 10 == 0) {
      print(new_f)
      print(diff)
    }
    num_steps <- num_steps + 1
    print(num_steps)
    curr_f <- new_f
    print(new_f)
    print(diff)
    print(result$ll)
  }

  list(new_f, ll)
}


#' Bounded skyline estimate of effective population size
#'
#' Computes a piecewise-constant estimate of effective population size under
#' the bounded coalescent.
#'
#' @param data A \code{phylo} object or a list/data.frame with
#'   \code{coal_times}, \code{samp_times}, and \code{n_sampled}.
#' @param bound Numeric upper bound on the TMRCA.
#'
#' @return A list with components \code{Ne} and \code{grid}.
#'
#' @export
bounded_skyline_ascent <- function(data, bound = 1) {
  if (inherits(data, "phylo")) {
    phy <- phylodyn::summarize_phylo(data)
  } else if (all(c("coal_times", "samp_times", "n_sampled") %in%
                 names(data))) {
    phy <- with(data, list(
      samp_times = samp_times, coal_times = coal_times,
      n_sampled = n_sampled
    ))
  } else {
    stop("data must be a phylo or a
  list/data.frame with coal_times, samp_times, n_sampled")
  }
  grid <- c(0, data$coal_times, bound + 1e-04)
  lik_init <- coal_lik_init(
    samp_times = phy$samp_times,
    n_sampled = phy$n_sampled, coal_times = phy$coal_times,
    grid = grid
  )
  f_init <- rep(0, lik_init$ng)
  eta <- 0.01
  eps <- 0.02

  fit <- n_e_gradient_ascent(f_init, lik_init, bound, eps, eta)

  n_e <- exp(c(fit[[1]][1], fit[[1]]))
  l <- length(n_e)
  list(Ne = n_e[-l], grid = grid[-l])
}


#' Simulate from an Inhomogeneous Bounded Coalescent
#'
#' Simulates coalescent times under an inhomogeneous bounded coalescent model
#' with contemporaneous or heterochronous sampling. The bound is interpreted as
#' an upper bound on the time to the most recent common ancestor (TMRCA).
#'
#' At present, this implementation is intended for contemporaneous sampling.
#' The function may not be correct for heterochronous sampling because the
#' `r_j` factors have not yet been updated for that case.
#'
#' @param samp_times Numeric vector of sampling times.
#' @param n_sampled Numeric vector giving the number of samples taken at each
#'   sampling time in `samp_times`.
#' @param traj Function returning the effective population size at time `t`.
#'   The default is `constant`.
#' @param bound Numeric scalar giving the upper bound on the TMRCA.
#' @param val_upper Numeric scalar giving the upper value used for root-finding
#'   in the transformation method.
#' @param ... Additional arguments passed to `traj` when `traj` is not a step
#'   function.
#'
#' @return A list with components:
#' \describe{
#'   \item{coal_times}{Numeric vector of coalescent times.}
#'   \item{lineages}{Numeric vector giving the number of active lineages just
#'   before each coalescent event.}
#'   \item{intercoal_times}{Numeric vector of intercoalescent times.}
#'   \item{samp_times}{The input sampling times.}
#'   \item{n_sampled}{The input numbers sampled at each sampling time.}
#' }
#'
#' @export
coalsim_bounded <- function(samp_times, n_sampled, traj, bound = 1, val_upper = 10, ...) {
  # This function is not functional for heterochronous sampling
  # the r_j factors need to be updated
  ## hazard and inv. hazard for dominating intensity
  if (stats::is.stepfun(traj)) {
    knots <- knots(traj)
    midpts <- c(
      min(knots) - 1, knots[-1] - diff(knots) / 2,
      max(knots) + 1
    )
    traj_inv <- stats::stepfun(x = knots, y = 1 / traj(midpts))
    hazard_simple <- function(t, start, target) integrate_step_fun(traj_inv, start, start + t) - target
    hazard_simple_bound <- hazard_simple(bound, 0, 0)
    is_stepfun <- TRUE
  } else {
    traj_inv <- function(t) 1 / traj(t, ...)
    hazard_simple <- function(t, start, target) stats::integrate(traj_inv, start, start + t)$value - target
    hazard_simple_bound <- hazard_simple(bound, 0, 0)
    is_stepfun <- FALSE
  }
  val_upper <- min(10, 2 * traj_inv(bound))

  ## hazard target
  r_func <- function(k, j) {
    if (j == 1) {
      return(1)
    }
    prod <- 1
    for (m in 1:(j - 1)) {
      prod <- prod * ((2 * m + 1) / (2 * m - 1)) * ((k - m) / (k + m))
    }
    return((-1)^(j - 1) * prod)
  }
  ntip <- sum(n_sampled)
  result_list <- lapply(seq_len(ntip), function(k) {
    sapply(seq_len(k), function(i) r_func(k, i))
  })
  com_vec <- choose(1:ntip, 2)

  coal_times <- NULL
  lineages <- NULL
  curr <- 1
  active_lineages <- n_sampled[curr]
  time <- samp_times[curr]
  while (time <= max(samp_times) || active_lineages > 1) {
    if (active_lineages == 1) {
      curr <- curr + 1
      active_lineages <- active_lineages + n_sampled[curr]
      time <- samp_times[curr]
    }
    w <- stats::rexp(1) / (0.5 * active_lineages * (active_lineages - 1))
    target <- hazard_simple_bound - log(1 + exp(-w) * (exp(hazard_simple_bound - hazard_simple(time, 0, 0)) - 1))
    if (is_stepfun) {
      y <- hazard_uniroot_stepfun(
        traj_inv_stepfun = hazard_simple,
        start = 0, target = target
      )
    } else {
      y <- stats::uniroot(hazard_simple,
                          start = 0, target = target, lower = 0, upper = val_upper,
                          extendInt = "upX"
      )$root
    }
    while (curr < length(samp_times) && y >= samp_times[curr +
                                                        1]) {
      curr <- curr + 1
      active_lineages <- active_lineages + n_sampled[curr]
      time <- samp_times[curr]
      w <- stats::rexp(1) / (0.5 * active_lineages * (active_lineages - 1))
      target <- hazard_simple_bound - log(1 + exp(-w) * (exp(hazard_simple_bound - hazard_simple(time, 0, 0)) - 1))
      if (is_stepfun) {
        y <- hazard_uniroot_stepfun(
          traj_inv_stepfun = hazard_simple,
          start = 0, target = target
        )
      } else {
        y <- stats::uniroot(hazard_simple,
                            start = 0, target = target, lower = 0, upper = val_upper,
                            extendInt = "upX"
        )$root
      }
    }
    time <- y
    ## now thinning here
    con <- hazard_simple(time, 0, 0) - hazard_simple_bound
    denom <- sum(result_list[[active_lineages]] * exp(com_vec[1:(active_lineages)] * con))
    numer <- sum(result_list[[active_lineages - 1]] * exp(com_vec[1:(active_lineages - 1)] * con))
    if (stats::runif(1) <= (numer / denom) * (1 - exp(con))) {
      coal_times <- c(coal_times, time)
      lineages <- c(lineages, active_lineages)
      active_lineages <- active_lineages - 1
    }
  }
  return(list(
    coal_times = coal_times, lineages = lineages,
    intercoal_times = c(coal_times[1], diff(coal_times)),
    samp_times = samp_times, n_sampled = n_sampled
  ))
}


#' Initialize Coalescent Likelihood Quantities
#'
#' Constructs the interval-level quantities needed to evaluate a piecewise
#' coalescent log-likelihood on a specified time grid.
#'
#' @param samp_times Numeric vector of sampling times.
#' @param n_sampled Numeric vector giving the number of samples taken at each
#'   sampling time in `samp_times`.
#' @param coal_times Numeric vector of coalescent times.
#' @param grid Numeric vector of grid breakpoints defining the piecewise
#'   intervals for the effective population size trajectory.
#'
#' @return A list containing interval-level quantities used in the coalescent
#'   likelihood calculation, including active lineage counts, interval lengths,
#'   coalescent indicators, grid replication counts, and bookkeeping indices.
#'
#' @examples
#' samp_times <- 0
#' n_sampled <- 5
#' coal_times <- c(0.1, 0.2, 0.4, 0.8)
#' grid <- c(0, 0.2, 0.5, 1)
#'
#' coal_lik_init(samp_times, n_sampled, coal_times, grid)
#'
#' @export
coal_lik_init = function(samp_times, n_sampled, coal_times, grid)
{
  ns = length(samp_times)
  nc = length(coal_times)
  ng = length(grid)-1


  if (length(samp_times) != length(n_sampled))
    stop("samp_times vector of differing length than n_sampled vector.")

  if (length(coal_times) != sum(n_sampled) - 1)
    stop("Incorrect length of coal_times: should be sum(n_sampled) - 1.")

  if (max(samp_times, coal_times) > max(grid))
    stop("Grid does not envelop all sampling and/or coalescent times.")

  t = sort(unique(c(samp_times, coal_times, grid)))
  l = rep(0, length(t))

  for (i in 1:ns)
    l[t >= samp_times[i]] = l[t >= samp_times[i]] + n_sampled[i]

  for (i in 1:nc)
    l[t >= coal_times[i]] = l[t >= coal_times[i]] - 1

  #print(l)

  if (sum((l < 1) & (t >= min(samp_times))) > 0)
    stop("Number of active lineages falls below 1 after the first sampling point.")

  mask = l > 0
  t = t[mask]
  l = utils::head(l[mask], -1)

  gridrep = rep(0, ng)
  for (i in 1:ng)
    gridrep[i] = sum(t > grid[i] & t <= grid[i+1])

  C = 0.5 * l * (l-1)
  D = diff(t)

  y = rep(0, length(D))
  y[t[-1] %in% coal_times] = 1

  bins = cut(x = samp_times, breaks = t,
             include.lowest = TRUE)
  tab <- stats::aggregate(n_sampled ~ bins, FUN = sum, labels = FALSE)
  count <- rep(0, length(D))
  count[as.numeric(tab$bins)] <- tab$n_sampled
  count[utils::head(t, -1) >= max(samp_times)] <- NA

  rep_idx = cumsum(gridrep)
  rep_idx = cbind(rep_idx-gridrep+1,rep_idx)

  return(list(t=t, l=l, C=C, D=D, y=y, count=count, gridrep=gridrep, ns=sum(n_sampled), nc=nc, ng=ng, rep_idx=rep_idx, args=list(samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times, grid=grid)))
}

#' Integrate a Step Function Over an Interval
#'
#' Computes the integral of a step function over the interval from `from` to
#' `to` by summing the areas of the constant segments.
#'
#' @param fun A step function, typically created with [stats::stepfun()].
#' @param from Numeric scalar giving the lower limit of integration.
#' @param to Numeric scalar giving the upper limit of integration.
#'
#' @return A numeric scalar giving the value of the integral of `fun` from
#'   `from` to `to`.
integrate_step_fun = function(fun, from, to) {
  breaks = stats::knots(fun)

  mask = breaks > from & breaks < to
  pts = c(from, breaks[mask], to)
  diffs = diff(pts)
  midpts = pts[-1] - diffs/2

  segments = diffs * fun(midpts)

  return(sum(segments))
}


#' Invert a Step-Function Hazard by Root Finding
#'
#' Finds the elapsed time needed for the cumulative hazard of a step-function
#' trajectory to reach a specified target value, given a fixed number of active
#' lineages.
#'
#' @param traj_inv_stepfun A step function giving the inverse effective
#'   population size trajectory, typically created with [stats::stepfun()].
#' @param lineages Numeric scalar giving the number of active lineages.
#' @param start Numeric scalar giving the starting time.
#' @param target Numeric scalar giving the target cumulative hazard to reach.
#'
#' @return A numeric scalar giving the elapsed time from `start` until the
#'   cumulative hazard reaches `target`.
hazard_uniroot_stepfun <- function(traj_inv_stepfun, lineages, start, target) {
  knots = knots(traj_inv_stepfun)
  lin_factor = 0.5 * lineages * (lineages - 1)
  t = start

  while (target > 0 && sum(knots > t) > 0)
  {
    next_knot = min(knots[knots > t])
    if ((next_knot - t) * lin_factor * traj_inv_stepfun(mean(c(t, next_knot))) > target)
    {
      result = t + target / (lin_factor * traj_inv_stepfun(mean(c(t, next_knot))))
      target = 0
    }
    else
    {
      target = target - (next_knot - t) * lin_factor * traj_inv_stepfun(mean(c(t, next_knot)))
      t = next_knot
    }
  }
  if (sum(knots > t) < 1)
    result = t + target / (lin_factor * traj_inv_stepfun(t + 1))

  return(result - start)
}
