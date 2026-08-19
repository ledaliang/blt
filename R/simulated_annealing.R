#' Adjacent Swap Proposal on a Tree Topology
#'
#' Proposes a new topology by swapping one of the two children of a coalescent
#' event with one of the two children of the next event back in time. This is
#' the proposal kernel used by the simulated annealing stage of BLT.
#'
#' @details
#' A direct port of \code{adjacent_swap()} in
#' \code{new_code/pipeline/sa_topology_fast.R}. One of the first \eqn{n - 2}
#' rows of the matching is chosen at random, and a child of that row is
#' exchanged with a child of the row above it. Only children whose labels are
#' smaller than the parent of the row above can be moved, so the result is
#' always a valid matching.
#'
#' @param pairs Matching representation of the topology (output of
#' \code{matching}).
#'
#' @returns A matching of the same shape as \code{pairs}, with one swap
#' applied.
#'
#' @export
adjacent_swap <- function(pairs) {
  n <- nrow(pairs) + 1
  r <- sample.int(n - 2, 1)
  pair1 <- pairs[r, 1:2]
  new_pairs <- pairs
  possible <- which(pair1 < pairs[r + 1, 3])
  s1 <- if (length(possible) == 1) possible else sample(possible, 1)
  s2 <- sample.int(2, 1)
  new_pairs[r, s1] <- pairs[r + 1, s2]
  new_pairs[r + 1, s2] <- pairs[r, s1]
  new_pairs
}


#' Simulated Annealing Acceptance Probability
#'
#' The Metropolis acceptance probability for a move which is maximizing the
#' log-likelihood. A proposal that improves the log-likelihood is always
#' accepted; one that lowers it is accepted with probability
#' \eqn{\exp((e_2 - e_1 + \epsilon) / T)}.
#'
#' @param e1 Log-likelihood of the current state.
#' @param e2 Log-likelihood of the proposed state.
#' @param temp Current temperature.
#' @param eps Small constant added to the difference, as in the original code.
#'
#' @returns The acceptance probability.
#'
#' @export
sa_accept_prob <- function(e1, e2, temp, eps = 1e-4) {
  if (e2 > e1) {
    return(1)
  }
  exp((e2 - e1 + eps) / temp)
}


#' Simulated Annealing Cooling Schedule
#'
#' @param schedule One of \code{"exp"} (geometric cooling by \code{alpha}),
#' \code{"log"}, or \code{"lin"}.
#' @param init_temp Initial temperature.
#' @param alpha Cooling factor, used by the exponential schedule.
#' @param temp Current temperature.
#' @param i Index of the current proposal.
#'
#' @returns The temperature for the next proposal.
#'
#' @export
sa_temp_update <- function(schedule, init_temp, alpha, temp, i) {
  switch(schedule,
    exp = alpha * temp,
    log = init_temp / (1 + log(1 + i)),
    lin = init_temp / (1 + i),
    stop("Unknown schedule: ", schedule)
  )
}


#' Build an ape Tree from a Matching and Coalescent Times
#'
#' Turns the internal representation used throughout the package, a matching
#' plus a vector of coalescent times measured from the present, into an
#' \code{ape} tree.
#'
#' @param pairs Matching representation of the topology (output of
#' \code{matching}).
#' @param coal_times Coalescent times measured from the leaves, of length
#' \eqn{n - 1}, in the order of the rows of \code{pairs}.
#'
#' @returns An object of class \code{phylo}.
#'
#' @export
tree_from_pairs <- function(pairs, coal_times) {
  n <- nrow(pairs) + 1L
  edges <- matrix(0L, nrow = 2 * n - 2, ncol = 2)
  edge_lengths <- rep(0, 2 * n - 2)
  for (r in seq_len(nrow(pairs))) {
    edges[2 * r - 1, 1] <- pairs[r, 3]
    edges[2 * r - 1, 2] <- pairs[r, 1]
    edge_lengths[2 * r - 1] <- if (pairs[r, 1] <= n) {
      coal_times[pairs[r, 3] - n]
    } else {
      coal_times[pairs[r, 3] - n] - coal_times[pairs[r, 1] - n]
    }
    edges[2 * r, 1] <- pairs[r, 3]
    edges[2 * r, 2] <- pairs[r, 2]
    edge_lengths[2 * r] <- if (pairs[r, 2] <= n) {
      coal_times[pairs[r, 3] - n]
    } else {
      coal_times[pairs[r, 3] - n] - coal_times[pairs[r, 2] - n]
    }
  }
  # relabel the inner nodes so that the root is n + 1, the ape convention
  for (r in seq_len(nrow(edges))) {
    if (edges[r, 1] > n) edges[r, 1] <- 2 * n - (edges[r, 1] - n)
    if (edges[r, 2] > n) edges[r, 2] <- 2 * n - (edges[r, 2] - n)
  }
  tree <- list(
    edge = edges, edge.length = edge_lengths,
    tip.label = paste0("t", seq_len(n)), Nnode = n - 1L
  )
  class(tree) <- "phylo"
  tree
}


#' Posterior Mean Mutation Rates Given Intercoalescent Times
#'
#' Refits the mutation rate matrix from the observed mutation counts and the
#' current intercoalescent times, under the model
#' \eqn{\theta \sim \mathrm{Dirichlet}(\alpha)},
#' \eqn{\lambda \sim \mathrm{Exp}(\beta)},
#' \eqn{C_i \mid \theta_i \sim \mathrm{Poisson}(\lambda L \theta_i)}, where
#' \eqn{L} is the total branch length of the tree. Because the number of
#' lineages in each interval of a coalescent tree is \eqn{n, n-1, \ldots, 2}
#' whatever the topology, \eqn{L} depends on the times but not on the topology.
#'
#' @param c_vec Observed mutation counts, flattened with
#' \code{flatten_upper}.
#' @param ict Intercoalescent times.
#' @param time_obs Duration of the experiment. The tree is rescaled to this
#' height before the branch lengths are totalled.
#' @param alpha_prior Dirichlet prior parameters for the relative rates.
#' @param beta_prior Rate of the exponential prior on the overall rate.
#'
#' @returns A list with the rate matrix \code{theta}, the overall rate
#' \code{lambda}, and the \code{total_branch_length} used.
#'
#' @export
refit_theta <- function(c_vec, ict, time_obs, alpha_prior, beta_prior = 1) {
  n <- length(ict) + 1L
  k <- n:2
  total_branch_length <- sum(k * ict) * (time_obs / sum(ict))
  post_mean <- (alpha_prior + c_vec) / (sum(c_vec) + sum(alpha_prior))
  lambda <- (sum(c_vec) + 1) / (beta_prior + total_branch_length)
  list(
    theta = unflatten_upper(lambda * post_mean),
    lambda = lambda,
    total_branch_length = total_branch_length
  )
}


#' Number of Alleles Seen at Each Site or Span
#'
#' Counts how many distinct alleles the data contains at each site and at each
#' span of sites cut simultaneously, and returns the matrix \code{m} of allele
#' counts used by the likelihood. Every entry starts at \code{m_default} and is
#' raised to the observed count plus \code{m_headroom} where that is larger, so
#' that alleles which were never observed still have some probability.
#'
#' @param d Matrix of observed alleles at the leaves (\eqn{n \times S}).
#' @param s Number of sites.
#' @param m_default Value every entry starts at.
#' @param m_headroom Number of unobserved alleles allowed on top of the
#' observed count.
#'
#' @returns An \eqn{S \times S} matrix of allele counts.
#'
#' @export
allele_span_counts <- function(d, s, m_default = 20, m_headroom = 5) {
  counts <- matrix(0L, nrow = s, ncol = s)
  seen <- list()
  for (i in seq_len(nrow(d))) {
    for (j in seq_len(s)) {
      value <- d[i, j]
      if (value == "0") {
        next
      }
      prefix <- as.integer(strsplit(value, ":")[[1]][1])
      if (prefix <= s) {
        from <- j
        to <- j
      } else {
        from <- prefix %/% 10
        to <- prefix %% 10
      }
      key <- paste(from, to, value, sep = "|")
      if (is.null(seen[[key]])) {
        seen[[key]] <- TRUE
        counts[from, to] <- counts[from, to] + 1L
      }
    }
  }
  m <- matrix(m_default, nrow = s, ncol = s)
  raised <- counts + m_headroom
  m[counts > 0] <- pmax(m[counts > 0], raised[counts > 0])
  m
}


#' Simulated Annealing over Tree Topologies
#'
#' Searches over topologies with the \code{adjacent_swap} proposal. Each
#' proposal is scored by refitting the intercoalescent times on the proposed
#' topology for a small number of gradient ascent steps and taking the
#' resulting log-likelihood, and is then accepted or rejected by the Metropolis
#' rule. The best topology and times seen along the way are returned.
#'
#' @details
#' A port of \code{sa_topology_swap_fast()} in
#' \code{new_code/pipeline/sa_topology_fast.R}, with the fast rescaled
#' likelihood replaced by \code{estimate_coal_times}.
#'
#' @param d Matrix of observed alleles at the leaves, or a list of such
#' matrices, one per integration barcode.
#' @param hap Output of \code{haplotypes} for \code{d}, or a list of such
#' outputs.
#' @param curr_loglik Log-likelihood of the current topology and times.
#' @param curr_pairs Matching representation of the current topology.
#' @param ict Current intercoalescent times.
#' @param eta Initial learning rate for the inner gradient ascent.
#' @param eps Stopping criterion for the inner gradient ascent.
#' @param q Rate matrix for the mutation states.
#' @param eigen_q Eigendecomposition of the rate matrix, with the extra
#' component \code{inv}.
#' @param states Hash table of mutation states.
#' @param states_matrix Matrix of mutation states.
#' @param m Matrix of allele counts.
#' @param approx Set to \code{0} for the exact likelihood and \code{1} for the
#' approximate likelihood.
#' @param t_prob Transition probability matrix over one time unit.
#' @param total_time Height of the tree.
#' @param inner_max_steps Gradient ascent steps allowed per proposal.
#' @param schedule Cooling schedule, one of \code{"exp"}, \code{"lin"},
#' \code{"log"}.
#' @param init_temp Initial temperature.
#' @param alpha Cooling factor.
#' @param max_iter Number of proposals.
#' @param verbose If \code{TRUE}, print one line per proposal.
#'
#' @returns A list with the best \code{pairs}, \code{ict} and \code{loglik}
#' found, and the number of proposals accepted, \code{n_accept}.
#'
#' @export
sa_topology_swap <- function(
  d,
  hap,
  curr_loglik,
  curr_pairs,
  ict,
  eta,
  eps,
  q,
  eigen_q,
  states,
  states_matrix,
  m,
  approx,
  t_prob,
  total_time = 10,
  inner_max_steps = 20,
  schedule = c("exp", "lin", "log"),
  init_temp = 1,
  alpha = 0.95,
  max_iter = 30,
  verbose = TRUE
) {
  schedule <- match.arg(schedule)
  config <- curr_pairs
  curr_ict <- ict
  e1 <- curr_loglik
  best_energy <- e1
  best_config <- config
  best_ict <- curr_ict
  temp <- init_temp
  n_accept <- 0L

  for (i in seq_len(max_iter)) {
    new_pairs <- adjacent_swap(config)
    inner <- estimate_coal_times(
      pairs = new_pairs, d = d, hap = hap, ict = curr_ict,
      eta = eta, eps = eps, q = q, eigen_q = eigen_q,
      states = states, states_matrix = states_matrix, m = m,
      t_prob = t_prob, total_time = total_time, approx = approx,
      max_steps = inner_max_steps
    )
    new_ict <- inner$ict
    e2 <- inner$logliks[length(inner$logliks)]

    prob <- sa_accept_prob(e1, e2, temp)
    accepted <- stats::runif(1) <= prob
    if (accepted) {
      config <- new_pairs
      curr_ict <- new_ict
      e1 <- e2
      n_accept <- n_accept + 1L
    }
    if (e1 > best_energy) {
      best_config <- config
      best_ict <- curr_ict
      best_energy <- e1
    }
    if (verbose) {
      cat(sprintf(
        "    iter %2d  lik=%.4f  prop=%.4f  prob=%.3g  temp=%.3g  %s\n",
        i, e1, e2, prob, temp, if (accepted) "ACC" else "rej"
      ))
    }
    temp <- sa_temp_update(schedule, init_temp, alpha, temp, i)
  }
  if (verbose) {
    cat(sprintf(
      "  [SA] accepted %d/%d, best loglik = %.4f\n",
      n_accept, max_iter, best_energy
    ))
  }
  list(
    pairs = best_config, ict = best_ict, loglik = best_energy,
    n_accept = n_accept
  )
}


#' Bayesian Lineage Tracing
#'
#' Estimates the tree topology, the coalescent times and the mutation rates
#' from observed alleles, by alternating simulated annealing over topologies
#' with refits of the rates and of the times.
#'
#' @details
#' The algorithm, referred to as BLT and as \code{sa_v1} in
#' \code{new_code/pipeline}, runs as follows.
#' \enumerate{
#'   \item A UPGMA tree on the Hamming distances between the lumped mutation
#'   states gives the starting topology and the starting intercoalescent
#'   times.
#'   \item The mutation rates are estimated from the observed mutation counts
#'   and the total branch length of that tree, by \code{refit_theta}.
#'   \item The intercoalescent times are fitted on the starting topology by
#'   \code{estimate_coal_times}, for at most \code{max_steps} steps.
#'   \item Then for \code{max_outer_iter} cycles: simulated annealing over
#'   topologies with \code{sa_max_iter} proposals, each scored by
#'   \code{inner_sa_steps} gradient ascent steps; the rates are refitted from
#'   the new times; and the times are refitted on the new topology for at most
#'   \code{inner_max_steps} steps.
#' }
#' The height of the tree is held at the height of the UPGMA tree throughout,
#' since the likelihood on its own cannot identify the overall scale.
#'
#' @param d Matrix of observed alleles at the leaves (\eqn{n \times S}), with
#' \code{"0"} for an uncut site.
#' @param s Number of sites.
#' @param time Duration of the experiment, used when refitting the rates.
#' @param eta Initial learning rate for gradient ascent.
#' @param eps Stopping criterion for gradient ascent.
#' @param max_steps Maximum steps for the initial gradient ascent.
#' @param inner_max_steps Maximum steps when refitting times at the end of a
#' cycle.
#' @param max_outer_iter Number of cycles.
#' @param sa_max_iter Number of proposals per simulated annealing stage.
#' @param inner_sa_steps Gradient ascent steps used to score a proposal.
#' @param sa_init_temp Initial temperature.
#' @param sa_alpha Cooling factor.
#' @param sa_schedule Cooling schedule.
#' @param alpha_prior Dirichlet prior for the relative rates. Defaults to
#' \code{rep(5, s * (s + 1) / 2)}, which keeps every rate away from zero.
#' @param beta_prior Rate of the exponential prior on the overall rate.
#' @param m_default Starting number of alleles per site or span.
#' @param m_headroom Unobserved alleles allowed on top of the observed count.
#' @param approx Set to \code{0} for the exact likelihood and \code{1} for the
#' approximate likelihood.
#' @param verbose If \code{TRUE}, report progress.
#'
#' @returns A list with components
#' \describe{
#'   \item{\code{tree}}{The estimated tree, of class \code{phylo}.}
#'   \item{\code{pairs}}{The estimated topology as a matching.}
#'   \item{\code{ict}}{The estimated intercoalescent times.}
#'   \item{\code{coal_times}}{The estimated coalescent times from the leaves.}
#'   \item{\code{theta}}{The estimated mutation rate matrix.}
#'   \item{\code{lambda}}{The estimated overall mutation rate.}
#'   \item{\code{upgma_tree}}{The starting UPGMA tree.}
#'   \item{\code{d}}{The allele matrix, with rows in tip order.}
#'   \item{\code{loglik_trace}}{Log-likelihood after the initial fit and after
#'   each cycle.}
#'   \item{\code{n_accept}}{Proposals accepted in each cycle.}
#' }
#'
#' @export
run_blt <- function(
  d,
  s,
  time,
  eta = 1e-3,
  eps = 1e-3,
  max_steps = 200,
  inner_max_steps = 100,
  max_outer_iter = 2,
  sa_max_iter = 30,
  inner_sa_steps = 20,
  sa_init_temp = 1,
  sa_alpha = 0.95,
  sa_schedule = c("exp", "lin", "log"),
  alpha_prior = NULL,
  beta_prior = 1,
  m_default = 20,
  m_headroom = 5,
  approx = 1,
  verbose = TRUE
) {
  sa_schedule <- match.arg(sa_schedule)
  if (is.null(alpha_prior)) alpha_prior <- rep(5, s * (s + 1) / 2)
  n <- nrow(d)
  if (is.null(rownames(d))) rownames(d) <- paste0("cell_", seq_len(n))
  cell_names <- rownames(d)

  # ---- starting topology and times from UPGMA ----
  lumped_d <- t(apply(d, 1, lumped_state))
  dist_matrix <- hamming_dist_matrix(lumped_d)
  dimnames(dist_matrix) <- list(cell_names, cell_names)
  upgma_tree <- phangorn::upgma(stats::as.dist(dist_matrix))

  # the likelihood indexes the data by tip number, so put the rows in tip order
  leaf_to_row <- match(upgma_tree$tip.label, cell_names)
  d <- d[leaf_to_row, , drop = FALSE]
  lumped_d <- lumped_d[leaf_to_row, , drop = FALSE]

  pairs <- matching(upgma_tree)
  intervals <- ape::coalescent.intervals(upgma_tree)
  tree_height <- intervals$total.depth
  ict <- pmax(intervals$interval.length, 1e-4)
  ict <- ict * (tree_height / sum(ict))
  if (verbose) {
    cat("[init] UPGMA tree:", n, "tips, height =", tree_height, "\n")
  }

  # ---- starting mutation rates ----
  c_vec <- flatten_upper(observed_mutations(d, s, lumped_d)$obsM)
  theta_fit <- refit_theta(c_vec, ict, time, alpha_prior, beta_prior)
  theta <- theta_fit$theta

  # ---- substitution model ----
  states <- state_space(s)
  states_matrix <- state_space_matrix(s)
  m <- allele_span_counts(d, s, m_default, m_headroom)
  hap <- haplotypes(d, n, s)
  build_sub_model <- function(theta) {
    q <- q_matrix(states_matrix, theta)
    eigen_q <- eigen(q)
    eigen_q$inv <- solve(eigen_q$vectors)
    t_prob <- eigen_q$vectors %*% diag(exp(eigen_q$values)) %*% eigen_q$inv
    list(q = list(q), eigen_q = list(eigen_q), t_prob = list(t_prob))
  }
  sub <- build_sub_model(theta)
  if (verbose) {
    cat("[init] lambda =", theta_fit$lambda,
        "  states:", nrow(states_matrix),
        "  haplotypes:", nrow(hap$alleles), "\n")
  }

  # ---- initial fit of the times ----
  if (verbose) cat("[init] gradient ascent on the intercoalescent times\n")
  fit <- estimate_coal_times(
    pairs = pairs, d = list(d), hap = list(hap), ict = ict,
    eta = eta, eps = eps, q = sub$q, eigen_q = sub$eigen_q,
    states = states, states_matrix = states_matrix, m = m,
    t_prob = sub$t_prob, total_time = tree_height, approx = approx,
    max_steps = max_steps
  )
  ict <- fit$ict
  curr_loglik <- fit$logliks[length(fit$logliks)]
  loglik_trace <- curr_loglik
  n_accept <- integer(0)
  if (verbose) cat("[init] loglik =", curr_loglik, "\n")

  # ---- cycles of simulated annealing, rates, times ----
  for (outer in seq_len(max_outer_iter)) {
    if (verbose) {
      cat(sprintf("\n[cycle %d/%d] simulated annealing over topologies\n",
                  outer, max_outer_iter))
    }
    sa_res <- sa_topology_swap(
      d = list(d), hap = list(hap), curr_loglik = curr_loglik,
      curr_pairs = pairs, ict = ict, eta = eta, eps = eps,
      q = sub$q, eigen_q = sub$eigen_q, states = states,
      states_matrix = states_matrix, m = m, approx = approx,
      t_prob = sub$t_prob, total_time = tree_height,
      inner_max_steps = inner_sa_steps, schedule = sa_schedule,
      init_temp = sa_init_temp, alpha = sa_alpha, max_iter = sa_max_iter,
      verbose = verbose
    )
    pairs <- sa_res$pairs
    ict <- sa_res$ict
    n_accept <- c(n_accept, sa_res$n_accept)

    theta_fit <- refit_theta(c_vec, ict, time, alpha_prior, beta_prior)
    theta <- theta_fit$theta
    sub <- build_sub_model(theta)
    if (verbose) cat(sprintf("[cycle %d] new lambda = %g\n", outer,
                             theta_fit$lambda))

    fit <- estimate_coal_times(
      pairs = pairs, d = list(d), hap = list(hap), ict = ict,
      eta = eta, eps = eps, q = sub$q, eigen_q = sub$eigen_q,
      states = states, states_matrix = states_matrix, m = m,
      t_prob = sub$t_prob, total_time = tree_height, approx = approx,
      max_steps = inner_max_steps
    )
    ict <- fit$ict
    curr_loglik <- fit$logliks[length(fit$logliks)]
    loglik_trace <- c(loglik_trace, curr_loglik)
    if (verbose) cat(sprintf("[cycle %d] loglik = %.4f\n", outer, curr_loglik))
  }

  coal_times <- cumsum(ict)
  tree <- tree_from_pairs(pairs, coal_times)
  tree$tip.label <- rownames(d)

  list(
    tree = tree, pairs = pairs, ict = ict, coal_times = coal_times,
    theta = theta, lambda = theta_fit$lambda, upgma_tree = upgma_tree,
    d = d, loglik_trace = loglik_trace, n_accept = n_accept
  )
}
