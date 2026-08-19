#' Project a Vector onto the L1 Ball
#'
#' Computes the Euclidean projection of a vector onto the L1 ball of a given
#' radius. This is used to keep the vector of intercoalescent times inside the
#' constraint set \eqn{\sum_i u_i \le T}, where \eqn{T} is the height of the
#' tree.
#'
#' Adapted from \url{https://rdrr.io/github/vguillemot/sparseMCA/src/R/projl1.R}
#'
#' @param v A numeric vector.
#' @param a The radius (> 0) of the L1 ball.
#'
#' @return The projection of \code{v} onto the L1 ball of radius \code{a}.
#'
#' @examples
#' proj_l1(1:10, 21)
#'
#' @export
proj_l1 <- function(v, a) {
  u <- sort(abs(v), decreasing = TRUE)
  n <- length(v)
  ukmaok <- (cumsum(u) - a) / (1:n)
  k <- max(which(ukmaok < u))
  tau <- ukmaok[k]

  sign(v) * pmax(abs(v) - tau, 0)
}


#' Rate Matrix Restricted to a Set of Alleles
#'
#' From a matrix of full allele states, builds the rate matrix specifically for
#' those alleles. Note that this is not a complete rate matrix (the rows will
#' not sum to zero), it only stores the rates for the alleles which have
#' non-zero probability in the likelihood.
#'
#' @details
#' This is a direct port of \code{Q_alleles()} in
#' \code{new_code/new_functions.R}. The only change is that the rate matrix is
#' passed in as \code{q} rather than picked up from the global environment.
#' The \code{eigenQ} argument of the original is dropped because it is unused
#' there.
#'
#' @param alleles A matrix of allele states, one per row (the first element of
#' the output of \code{pq_vectors}).
#' @param q Rate matrix for the mutation states (output of \code{q_matrix}).
#' @param states A hash table of mutation states (output of
#' \code{state_space}).
#' @param m Symmetric matrix which stores the number of alleles per site or
#' site group.
#' @param t_prob Transition probability matrix over one time unit.
#'
#' @returns An \eqn{A \times A} matrix, where \eqn{A} is the number of rows of
#' \code{alleles}, holding the mutation state rate for every pair of alleles
#' with a non-zero transition probability.
#'
#' @export
q_alleles <- function(alleles, q, states, m, t_prob) {
  a <- nrow(alleles)
  q_new <- matrix(0, nrow = a, ncol = a)

  for (i in 1:a) {
    allele1 <- alleles[i, ]
    for (j in 1:a) {
      allele2 <- alleles[j, ]
      if (transition_prob_finite_alleles(
        allele1, allele2, t_prob, states, m
      ) > 0) {
        lumped1 <- lumped_state(allele1)
        lumped2 <- lumped_state(allele2)
        i1 <- states[[paste(lumped1, collapse = "")]]
        i2 <- states[[paste(lumped2, collapse = "")]]
        q_new[i, j] <- q[i1, i2]
      }
    }
  }
  q_new
}


#' Gradient of the Log-Likelihood with Respect to Branch Lengths
#'
#' Returns the vector of gradients of the log-likelihood with respect to the
#' branch lengths of the tree. \code{p} and \code{q} are the forward and
#' backward probability matrices returned by \code{pq_vectors}: each column
#' corresponds to a node, and the last column corresponds to the root, which
#' has no branch above it and is therefore dropped.
#'
#' @details
#' This is a direct port of \code{gradient()} in
#' \code{new_code/new_functions.R}.
#'
#' @param p Forward probability matrix (second element of the output of
#' \code{pq_vectors}).
#' @param q Backward probability matrix (third element of the output of
#' \code{pq_vectors}).
#' @param q_allele Rate matrix restricted to the observed alleles (output of
#' \code{q_alleles}).
#'
#' @returns A numeric vector of length \eqn{2n - 2}, where entry \eqn{j} is the
#' derivative of the log-likelihood with respect to the branch which ends at
#' node \eqn{j}.
#'
#' @export
gradient_branches <- function(p, q, q_allele) {
  like <- t(q[, 1]) %*% p[, 1]
  # exclude the last column, the root, which does not have a branch above it
  n_branches <- ncol(p) - 1
  g <- rep(0, n_branches)
  for (c in 1:n_branches) {
    g[c] <- t(q[, c]) %*% q_allele %*% p[, c]
  }

  g / c(like)
}


# Wraps single-barcode arguments into length-one lists, so that the rest of the
# code can always assume a list with one entry per integration barcode.
as_barcode_list <- function(x, k) {
  if (is.list(x) && length(x) == k && (is.list(x[[1]]) || is.matrix(x[[1]]))) {
    return(x)
  }
  replicate(k, x, simplify = FALSE)
}


# Runs pq_vectors for every integration barcode at the given branch lengths and
# returns the summed log-likelihood together with the forward and backward
# probability matrices.
pq_all <- function(
  pairs, d, hap, b, eigen_q, states, states_matrix, m, t_prob, approx, mu
) {
  k <- length(d)
  p_list <- list()
  q_list <- list()
  loglik <- 0
  for (i in seq_len(k)) {
    pq <- pq_vectors(
      d[[i]], hap[[i]], pairs, b, eigen_q[[i]], states, states_matrix, m,
      approx, t_prob[[i]],
      mu = if (is.null(mu)) NULL else mu[[i]]
    )
    p_list[[i]] <- pq[[2]]
    q_list[[i]] <- pq[[3]]
    loglik <- loglik + log(pq[[2]][1, ncol(pq[[2]])])
  }
  list(loglik = loglik, p = p_list, q = q_list)
}


#' Log-Likelihood of the Data Given Branch Lengths
#'
#' Computes the log-likelihood of the observed allele data on a tree with a
#' fixed topology and fixed branch lengths. When several independent
#' integration barcodes are observed on the same tree, the log-likelihoods are
#' summed.
#'
#' @param pairs Matching representation of the topology (output of
#' \code{matching}).
#' @param d Matrix of observed alleles at the leaves (\eqn{n \times S}), or a
#' list of such matrices, one per integration barcode.
#' @param hap Output of \code{haplotypes} for \code{d}, or a list of such
#' outputs.
#' @param b Numeric vector of branch lengths of length \eqn{2n - 2}, in the
#' node labelling used by \code{pairs}.
#' @param eigen_q Eigendecomposition of the rate matrix, with the extra
#' component \code{inv} holding the inverse of the eigenvectors. May be a list,
#' one per integration barcode.
#' @param states Hash table of mutation states (output of \code{state_space}).
#' @param states_matrix Matrix of mutation states (output of
#' \code{state_space_matrix}).
#' @param m Symmetric matrix which stores the number of alleles per site or
#' site group.
#' @param t_prob Transition probability matrix over one time unit. May be a
#' list, one per integration barcode.
#' @param approx Set to \code{0} for the exact likelihood and \code{1} for the
#' approximate likelihood which only considers alleles seen at the leaves.
#' @param mu Optional list of conditional mutation probabilities, one per
#' integration barcode. If \code{NULL} a uniform distribution is assumed.
#'
#' @returns The log-likelihood, summed over integration barcodes.
#'
#' @export
branch_loglik <- function(
  pairs,
  d,
  hap,
  b,
  eigen_q,
  states,
  states_matrix,
  m,
  t_prob,
  approx = 1,
  mu = NULL
) {
  if (!is.list(d)) {
    d <- list(d)
    hap <- list(hap)
  }
  k <- length(d)
  eigen_q <- as_barcode_list(eigen_q, k)
  t_prob <- as_barcode_list(t_prob, k)

  pq_all(
    pairs, d, hap, b, eigen_q, states, states_matrix, m, t_prob, approx, mu
  )$loglik
}


#' Estimate Coalescent Times by Gradient Ascent
#'
#' Maximizes the likelihood of the observed allele data over the vector of
#' intercoalescent times, holding the tree topology, the mutation rate matrix
#' and the number of alleles fixed. Only the branch lengths, that is the
#' coalescent times, are estimated.
#'
#' @details
#' Write \eqn{u = (u_1, \ldots, u_{n-1})} for the intercoalescent times, where
#' \eqn{u_k} is the length of the interval during which the tree has \eqn{k + 1}
#' lineages. Every branch length is a sum of intercoalescent times, so
#' \eqn{b = A u} where \eqn{A} is the \eqn{(2n - 2) \times (n - 1)} matrix
#' returned by \code{decompose_branches}. The chain rule then gives
#' \deqn{\frac{\partial \ell}{\partial u} = A^{T}
#'       \frac{\partial \ell}{\partial b},}
#' where the gradient with respect to the branch lengths is computed by
#' \code{gradient_branches}.
#'
#' The step size is chosen by AdaGrad. After every step the times are projected
#' back onto the L1 ball of radius \code{total_time} whenever they overshoot the
#' height of the tree, and the learning rate is repeatedly divided by five
#' whenever a step produces a non-positive time or decreases the
#' log-likelihood. Note that the learning rate is only ever decreased, never
#' restored.
#'
#' This is a direct port of \code{gradient_ascent_intercoaltimes_adagrad()} in
#' \code{new_code/new_functions.R} and returns the same times for the same
#' inputs. Two things differ outside the optimization itself: a single data
#' matrix is wrapped into a list of one rather than taking the separate
#' single-barcode code path of the original, and the return value is a named
#' list holding log-likelihoods instead of \code{list(ps, currICT)}.
#'
#' @param pairs Matching representation of the topology (output of
#' \code{matching}). The topology is held fixed.
#' @param d Matrix of observed alleles at the leaves (\eqn{n \times S}), or a
#' list of such matrices, one per integration barcode.
#' @param hap Output of \code{haplotypes} for \code{d}, or a list of such
#' outputs.
#' @param ict Numeric vector of length \eqn{n - 1} of starting intercoalescent
#' times.
#' @param eta Initial learning rate.
#' @param eps Stopping criterion on the change in log-likelihood.
#' @param q Rate matrix for the mutation states (output of \code{q_matrix}).
#' May be a list, one per integration barcode.
#' @param eigen_q Eigendecomposition of the rate matrix, with the extra
#' component \code{inv} holding the inverse of the eigenvectors. May be a list,
#' one per integration barcode.
#' @param states Hash table of mutation states (output of \code{state_space}).
#' @param states_matrix Matrix of mutation states (output of
#' \code{state_space_matrix}).
#' @param m Symmetric matrix which stores the number of alleles per site or
#' site group.
#' @param t_prob Transition probability matrix over one time unit. May be a
#' list, one per integration barcode.
#' @param total_time Height of the tree. The intercoalescent times are
#' constrained to sum to at most this value.
#' @param approx Set to \code{0} for the exact likelihood and \code{1} for the
#' approximate likelihood which only considers alleles seen at the leaves.
#' @param max_steps Maximum number of gradient ascent steps.
#' @param mu Optional list of conditional mutation probabilities, one per
#' integration barcode. If \code{NULL} a uniform distribution is assumed.
#' @param verbose If \code{TRUE}, print the current times and log-likelihood
#' every ten steps.
#'
#' @returns A list with components
#' \describe{
#'   \item{\code{ict}}{The estimated intercoalescent times, a numeric vector of
#' length \eqn{n - 1}.}
#'   \item{\code{coal_times}}{The estimated coalescent times measured from the
#' leaves, \code{cumsum(ict)}.}
#'   \item{\code{logliks}}{The log-likelihood at every accepted step.}
#'   \item{\code{steps}}{The number of steps taken.}
#' }
#'
#' @export
estimate_coal_times <- function(
  pairs,
  d,
  hap,
  ict,
  eta,
  eps,
  q,
  eigen_q,
  states,
  states_matrix,
  m,
  t_prob,
  total_time = 10,
  approx = 1,
  max_steps = 1000,
  mu = NULL,
  verbose = FALSE
) {
  n <- nrow(pairs) + 1

  if (!is.list(d)) {
    d <- list(d)
    hap <- list(hap)
  }
  k <- length(d)
  q <- as_barcode_list(q, k)
  eigen_q <- as_barcode_list(eigen_q, k)
  t_prob <- as_barcode_list(t_prob, k)

  # branch lengths are sums of intercoalescent times: b = A %*% ict
  a <- decompose_branches(pairs)
  curr_ict <- matrix(ict, nrow = n - 1, ncol = 1)
  curr_b <- a %*% curr_ict

  # the rate matrix restricted to the observed alleles only depends on the
  # data, so it is computed once
  q_allele_list <- list()
  p_list <- list()
  q_list <- list()
  curr_loglik <- 0
  for (i in seq_len(k)) {
    pq <- pq_vectors(
      d[[i]], hap[[i]], pairs, curr_b, eigen_q[[i]], states, states_matrix, m,
      approx, t_prob[[i]],
      mu = if (is.null(mu)) NULL else mu[[i]]
    )
    q_allele_list[[i]] <- q_alleles(pq[[1]], q[[i]], states, m, t_prob[[i]])
    p_list[[i]] <- pq[[2]]
    q_list[[i]] <- pq[[3]]
    curr_loglik <- curr_loglik + log(pq[[2]][1, ncol(pq[[2]])])
  }
  logliks <- curr_loglik

  diff <- 1
  num_steps <- 1
  g_t <- rep(0, n - 1)

  while (diff > eps) {
    grad_b <- 0
    for (i in seq_len(k)) {
      grad_b <- grad_b +
        gradient_branches(p_list[[i]], q_list[[i]], q_allele_list[[i]])
    }
    # (dl/d ict) = t(A) %*% (dl/d b), since b = A %*% ict
    grad_ict <- t(a) %*% matrix(grad_b, nrow = 2 * n - 2, ncol = 1)

    # adagrad step
    g_t <- g_t + grad_ict * grad_ict
    rho <- if (num_steps == 1) eta else eta / sqrt(g_t)
    old_ict <- curr_ict
    curr_ict <- curr_ict + rho * grad_ict

    # keep the times inside the constraint set sum(ict) <= total_time
    if (sum(curr_ict) > total_time && sum(curr_ict > 0) == (n - 1)) {
      curr_ict <- proj_l1(curr_ict, total_time)
    }

    # a step that produces a non-positive time is retried with a smaller rate
    while (sum(curr_ict <= 0) > 0) {
      eta <- eta / 5
      curr_ict <- old_ict
      rho <- if (num_steps == 1) eta else eta / sqrt(g_t)
      curr_ict <- curr_ict + rho * grad_ict
      if (sum(curr_ict) > total_time && sum(curr_ict > 0) == (n - 1)) {
        curr_ict <- proj_l1(curr_ict, total_time)
      }
    }

    curr_b <- a %*% curr_ict
    fit <- pq_all(
      pairs, d, hap, curr_b, eigen_q, states, states_matrix, m, t_prob,
      approx, mu
    )
    diff <- abs(logliks[length(logliks)] - fit$loglik)
    real_diff <- fit$loglik - logliks[length(logliks)]

    if (real_diff >= 0) {
      logliks <- c(logliks, fit$loglik)
    }

    # a step that decreases the log-likelihood is retried with a smaller rate
    while (real_diff < 0) {
      eta <- eta / 5
      curr_ict <- old_ict
      rho <- if (num_steps == 1) eta else eta / sqrt(g_t)
      curr_ict <- curr_ict + rho * grad_ict
      if (sum(curr_ict) > total_time && sum(curr_ict > 0) == (n - 1)) {
        curr_ict <- proj_l1(curr_ict, total_time)
      }
      curr_b <- a %*% curr_ict
      fit <- pq_all(
        pairs, d, hap, curr_b, eigen_q, states, states_matrix, m, t_prob,
        approx, mu
      )
      diff <- abs(logliks[length(logliks)] - fit$loglik)
      real_diff <- fit$loglik - logliks[length(logliks)]
    }

    p_list <- fit$p
    q_list <- fit$q

    if (verbose && num_steps %% 10 == 0) {
      print(c(curr_ict))
      print("Log-likelihood")
      print(logliks[length(logliks)])
      print(diff)
      print(num_steps)
    }

    num_steps <- num_steps + 1
    if (num_steps > max_steps) {
      break
    }
  }

  curr_ict <- as.vector(curr_ict)
  list(
    ict = curr_ict,
    coal_times = cumsum(curr_ict),
    logliks = logliks,
    steps = num_steps
  )
}
