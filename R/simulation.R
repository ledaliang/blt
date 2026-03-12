#' State Space Hash Table
#'
#' This generates the state space hash table for a given number of sites
#'
#' @param s Number of sites
#'
#' @return A hash table mapping state strings to their indices
#' @examples
#' state_space(3)
#' # <hash> containing 13 key-value pair(s).
#' # 000 : 13
#' # 001 : 5
#' # 010 : 10
#' # 011 : 2
#' # 033 : 7
#' # 100 : 12
#' # 101 : 4
#' # 110 : 9
#' # 111 : 1
#' # 133 : 6
#' # 220 : 11
#' # 221 : 3
#' # 222 : 8
#'
#' @export
state_space <- function(s) {
  states <- state_space_list(s)
  states <- relabel(states) # matrix
  h <- hash::hash()
  for (i in seq_len(nrow(states))) {
    state <- paste(states[i, ], collapse = "")
    h[[state]] <- i
  }
  h
}

#' State Space List
#'
#' The state space contains all possible vectors of length S,
#' where 0 indicates no mutation, 1 indicates mutation just
#' at that site, and larger numbers i are mutations that are
#' shared between multiple adjacent sites
#'
#' @param s Number of sites
#'
#' @return A list of all possible states
#'
#' @examples
#' state_space_list(3)
#' # [[111],[011],[221],[101],[001],[122],[022],[222],
#' # [110],[010],[220],[100],[000]]
#' @export
state_space_list <- function(s) {
  # recursive function
  if (s == 1) {
    list(c(1), c(0))
  } else {
    ss_1 <- state_space_list(s - 1)

    new <- list()
    for (state in ss_1) {
      new[[length(new) + 1]] <- c(1, state)
      new[[length(new) + 1]] <- c(0, state)

      if (state[1] == 1) { # can merge
        m <- max(state) + 1
        new[[length(new) + 1]] <- c(m, m, state[-1])
      }
      if (state[1] > 1) { # 1 is already part of an overlap, merge with that
        new[[length(new) + 1]] <- c(state[1], state)
      }
    }
    new
  }
}

#' Relabel States
#'
#' Helper function for state_space and state_space_matrix.
#' Renames overlapping mutations to (start + 1) where start
#' is the leftmost index. For example, state_space_list(3)
#' will return a list containing states 122 and 022 that will be
#' relabeled to 133 and 033 respectively.
#'
#' @param states A list of states
#'
#' @return A matrix of relabeled states
#'
#' @examples
#' states <- list(
#'   c(1, 1, 1),
#'   c(0, 1, 1),
#'   c(2, 2, 1),
#'   c(1, 0, 1),
#'   c(0, 0, 1),
#'   c(1, 2, 2),
#'   c(0, 2, 2),
#'   c(2, 2, 2),
#'   c(1, 1, 0),
#'   c(0, 1, 0),
#'   c(2, 2, 0),
#'   c(1, 0, 0),
#'   c(0, 0, 0)
#' )
#'
#' relabel(states)
#' #       [,1] [,2] [,3]
#' #  [1,]    1    1    1
#' #  [2,]    0    1    1
#' #  [3,]    2    2    1
#' #  [4,]    1    0    1
#' #  [5,]    0    0    1
#' #  [6,]    1    3    3
#' #  [7,]    0    3    3
#' #  [8,]    2    2    2
#' #  [9,]    1    1    0
#' # [10,]    0    1    0
#' # [11,]    2    2    0
#' # [12,]    1    0    0
#' # [13,]    0    0    0
#'
#' @export
relabel <- function(states) {
  l <- length(states) # states is a list
  new <- matrix(0, nrow = l, ncol = length(states[[1]]))

  for (i in seq_len(l)) {
    s <- states[[i]]
    new_s <- s
    if (max(s) > 1) {
      for (k in 2:max(s)) {
        overlap <- which(s == k)
        first <- overlap[1]
        new_s[overlap] <- first + 1
      }
    }
    new[i, ] <- new_s
  }
  new
}

#' State Space Matrix
#'
#' Turns the state space list into a matrix and relabels
#'
#' @param s Number of sites
#'
#' @returns An nxs matrix where each row represents a state
#' with one column for each site.
#'
#' @examples
#'
#' state_space_matrix(3)
#' #       [,1] [,2] [,3]
#' #  [1,]    1    1    1
#' #  [2,]    0    1    1
#' #  [3,]    2    2    1
#' #  [4,]    1    0    1
#' #  [5,]    0    0    1
#' #  [6,]    1    3    3
#' #  [7,]    0    3    3
#' #  [8,]    2    2    2
#' #  [9,]    1    1    0
#' # [10,]    0    1    0
#' # [11,]    2    2    0
#' # [12,]    1    0    0
#' # [13,]    0    0    0
#'
#' @export
state_space_matrix <- function(s) {
  states <- state_space_list(s)
  relabel(states)
}


#' Q Matrix
#'
#' Q is the generator matrix for a continuous time Markov chain over the list
#' of states. Theta is an SxS matrix with theta[i, j] is the rate of a cut from
#' site i to site j. Theta[i, i] is the rate of a single cut at site i.
#'
#'
#' @param states A matrix of states (output of state_space_matrix)
#'                with length n
#' @param theta An SxS matrix of cut rates
#'
#' @returns A nxn generator matrix
#'
#' @examples
#' states <- state_space_matrix(3)
#' #       [,1] [,2] [,3]
#' #  [1,]    1    1    1
#' #  [2,]    0    1    1
#' #  [3,]    2    2    1
#' #  [4,]    1    0    1
#' #  [5,]    0    0    1
#' #  [6,]    1    3    3
#' #  [7,]    0    3    3
#' #  [8,]    2    2    2
#' #  [9,]    1    1    0
#' # [10,]    0    1    0
#' # [11,]    2    2    0
#' # [12,]    1    0    0
#' # [13,]    0    0    0
#'
#' theta <- matrix(1, nrow=3, ncol=3)
#' #      [,1] [,2] [,3]
#' # [1,]    1    1    1
#' # [2,]    1    1    1
#' # [3,]    1    1    1
#'
#' q_matrix(states, theta)
#' #       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#' #  [1,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#' #  [2,]    1   -1    0    0    0    0    0    0    0     0     0     0     0
#' #  [3,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#' #  [4,]    1    0    0   -1    0    0    0    0    0     0     0     0     0
#' #  [5,]    0    1    1    1   -3    0    0    0    0     0     0     0     0
#' #  [6,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#' #  [7,]    0    0    0    0    0    1   -1    0    0     0     0     0     0
#' #  [8,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#' #  [9,]    1    0    0    0    0    0    0    0   -1     0     0     0     0
#' # [10,]    0    1    0    0    0    0    0    1    1    -3     0     0     0
#' # [11,]    0    0    1    0    0    0    0    0    0     0    -1     0     0
#' # [12,]    0    0    0    1    0    1    0    0    1     0     0    -3     0
#' # [13,]    0    0    0    0    1    0    1    1    0     1     1     1    -6
#'
#' @export
q_matrix <- function(states, theta) {
  l <- nrow(states)
  s <- ncol(states)
  m <- matrix(0, nrow = l, ncol = l)

  for (i in 1:l) {
    s1 <- states[i, ]
    active <- which(s1 == 0)
    if (length(active) > 0) {
      # First, all possible single mutations
      for (k in seq_along(active)) {
        site <- active[k]
        s2 <- s1
        s2[site] <- 1
        # This tells me which index the new state corresponds to,
        # so I can update the appropriate entry in the rate matrix
        # (Like a cheat version of a dictionary)
        j <- which(apply(
          states, 1,
          function(x) all.equal(x[1:s], s2)
        ) == "TRUE")
        m[i, j] <- theta[site, site]

        if (length(active) - k > 0) {
          # now all pairs of mutations
          for (l in (k + 1):length(active)) {
            site2 <- active[l]
            s2 <- s1
            # This mutation assignment is consistent with what we did before
            s2[site:site2] <- site + 1 # we should always have site2 > site
            j <- which(apply(
              states, 1,
              function(x) all.equal(x[1:s], s2)
            ) == "TRUE")
            m[i, j] <- theta[site, site2]
          }
        }
      }
    }
  }
  # Add negatives along diagonal so rows sum to 1
  rows <- rowSums(m)
  diagonal <- diag(rows, nrow = nrow(m), ncol = ncol(m))
  m <- m - diagonal
  m
}
#' Simulate Data Finite Alleles
#'
#' This function simulates allele data for all nodes of the tree
#' (the leaves are the first n rows), by starting with the state (0, 0, ...0)
#' at the root and simulating the Markov mutation process down the tree.
#'
#' @param tree An ape tree object
#' @param q Rate matrix for mutation states (output of q_matrix)
#' @param states_matrix Matrix storing possible mutation states
#' (output of state_space_matrix). M(i, i) is the number of alleles at site i,
#' M(i, j) is the number of alleles from a simultaneous cut at sites i, j.
#' In principle these numbers could be the same and we could input a single
#' parameter, but we will leave it as matrix for full generality.
#' @param m Symmetric matrix which stores the number of alleles per site or
#' site group
#'
#' @returns A list where the first entry is the matrix of mutation states
#' for the nodes and second entry is matrix of allele states.
#'
#' @examples
#' simul1 <- phylodyn::coalsim(samp_times = 0, n_sampled = 5, traj = exp_traj,method="tt",val_upper=11)
#' mytree <- phylodyn::generate_newick((simul1))$newick
#' # Phylogenetic tree with 5 tips and 4 internal nodes.
#' #
#' # Tip labels:
#' #   t2_0, t3_0, t4_0, t5_0, t1_0
#' #
#' #  Rooted; includes branch lengths.
#'
#' theta <- matrix(c(1, 1, 1, 0, 1, 1, 0, 0, 1), nrow = 3, byrow = TRUE)
#' S <- 3
#' states_matrix <- state_space_matrix(S)
#' Q <- q_matrix(states_matrix, theta)
#' #       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#' #  [1,]    0    0  0.0    0  0.0  0.0  0.0  0.0    0   0.0   0.0   0.0   0.0
#' #  [2,]    1   -1  0.0    0  0.0  0.0  0.0  0.0    0   0.0   0.0   0.0   0.0
#' #  [3,]    0    0  0.0    0  0.0  0.0  0.0  0.0    0   0.0   0.0   0.0   0.0
#' #  [4,]    1    0  0.0   -1  0.0  0.0  0.0  0.0    0   0.0   0.0   0.0   0.0
#' #  [5,]    0    1  0.5    1 -2.5  0.0  0.0  0.0    0   0.0   0.0   0.0   0.0
#' #  [6,]    0    0  0.0    0  0.0  0.0  0.0  0.0    0   0.0   0.0   0.0   0.0
#' #  [7,]    0    0  0.0    0  0.0  1.0 -1.0  0.0    0   0.0   0.0   0.0   0.0
#' #  [8,]    0    0  0.0    0  0.0  0.0  0.0  0.0    0   0.0   0.0   0.0   0.0
#' #  [9,]    1    0  0.0    0  0.0  0.0  0.0  0.0   -1   0.0   0.0   0.0   0.0
#' # [10,]    0    1  0.0    0  0.0  0.0  0.0  0.2    1  -2.2   0.0   0.0   0.0
#' # [11,]    0    0  1.0    0  0.0  0.0  0.0  0.0    0   0.0  -1.0   0.0   0.0
#' # [12,]    0    0  0.0    1  0.0  0.5  0.0  0.0    1   0.0   0.0  -2.5   0.0
#' # [13,]    0    0  0.0    0  1.0  0.0  0.5  0.2    0   1.0   0.5   1.0  -4.2
#'
#' #       [,1] [,2] [,3]
#' #  [1,]    1    1    1
#' #  [2,]    0    1    1
#' #  [3,]    2    2    1
#' #  [4,]    1    0    1
#' #  [5,]    0    0    1
#' #  [6,]    1    3    3
#' #  [7,]    0    3    3
#' #  [8,]    2    2    2
#' #  [9,]    1    1    0
#' # [10,]    0    1    0
#' # [11,]    2    2    0
#' # [12,]    1    0    0
#' # [13,]    0    0    0
#'
#' M <- matrix(20, nrow = S, ncol = S)
#' simulate_data_finite_alleles(mytree, Q, states_matrix, M)
#' # [[1]]
#' #      [,1] [,2] [,3]
#' # [1,]    2    2    1
#' # [2,]    0    0    1
#' # [3,]    1    1    0
#' # [4,]    1    1    0
#' # [5,]    0    0    1
#' # [6,]    0    0    0
#' # [7,]    0    0    0
#' # [8,]    0    0    0
#' # [9,]    1    1    0
#' #
#' # [[2]]
#' #      [,1]    [,2]    [,3]
#' # [1,] "12:17" "12:17" "3:14"
#' # [2,] "0"     "0"     "3:16"
#' # [3,] "1:5"   "2:5"   "0"
#' # [4,] "1:5"   "2:5"   "0"
#' # [5,] "0"     "0"     "3:13"
#' # [6,] "0"     "0"     "0"
#' # [7,] "0"     "0"     "0"
#' # [8,] "0"     "0"     "0"
#' # [9,] "1:5"   "2:5"   "0"
#'
#' @export
simulate_data_finite_alleles <- function(tree, q, states_matrix, m) {
  n <- tree$Nnode + 1
  s <- ncol(states_matrix)
  alleles_lumped <- matrix(0, 2 * n - 1, s) # To store alleles of all nodes
  alleles_full <- matrix(0, 2 * n - 1, s)

  # The leaves are labeled 1, ..., n and the interior nodes are
  # labeled n+1, ..2n - 1, with n+1 as the root
  edges <- tree$edge # Returns a list of (parent, child)
  edge_lengths <- tree$edge.length # 2n - 2 edges
  # loop through all edges
  for (i in 1:(2 * n - 2)) {
    parent <- edges[i, 1]
    child <- edges[i, 2]
    t <- edge_lengths[i]

    # Assuming that we're going down in the tree,
    # the parent allele should already be set
    parent_lumped <- alleles_lumped[parent, ]
    child_lumped <- transition2(parent_lumped, t, q, states_matrix)
    alleles_lumped[child, ] <- child_lumped

    # Now fill in the specific mutations in a consistent way
    parent_full <- alleles_full[parent, ]
    child_full <- rep(0, s)
    j <- 1
    while (j < s + 1) {
      if (child_lumped[j] > 0) {
        if (child_lumped[j] == parent_lumped[j]) {
          child_full[j] <- parent_full[j]
          j <- j + 1
        } else { # A new  mutation
          if (child_lumped[j] == 1) { # a new mutation at a single site
            # sample uniformly from the number of possible mutations
            mut <- sample(1:m[j, j], 1)

            # Just want a way to indicate this is a single mutation at site j
            # One way of doing that is to store as a string j:m
            child_full[j] <- paste(paste(j, ":", sep = ""), mut, sep = "")
            j <- j + 1
          } else { # an overlapping mutation
            # To assign all sites in the overlap the same mutation
            all <- which(child_lumped == child_lumped[j])
            # assume the mutation arose fro a simultaneous cut at sites i and j
            i <- max(all)
            mut <- sample(1:m[i, j], 1)

            # The overlap mutation will be stored ji:m
            # (where j is start site, i is end site)
            child_full[all] <- paste(paste(j * 10 + i, ":", sep = ""), mut,
              sep = ""
            )
            j <- j + length(all)
          }
        }
      } else {
        j <- j + 1
      }
    }
    alleles_full[child, ] <- child_full
  }
  list(alleles_lumped, alleles_full)
}


#' Transition Function
#'
#' Computes the transition probabilities given a parent state, time t,
#' and rate matrix q, and samples a new state.
#'
#' @param parent Parent state vector
#' @param t Time duration
#' @param q Rate matrix (output from q_matrix)
#' @param states Matrix of possible states (output from state_space_matrix)
#'
#' @returns A transition in time t from a given parent and rate matrix q
transition2 <- function(parent, t, q, states) {
  s <- ncol(states)
  i <- which(apply(states, 1, function(x) all.equal(x[1:s], parent)) == "TRUE")
  t_prob <- expm::expm(t * q)
  j <- sample(seq_len(nrow(states)), 1, replace = FALSE, prob = t_prob[i, ])

  states[j, ]
}


#' Lumped State
#'
#' Removes mutation labels from the full state to give a lumped state
#'
#' @param state sx1 vector of allele states
#'
#' @returns sx1 vector of lumped states using mutation state notation
#'
#' @examples
#' simul1 <- phylodyn::coalsim(samp_times = 0, n_sampled = 5, traj = exp_traj,method="tt",val_upper=11)
#' mytree <- phylodyn::generate_newick((simul1))$newick
#' theta <- matrix(c(1, 1, 1, 0, 1, 1, 0, 0, 1), nrow = 3, byrow = TRUE)
#' states_matrix <- state_space_matrix(3)
#' q <- q_matrix(states_matrix, theta)
#' M <- matrix(20, nrow = 3, ncol = 3)
#' datalist <- simulate_data_finite_alleles(mytree, q, states_matrix, M)
#' datalist[[2]]
#' #      [,1]   [,2]   [,3]
#' # [1,] "1:18" "0"    "3:7"
#' # [2,] "0"    "0"    "0"
#' # [3,] "0"    "2:8"  "0"
#' # [4,] "0"    "0"    "0"
#' # [5,] "1:20" "23:8" "23:8"
#' # [6,] "0"    "0"    "0"
#' # [7,] "0"    "0"    "0"
#' # [8,] "0"    "0"    "0"
#' # [9,] "0"    "0"    "0"
#'
#' lumped_state(datalist[[2]][1, ])
#' # 1 0 1
#' lumped_state(datalist[[2]][5, ])
#' # 1 3 3
#'
#' @export
lumped_state <- function(state) {
  s <- length(state)
  lumped <- rep(0, s)

  i <- 1

  if (length(state) == 0) {
    print("lumped state error")
    print(state)
  }
  while (i < s + 1) {
    if (state[i] == "None") {
      state[i] <- 0
    }
    if (state[i] == "0") {
      state[i] <- 0
    }
    if (state[i] != 0) {
      if (i < s) {
        if (state[i] != state[i + 1]) {
          lumped[i] <- 1
          i <- i + 1
        } else { # overlapping to mark
          all <- which(state == state[i])
          lumped[all] <- i + 1
          i <- i + length(all)
        }
      } else {
        lumped[i] <- 1
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  lumped
}


#' Observed Mutations
#'
#' This function counts the number of observed mutations from the matrix D,
#' which stores alleles data from the tips. Note that this function will return
#' counts only based on the lumped mutation state, not the allele states
#'
#' @param d Matrix of allele states at the tips
#' @param s Number of sites
#' @param lumped_d (Optional) Matrix of lumped states at the tips
#'
#' @returns An sxs matrix where s[i,i] is number of cells with a single mutation
#'          in site i, s[i, j], i < j, is number of cells with an overlapping
#'          mutation going from i to j
#'
#' @examples
#' simul1 <- phylodyn::coalsim(samp_times = 0, n_sampled = 5, traj = exp_traj,method="tt",val_upper=11)
#' mytree <- phylodyn::generate_newick((simul1))$newick
#' theta <- matrix(c(1, 1, 1, 0, 1, 1, 0, 0, 1), nrow = 3, byrow = TRUE)
#' states_matrix <- state_space_matrix(3)
#' q <- q_matrix(states_matrix, theta)
#' M <- matrix(20, nrow = 3, ncol = 3)
#' datalist <- simulate_data_finite_alleles(mytree, q, states_matrix, M)
#' observed_mutations(datalist[[2]], 3)
#' # $obsM
#' #      [,1] [,2] [,3]
#' # [1,]    2    0    0
#' # [2,]    0    1    1
#' # [3,]    0    0    1
#' #
#' # $obsA
#' # [1] 2 1 1 0 1 0
#' @export
observed_mutations <- function(d, s, lumped_d = NA) {
  obs_m <- matrix(0, nrow = s, ncol = s)
  obs_a <- matrix(0, nrow = s, ncol = s)
  n <- nrow(d)

  # check if lumped states provided, otherwise have to get that now
  if (max(is.na(lumped_d)) > 0) {
    lumped_d <- matrix(0, nrow = n, ncol = s)
    for (r in 1:n) {
      lumped_d[r, ] <- lumped_state(d[r, ])
    }
  }

  for (i in 1:s) {
    have_mut <- which(lumped_d[, i] == 1)
    obs_m[i, i] <- length(have_mut)
    obs_a[i, i] <- length(unique(d[have_mut, i]))
  }
  for (i in 1:(s - 1)) {
    for (j in (i + 1):s) {
      have_mut <- which(lumped_d[, i] == (i + 1) & lumped_d[, j] == (i + 1))
      if (j < s) {
        have_mut <- intersect(have_mut, which(lumped_d[, j + 1] != i + 1))
        # this identifies the rows which have an overlap exactly between i, j
      }
      obs_m[i, j] <- length(have_mut)
      obs_a[i, j] <- length(unique(d[have_mut, i]))
    }
  }
  vec_a <- 0
  for (j in 1:s) {
    vec_a <- c(vec_a, obs_a[cbind(1:(s - j + 1), j:s)])
  }
  list(obsM = obs_m, obsA = vec_a[-1])
}
#' Probability of Mutations
#'
#' This function computes the probability of observing a mutation of a given
#' type given the rate vector theta and time t.
#'
#' @param theta (s(s+1)/2)x1 vector, where first s entries are single cut rates
#'               (the diagonal of theta matrix), the next s-1 entries are the
#'               2 site mutation rates (off diagonal), and so on.
#' @param observed_muts Observed mutations matrix
#'                      (output of observed_mutations2$obsM))
#' @param t Time duration
#' @param s Number of sites
#' @param states List of states (output of state_space_list)
#' @param states_matrix Matrix of states (output of state_space_matrix)
#'
#' @returns An sxs matrix with (i, i) the probability of a single cut at site i,
#'          and (i, j), i < j, the probability of an overlapping cut from
#'          sites i to j.
#' @export
probs_function <- function(theta, observed_muts, t, s, states, states_matrix) {
  theta_matrix <- matrix(0, nrow = s, ncol = s)
  pos <- c(0, 0)
  for (j in 1:s) {
    pos <- rbind(pos, cbind(1:(s - j + 1), j:s))
  }
  theta_matrix[pos] <- theta
  q <- q_matrix(states_matrix, theta_matrix)
  eigen_q <- eigen(q)
  eigen_q$inv <- solve(eigen_q$vectors)
  p_mat <- eigen_q$vectors %*% diag(exp(eigen_q$values * t)) %*% eigen_q$inv

  probs <- matrix(0, nrow = s, ncol = s)

  # Should always be the last entry in the states list, but just to check
  zero_state_index <- states[[paste(rep(0, s), collapse = "")]]
  for (i in 1:s) {
    # Find the indices of the states which have 1 in ith position
    indices <- which(states_matrix[, i] == 1)
    probs[i, i] <- sum(p_mat[zero_state_index, indices])
  }

  for (i in 1:(s - 1)) {
    for (j in (i + 1):s) {
      # Based on how we assign mutation states, these are the states with an
      # overlap going from i to j (exactly)
      indices <- which(
        states_matrix[, i] == (i + 1) &
          states_matrix[, j] == (i + 1)
      )
      probs[i, j] <- sum(p_mat[zero_state_index, indices])
    }
  }
  diff <- probs - observed_muts
  sum(diff^2)
}


#' Tree mapping representation
#'
#' This function returns the matching representation of an ape tree as a matrix.
#' Each row is a merge and the third entry in the row is the label of the
#' interior node for th emerge. Leaves are 1,...,n and interior nodes are
#' n+1, .., 2n-1 (Note that ape labels with te root as n+1)
#' Warning: this may give invalid tree if an edge length is 0
#' (which can happen with upgma)
#'
#' @param tree An ape tree object
#'
#' @returns (n-1)x3 matrix where each row is a merge event. The first two
#' coordinates represent the labels of the children and the third coordinate is
#' the label of the parent node
#' @examples
#' simul1 <- phylodyn::coalsim(samp_times = 0, n_sampled = 5, traj = exp_traj,method="tt",val_upper=11)
#' mytree <- phylodyn::generate_newick((simul1))$newick
#' matching(mytree)
#' #      [,1] [,2] [,3]
#' # [1,]    1    8    9
#' # [2,]    2    7    8
#' # [3,]    6    5    7
#' # [4,]    3    4    6
#'
#' @export
matching <- function(tree) {
  n <- tree$Nnode + 1

  # if any edges are 0, then interior nodes might not be ordered which would
  # give an invalid tree
  # a quick fix, just add tiny number to any 0s
  tree$edge.length[which(tree$edge.length == 0)] <- 0.0000000001

  t_tot <- max(ape::node.depth.edgelength(tree)) # the total tree length
  # distance from leaves to each node
  n_t <- t_tot - ape::node.depth.edgelength(tree)


  # Sorts according to distance from leaves ... You would hope that the nodes
  # are already labeled in logical order according to their coal time,
  # e.g. interior nodes are ordered (2n-1), (2n-2), ..., (n+1)
  # but apparently this doesn't always happen so we need to check
  xx <- sort(n_t, index.return = TRUE)
  newlabel <- seq(n + 1, 2 * n - 1)
  oldlabel <- xx$ix[xx$ix > n]

  for (j in seq_len(nrow(tree$edge))) {
    # the first entry in each edge is always an interior node
    tree$edge[j, 1] <- newlabel[oldlabel == tree$edge[j, 1]]
    if (tree$edge[j, 2] > n) { # if entry is a leaf, don't need to change
      tree$edge[j, 2] <- newlabel[oldlabel == tree$edge[j, 2]]
    }
  }

  # the inner nodes are already relabeled, now form the matches
  match <- matrix(0, nrow = n - 1, ncol = 3)
  for (r in 1:(n - 1)) {
    # find which children form 2*n - r merge,
    child <- tree$edge[tree$edge[, 1] == (2 * n - r), 2]
    # this will be a vector with two entries
    match[r, ] <- c(child, 2 * n - r)
  }
  match
}


#' Decompose Branches
#'
#' This function decomposes the branch lengths of the tree into a matrix of
#' intercalescent times.
#'
#' @param pairs Matching representation of the tree (output of matching
#' function)
#'
#' @returns A (2n -2)x(n-1) matrix where each row represents a branch and
#'          each column represents an intercoalescent time. Each branch length
#' b_j can be written as a sum of intercoalescent times. A[j, i] = 1 if branch j
#'  uses intercoalescent time i.
#'
#' @examples
#' simul1 <- phylodyn::coalsim(samp_times = 0, n_sampled = 5, traj = exp_traj,method="tt",val_upper=11)
#' mytree <- phylodyn::generate_newick((simul1))$newick
#' plot(mytree)
#' ape::nodelabels()
#' pairs <- matching(mytree)
#' #      [,1] [,2] [,3]
#' # [1,]    1    8    9
#' # [2,]    2    7    8
#' # [3,]    6    5    7
#' # [4,]    3    4    6
#'
#' decompose_branches(pairs)
#' #     [,1] [,2] [,3] [,4]
#' # [1,]    1    1    1    1
#' # [2,]    1    1    1    0
#' # [3,]    1    0    0    0
#' # [4,]    1    0    0    0
#' # [5,]    1    1    0    0
#' # [6,]    0    1    0    0
#' # [7,]    0    0    1    0
#' # [8,]    0    0    0    1
#'
#' @export
decompose_branches <- function(pairs) {
  n <- nrow(pairs) + 1

  a <- matrix(0, nrow = n - 1, ncol = 2 * n - 2)

  for (r in 1:(n - 1)) { # loop through pairs
    parent <- pairs[r, 3]

    child1 <- pairs[r, 1]
    child2 <- pairs[r, 2]

    if (child1 <= n) {
      a[1:(parent - n), child1] <- 1
    } else {
      a[(child1 - n + 1):(parent - n), child1] <- 1
    }

    if (child2 <= n) {
      a[1:(parent - n), child2] <- 1
    } else {
      a[(child2 - n + 1):(parent - n), child2] <- 1
    }
  }
  t(a)
}

#' Compute p and q Vectors for Gradient Calculations under a
#' Finite Alleles Model
#'
#' @description
#' Computes the forward (\eqn{p}) and backward (\eqn{q}) probability vectors
#' for each node in a phylogenetic tree under a finite-alleles mutation model.
#' The \eqn{p}-vectors represent the conditional probability of observing all
#' descendant data given a specific allele at a node, while the \eqn{q}-vectors
#' represent the probability of observing all non-descendant data given a
#' specific allele at that node. These vectors can be used in gradient-based
#' optimization of likelihoods.
#'
#' @param d A numeric matrix of dimension \eqn{n \times S}, where \eqn{n} is the
#'  number of haplotypes (tips) and \eqn{S} is the number of sites. Each row
#' corresponds to observed alleles for a haplotype.
#' @param hap A list containing:
#'   \describe{
#'     \item{\code{p}}{Initial \eqn{p}-matrix with leaf initialization.}
#'     \item{\code{alleles}}{Matrix of observed or possible alleles, each row
#' representing an allele vector.}
#'   }
#' @param pairs A matrix where each row corresponds to a pair of child node
#' indices and their parent node index: \code{c(child1, child2, parent)}.
#' @param edge_lengths A numeric vector of edge lengths for the tree (in the
#' same node ordering as `pairs`).
#' @param eigen_q A list containing the eigendecomposition of the mutation rate
#' matrix \eqn{Q}: \code{vectors}, \code{values}, and \code{inv}.
#' @param states A named list of possible allele states indexed by string
#' representations.
#' @param states_matrix A matrix representation of the possible allele states
#' for faster lookup.
#' @param m An integer matrix defining the number of alleles allowed at each
#' site and overlap: \code{m[i, i]} gives the number of alleles allowed at site
#' \eqn{i} and \code{m[i, j]} gives the number of alleles for an overlap from
#' site \eqn{i} to \eqn{j}.
#' @param approx Integer flag. Set to `0` for exact likelihood mode (all allele
#' combinations considered), or `1` for approximate likelihood mode (only
#' alleles present in the leaves or wildcard).
#' @param t_prob Transition probability matrix (if precomputed), otherwise
#' constructed from `eigen_q`.
#' @param mu Optional list structure containing conditional mutation probability
#' distributions for each overlap. If `NULL`, a uniform distribution over
#' allowed alleles in `m` is assumed.
#'
#' @return A list with three components:
#' \describe{
#'   \item{\code{alleles}}{A matrix of all unique alleles encountered, one per
#' row.}
#'   \item{\code{p}}{An \eqn{A \times (2n - 1)} matrix of forward
#' probabilities.}
#'   \item{\code{q}}{An \eqn{A \times (2n - 1)} matrix of backward
#' probabilities.}
#' }
pq_vectors <- function(
  d,
  hap,
  pairs,
  edge_lengths,
  eigen_q,
  states,
  states_matrix,
  m,
  approx,
  t_prob,
  mu = NULL
) {
  n <- nrow(d)
  s <- ncol(d)

  p <- hap$p
  alleles <- hap$alleles

  # Now go up the tree using the pairs
  if (approx == 0) {
    for (j in 1:(n - 1)) {
      node <- n + j
      row <- n - j
      c1 <- pairs[row, 1]
      c2 <- pairs[row, 2]
      pt1 <- eigen_q$vectors %*%
        diag(exp(eigen_q$values * edge_lengths[c1])) %*%
        eigen_q$inv
      pt2 <- eigen_q$vectors %*%
        diag(exp(eigen_q$values * edge_lengths[c2])) %*%
        eigen_q$inv


      # There must be a better way to do this??
      pos_child1 <- which(p[, c1] > 0)
      pos_child2 <- which(p[, c2] > 0)

      for (i in seq_along(pos_child1)) {
        child1 <- alleles[pos_child1[i], ]
        c1_lumped <- lumped_state(child1)
        pc1 <- p[pos_child1[i], c1]
        parents1 <- possible_parents_fa(
          child1, states, states_matrix, t_prob, m
        )
        for (l in seq_along(pos_child2)) {
          child2 <- alleles[pos_child2[l], ]
          c2_lumped <- lumped_state(child2)
          pc2 <- p[pos_child2[l], c2]
          parents2 <- possible_parents_fa(
            child2, states, states_matrix, t_prob, m
          )

          possible <- Reduce(intersect, list(parents1, parents2))
          for (k in seq_along(possible)) {
            allele <- possible[[k]]
            allele_lumped <- lumped_state(allele)

            p1 <- transition_prob_finite_alleles(
              allele, child1, pt1, states, m, allele_lumped, c1_lumped,
              mu = mu
            )
            p2 <- transition_prob_finite_alleles(
              allele, child2, pt2, states, m, allele_lumped, c2_lumped,
              mu = mu
            )

            allele_index <- which(apply(
              alleles, 1,
              function(x) all.equal(x[1:s], allele)
            ) == "TRUE")
            # The allele has already been added to list
            if (length(allele_index) > 0) {
              p[allele_index, node] <- p[allele_index, node] +
                p1 * p2 * pc1 * pc2
            } else {
              alleles <- rbind(alleles, allele)
              p <- rbind(p, rep(0, ncol(p)))
              p[nrow(p), node] <- p1 * p2 * pc1 * pc2
            }
          }
        }
      }
    }
  } else {
    for (j in 1:(n - 1)) {
      node <- n + j

      # The pairs should be stored in decreasing order, but just to be safe...
      row <- which(pairs[, 3] == node)
      c1 <- pairs[row, 1]
      c2 <- pairs[row, 2]
      pt1 <- eigen_q$vectors %*%
        diag(exp(eigen_q$values * edge_lengths[c1])) %*%
        eigen_q$inv
      pt2 <- eigen_q$vectors %*%
        diag(exp(eigen_q$values * edge_lengths[c2])) %*%
        eigen_q$inv

      pos_child1 <- which(p[, c1] > 0)
      pos_child2 <- which(p[, c2] > 0)

      for (i in seq_along(pos_child1)) {
        child1 <- alleles[pos_child1[i], ]
        c1_lumped <- lumped_state(child1)
        c1_index <- states[[stringi::stri_join(c1_lumped, collapse = "")]]

        pc1 <- p[pos_child1[i], c1]

        parents1 <- possible_parents_approximate(
          child1, states, states_matrix, t_prob, d
        )
        for (l in seq_along(pos_child2)) {
          child2 <- alleles[pos_child2[l], ]
          c2_lumped <- lumped_state(child2)
          c2_index <- states[[stringi::stri_join(c2_lumped, collapse = "")]]

          pc2 <- p[pos_child2[l], c2]

          parents2 <- possible_parents_approximate(
            child2, states, states_matrix, t_prob, d
          )
          possible <- Reduce(intersect, list(parents1, parents2))

          for (k in seq_along(possible)) {
            allele <- possible[[k]]
            allele_lumped <- lumped_state(allele)
            parent_index <- states[[
              stringi::stri_join(allele_lumped, collapse = "")
            ]]

            p1 <- transition_prob_finite_alleles(
              allele, child1, pt1, states, m, allele_lumped, c1_lumped,
              parent_index, c1_index,
              mu = mu
            )
            p2 <- transition_prob_finite_alleles(allele, child2, pt2, states,
              m, allele_lumped, c2_lumped, parent_index, c2_index,
              mu = mu
            )

            allele_index <- which(apply(
              alleles, 1,
              function(x) all.equal(x[1:s], allele)
            ) == "TRUE")
            # The allele has already been added to list
            if (length(allele_index) > 0) {
              p[allele_index, node] <- p[allele_index, node] +
                p1 * p2 * pc1 * pc2
            } else {
              alleles <- rbind(alleles, allele)
              p <- rbind(p, rep(0, ncol(p)))
              p[nrow(p), node] <- p1 * p2 * pc1 * pc2
            }
          }
        }
      }
    }
  }

  # The q vectors are discovered from the root down
  q <- matrix(0, nrow = nrow(alleles), ncol = 2 * n - 1)
  # This condition indicates the root equals ancestral, which should be the
  # first entry in the alleles matrix anyway
  q[1, 2 * n - 1] <- 1

  for (k in 1:(2 * n - 2)) {
    node <- (2 * n - 1) - k
    pt1 <- eigen_q$vectors %*%
      diag(exp(eigen_q$values * edge_lengths[node])) %*%
      eigen_q$inv
    row <- c(which(pairs[, 1] == node), which(pairs[, 2] == node))
    parent <- pairs[row, 3]
    sibling <- pairs[row, c(which(pairs[row, 1:2] != node))]
    pt2 <- eigen_q$vectors %*%
      diag(exp(eigen_q$values * edge_lengths[sibling])) %*%
      eigen_q$inv

    for (a in seq_len(nrow(alleles))) {
      allele <- alleles[a, ]
      q_total <- 0
      for (b in seq_len(nrow(alleles))) {
        allele2 <- alleles[b, ]
        qb <- q[b, parent]
        p1 <- transition_prob_finite_alleles(
          allele2, allele, pt1, states, m,
          mu = mu
        )
        for (c in seq_len(nrow(alleles))) {
          allele3 <- alleles[c, ]
          p2 <- transition_prob_finite_alleles(
            allele2, allele3, pt2, states, m,
            mu = mu
          )
          pc <- p[c, sibling]
          q_total <- q_total + qb * p1 * p2 * pc
        }
      }
      q[a, node] <- q_total
    }
  }
  list(alleles, p, q)
}


#' Initialize Haplotypes and Leaf Probabilities for Finite Alleles Model
#'
#' @description
#' Constructs a summary of the observed haplotypes in the data and initializes
#' the corresponding forward probability matrix (\eqn{p}) for use in the
#' \code{pq_vectors()} function. Each unique haplotype is treated as a distinct
#' allele and is assigned a row in both the allele matrix and the
#' \eqn{p}-matrix.
#'
#' @details
#' The function initializes two objects:
#' \describe{
#' \item{\code{alleles}}{A matrix of dimension \eqn{A \times S}, where each
#' row corresponds to one of the \eqn{A} unique alleles observed in the data,
#' and each column corresponds to a genetic site (of length \eqn{S}).}
#' \item{\code{p}}{A matrix of dimension \eqn{A \times (2n - 1)} containing
#' the initial forward probabilities for each allele at each tree node.
#' The first \eqn{n} columns correspond to the observed leaf nodes and are
#' initialized as indicator vectors for the observed alleles.}
#' }
#'
#' @param d A character or numeric matrix of dimension \eqn{n \times S}
#' containing the observed alleles (haplotypes) for each of \eqn{n} leaves
#' across \eqn{S} sites. Each row corresponds to one haplotype.
#' @param n Integer. The number of haplotypes (number of rows in \code{d}).
#' @param s Integer. The number of sites (number of columns in \code{d}).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{alleles}}{A matrix of unique observed alleles, one per row.}
#'   \item{\code{p}}{An \eqn{A \times (2n - 1)} matrix of initialized forward
#' probabilities.}
#' }
haplotypes <- function(d, n, s) {
  # Because it is function of the data only, we only need to run it once
  # Initialize the p-vectors with the data
  alleles <- matrix("0", nrow = 1, ncol = s)
  p <- matrix(0, ncol = 2 * n - 1, nrow = 1)
  for (r in 1:n) {
    node <- d[r, ]
    allele_index <- which(
      apply(alleles, 1, function(x) all.equal(x[1:s], node)) == "TRUE"
    )
    if (length(allele_index) > 0) { # The allele has already been added to list
      p[allele_index, r] <- 1
    } else {
      alleles <- rbind(alleles, node)
      p <- rbind(p, rep(0, ncol(p)))
      p[nrow(p), r] <- 1
    }
  }
  list(alleles = alleles, p = p)
}

#' Possible Parental Alleles Under Finite Alleles Model
#'
#' @description
#' Given an observed allele state at a node (child), this function returns all
#' possible allele states of its parent under a finite alleles mutation model
#' with potential overlapping mutations. The function accounts for cases where
#' mutations at multiple sites overlap in state space and may mask earlier
#' mutations.
#'
#' @details
#' The output is a list where each element is a character vector of length
#' \eqn{S} representing a valid allele configuration in the parent that could
#' lead to the child allele through mutation and masking events. The number of
#' such configurations can grow exponentially in the case of overlapping or
#' masked mutations.
#'
#' Parent states are determined according to lumped allele categories (i.e.\
#' allele states grouped based on their mutation class). The mutation landscape
#' is defined by:
#' \itemize{
#'   \item \code{M[i, i]}: number of alleles allowed at site \eqn{i}
#'   \item \code{M[i, j]}: number of alleles allowed for an overlap spanning
#'   sites \eqn{i} through \eqn{j} (inclusive)
#' }
#'
#' If the child is the ancestral type (i.e., no derived mutations), then only
#' the ancestral allele is returned. For overlapping alleles, the function calls
#' \code{all_allele_states()} to enumerate fully expanded/masked possibilities.
#'
#' @param child A character vector of length \eqn{S} representing the allele
#' state of the child at \eqn{S} sites.
#' @param states A named list of valid lumped allele states indexed by their
#' string representations.
#' @param states_matrix A numeric matrix encoding allele state transitions
#' for efficient lookup.
#' @param t_prob A precomputed transition probability structure used for
#' validating compatible parent states. Should correspond to the same mutation
#' model used throughout.
#' @param m An integer matrix where \code{m[i, i]} gives the number of possible
#' alleles at site \eqn{i}, and \code{m[i, j]} (for \eqn{i < j}) gives the
#' number of alleles for a mutation overlap spanning from site \eqn{i} to
#' \eqn{j}.
#'
#' @return A list where each element is a character vector representing a valid
#' parental allele configuration that could give rise to the child allele
#' through mutation or overlapping events.
possible_parents_fa <- function(child, states, states_matrix, t_prob, m) {
  # This a matrix, rows are states
  lumped_parents <- possible_parents_lumped(
    lumped_state(child), states, states_matrix, t_prob
  )
  lumped_child <- lumped_state(child)
  s <- length(child)

  parents <- list()

  # If the child is ancestral type, only the ancestral type is possible
  if (sum(lumped_child) == 0) {
    # Note: important because alleles are characters, mutations are integers
    return(list(as.character(lumped_parents)))
  } else if (nrow(lumped_parents) > 0) {
    for (r in seq_len(nrow(lumped_parents))) {
      state <- lumped_parents[r, ]

      # This could be a matrix, if there is more than one
      # possible parent due to masking mutations
      parent <- state
      s <- 1
      while (s < s + 1) {
        if (lumped_child[s] == state[s]) {
          # already more possibilities from a previous overlap
          if (!is.null(nrow(parent))) {
            parent[, s] <- rep(child[s], nrow(parent))
          } else {
            parent[s] <- child[s]
          }
          news <- s + 1
        } else {
          if (lumped_child[s] == 1) {
            news <- s + 1
          } else {
            # child has a new overlap mutation
            # This case is if the parent has a mutation which vanished in the
            # child due to an overlapping mutation. We will generate all
            # possible mutations, using M
            endstate <- max(which(lumped_child == lumped_child[s]))
            news <- endstate + 1
            if (sum(state[s:endstate]) > 0) { # indicating an overlap

              temp_parent <- rep(0, s)
              temp_parent[s:endstate] <- state[s:endstate]
              possible_masked <- all_allele_states(temp_parent, m)

              # already more possibilities from a previous overlap
              if (!is.null(nrow(parent))) {
                results <- c()
                for (r in seq_len(nrow(parent))) {
                  for (k in seq_len(nrow(possible_masked))) {
                    new_parent <- parent[r, ]
                    new_parent[s:endstate] <- possible_masked[k, s:endstate]
                    results <- c(results, new_parent)
                  }
                }
                parent <- matrix(results,
                  nrow = nrow(possible_masked) * nrow(parent),
                  ncol = s, byrow = TRUE
                )
              } else {
                results <- c()
                for (r in seq_len(nrow(possible_masked))) {
                  new_parent <- parent
                  new_parent[s:endstate] <- possible_masked[r, s:endstate]
                  results <- c(results, new_parent)
                }
                parent <- matrix(results,
                  nrow = nrow(possible_masked), ncol = s,
                  byrow = TRUE
                )
              }
            }
          }
        }
        s <- news
      }

      if (!is.null(nrow(parent))) { # more than one option
        for (r in seq_len(nrow(parent))) {
          parents[[length(parents) + 1]] <- as.character(parent[r, ])
        }
      } else {
        parents[[length(parents) + 1]] <- as.character(parent)
      }
    }
  }
  parents
}


#' Compute Approximate Parental Allele States Using Leaf Alleles
#'
#' @description
#' Computes all valid parental allele configurations for a given child allele
#' under an approximate finite-alleles mutation model. Unlike the exact version
#' implemented in \code{\link{possible_parents_fa}}, this function restricts
#' overlap-resolving mutations to only those alleles that appear in the observed
#' leaf data (plus one wildcard for unseen alleles).
#'
#' @details
#' The function enumerates possible parent allele states for the child based on
#' the lumped allele transition model and partial masking from overlapping
#' mutations. For overlapping sites, rather than generating all possible
#' mutations from the mutation count matrix \code{M}, this approximate approach
#' limits possible alleles at each masked site to those observed in the leaf
#' set (plus a single wildcard to model any unseen allele).
#'
#' This can drastically reduce the number of possible masked configurations,
#' making computation practical for large trees, while still preserving the
#' major allele patterns.
#'
#' @param child A character vector of length \eqn{S} representing the observed
#' allele states of the child node across \eqn{S} sites.
#' @param states A named list of valid lumped allele states, indexed by string
#' encodings.
#' @param states_matrix A matrix encoding valid combinations of lumped allele
#' states, for efficient lookup and indexing during recombination.
#' @param t_prob A precomputed list or structure of transition probabilities,
#' governing lumped allele transitions.
#' @param leaves A matrix where each row is a valid allele configuration
#' observed in the tree leaves (tip nodes). These are used to restrict masked
#' overlaps.
#'
#' @return
#' A list of character vectors, each of length \eqn{S}, corresponding to
#' possible parental allele configurations for the given child under the
#' approximate scheme.
possible_parents_approximate <- function(
  child, states, states_matrix, t_prob, leaves
) {
  # This a matrix, rows are states
  lumped_parents <- possible_parents_lumped(
    lumped_state(child), states, states_matrix, t_prob
  )
  lumped_child <- lumped_state(child)
  s <- length(child)

  parents <- list()

  # If the child is ancestral type, only the ancestral type is possible
  if (sum(lumped_child) == 0) {
    # Note: important because alleles are characters, mutations are integers
    return(list(as.character(lumped_parents)))
  } else if (nrow(lumped_parents) > 0) {
    for (r in seq_len(nrow(lumped_parents))) {
      state <- lumped_parents[r, ]

      # This could be a matrix, if there is more than one
      # possible parent due to masking mutations
      parent <- state
      s <- 1
      while (s < s + 1) {
        if (lumped_child[s] == state[s]) {
          # already more possibilities from a previous overlap
          if (!is.null(nrow(parent))) {
            parent[, s] <- rep(child[s], nrow(parent))
          } else {
            parent[s] <- child[s]
          }
          news <- s + 1
        } else {
          if (lumped_child[s] == 1) {
            news <- s + 1
          } else { # child has a new overlap mutation
            # This case is if the parent has a mutation which vanished in the
            # child due to an overlapping mutation. We will generate all
            # possible configurations, but only from alleles which occurred in
            # the leaves, and so show up elsewhere in the tree. We also add on
            # one `wildcard' indicating something not seen in the tree
            # !!!!NOTE!!! There could be more than one WC introduced at the same
            # site. Currently the code treats all WC as the same allele,
            # ignoring the possibility (which seems more likely). That the WC
            # are different alleles ... could modify this later to treat all WC
            # as unique
            endstate <- max(which(lumped_child == lumped_child[s]))
            news <- endstate + 1
            if (sum(state[s:endstate]) > 0) { # indicating an overlap

              temp_parent <- rep(0, s)
              temp_parent[s:endstate] <- state[s:endstate]
              possible_masked <- all_allele_states_leaves(temp_parent, leaves)

              # already more possibilities from a previous overlap
              if (!is.null(nrow(parent))) {
                results <- c()
                for (r in seq_len(nrow(parent))) {
                  for (k in seq_len(nrow(possible_masked))) {
                    new_parent <- parent[r, ]
                    new_parent[s:endstate] <- possible_masked[k, s:endstate]
                    results <- c(results, new_parent)
                  }
                }
                parent <- matrix(results,
                  nrow = nrow(possible_masked) * nrow(parent),
                  ncol = s, byrow = TRUE
                )
              } else {
                results <- c()
                for (r in seq_len(nrow(possible_masked))) {
                  new_parent <- parent
                  new_parent[s:endstate] <- possible_masked[r, s:endstate]
                  results <- c(results, new_parent)
                }
                parent <- matrix(results,
                  nrow = nrow(possible_masked),
                  ncol = s, byrow = TRUE
                )
              }
            }
          }
        }
        s <- news
      }

      if (!is.null(nrow(parent))) { # more than one option
        for (r in seq_len(nrow(parent))) {
          parents[[length(parents) + 1]] <- as.character(parent[r, ])
        }
      } else {
        parents[[length(parents) + 1]] <- as.character(parent)
      }
    }
  }
  parents
}


#' Compute Lumped Parental Allele States
#'
#' @description
#' Given the lumped allele state of a child, this function returns all valid
#' lumped parental allele states with non-zero transition probability under a
#' specified mutation model. This is a preliminary step used by both exact and
#' approximate parental expansion functions to efficiently prune impossible
#' transitions.
#'
#' @details
#' The \emph{lumped} allele model groups allele states by mutation class (or
#' mutational load), reducing the full mutation model to a smaller state space.
#' A parental state is included if the corresponding entry in the transition
#' probability matrix \code{t_prob} (representing \eqn{P(child | parent}) is
#' strictly positive.
#'
#' The ancestral (zero-mutation) state is automatically included as a fallback,
#' even if not explicitly allowed by \code{t_prob}. This ensures numerical
#' stability in downstream expansion of masked or overlapping allele states.
#'
#' @param child A numeric or integer vector of length \eqn{S} representing the
#' lumped state of the child across \eqn{S} sites. Each element typically
#' encodes mutation class at that site.
#' @param states A named list mapping string representations of lumped states to
#' integer indices in \code{states_matrix}.
#' @param states_matrix A matrix of dimension \eqn{K \times S} where each row
#' corresponds to a valid lumped allele configuration.
#' @param t_prob A matrix of dimension \eqn{K \times K} giving the transition
#' probabilities between lumped states, where entry \code{t_prob[i, j]} encodes
#' \eqn{P(\text{child state } j \mid \text{parent state } i)}.
#'
#' @return
#' A matrix of dimension \eqn{L \times S}, where each of the \eqn{L} rows is a
#' lumped parental allele state that can give rise to the specified child state
#' with non-zero probability.
possible_parents_lumped <- function(child, states, states_matrix, t_prob) {
  j <- states[[paste(child, collapse = "")]]

  parents_indices <- which(t_prob[, j] > 0)
  parents_indices <- unique(c(parents_indices, nrow(states_matrix)))

  states_matrix[parents_indices, ]
}

#' Generate Site Allele Labels
#'
#' Creates a vector of allele labels of the form `"s:j"` for
#' indices from 1 through `num`.
#'
#' @param s A character or numeric identifier for the site.
#' @param num Integer giving the number of alleles to generate.
#'
#' @return A character vector of length `num` with elements
#' `"s:1"`, `"s:2"`, ..., `"s:num"`.
#'
#' @examples
#' site_alleles(3, 3)
#' # c("3:1", "3:2", "3:3")
#'
#' @export
site_alleles <- function(s, num) {
  paste0(s, ":", seq_len(num))
}


#' Enumerate All Possible Allele States for a Given Mutation State
#'
#' @description
#' Computes all possible allele configurations corresponding to a given mutation
#' state under the finite-alleles mutation model. Allele states are expanded
#' recursively based on the number of allowable alleles per site or overlapping
#' region as specified by the matrix \code{m}.
#'
#' @details
#' The mutation state vector \code{m_state} encodes the mutation status at each
#' of \eqn{S} sites:
#' \itemize{
#'   \item \code{0}: no mutation at that site
#'   \item \code{1}: a single-site mutation (derived allele)
#'   \item \code{>1}: marker for a multi-site overlapping mutation involving a
#'     block of adjacent sites (from site \eqn{s} to \eqn{e} where
#'     \code{m_state[s] == m_state[e]}).
#' }
#'
#' For each segment or block, the number of distinct allowable alleles is
#' determined by the matrix \code{M}, where:
#' \itemize{
#'   \item \code{m[i, i]} gives the number of alleles with mutation at site
#'         \eqn{i},
#'   \item \code{m[i, j]} (for \eqn{i < j}) gives the number of alleles for a
#'         mutation block spanning sites \eqn{i} to \eqn{j}.
#' }
#'
#' Expansion is handled recursively site by site, and the resulting allele
#' configurations are returned as a matrix, one per row.
#'
#' @param m_state An integer vector of length \eqn{S} representing the mutation
#' state of an allele across \eqn{S} genomic sites.
#' @param m An integer matrix of dimension \eqn{S \times S}, specifying the
#' number of allowable allele states per site or per overlap segment (as
#' described above).
#'
#' @return
#' A matrix of dimension \eqn{K \times S}, where \eqn{K} is the total number of
#' valid allele configurations corresponding to the input mutation state. Each
#' row is one possible allele state.
#' #For the mutation state, returns all possible
all_allele_states <- function(m_state, m) {
  s <- length(m_state)

  alleles <- c()
  s <- 1
  while (s < s + 1) {
    if (m_state[s] == 0) {
      site_possible <- c(0)
      news <- s + 1
      end_site <- s
    } else if (m_state[s] == 1) {
      site_possible <- site_alleles(s, m[s, s])
      news <- s + 1
      end_site <- s
    } else if (m_state[s] > 1) {
      end_site <- max(which(m_state == m_state[s]))
      site_possible <- site_alleles(s * 10 + end_site, m[s, end_site])
      news <- end_site + 1
    }

    p <- length(site_possible)
    if (s == 1) { # first site
      new_matrix <- matrix(0, nrow = p, ncol = (news - 1))
      for (j in 1:p) {
        add_on <- rep(site_possible[j], end_site - s + 1)
        new_matrix[j, ] <- add_on
      }
    } else {
      new_matrix <- matrix(0, nrow = nrow(alleles) * p, ncol = news - 1)
      for (r in seq_len(nrow(alleles))) {
        row <- alleles[r, ]
        for (j in 1:p) {
          add_on <- rep(site_possible[j], end_site - s + 1)
          new_state <- c(row, add_on)
          new_matrix[p * (r - 1) + j, ] <- new_state
        }
      }
    }
    alleles <- new_matrix
    s <- news
  }
  alleles
}





#' Enumerate Possible Allele States Using Leaf-Observed Alleles
#'
#' @description
#' Computes all possible allele configurations corresponding to a given mutation
#' state, but restricts the allowable allele values at each site to only those
#' observed in the leaves (i.e., tip nodes) of the phylogenetic tree. A special
#' wildcard allele is added to model the effect of any allele not observed in
#' the leaves.
#'
#' @details
#' This function is used as part of the approximate likelihood calculation in
#' \code{\link{pq_vectors}}. Instead of expanding all possible alleles at each
#' site or overlapping block (as done in \code{\link{all_allele_states}}),
#' this function considers only alleles already present in the observed leaf
#' data, plus a suffix `"WC"` denoting a wildcard allele (for unseen alleles).
#'
#' The mutation state vector \code{m_state} has the same interpretation as in
#' \code{\link{all_allele_states}}: positions with value \code{0} indicate no
#' mutation, \code{1} indicate a site-specific mutation, and values \code{>1}
#' indicate an overlapping block. The function expands these blocks by sampling
#' allele configurations from the leaf data or assigning a wildcard constant.
#'
#' @param m_state A numeric or character vector of length \eqn{S} representing
#' the mutation state of an allele across \eqn{S} sites.
#' @param leaves A character matrix of dimension \eqn{n \times S}, where each
#' row is a leaf allele state across all sites. Used to restrict the possible
#' allele values at each site to those already observed.
#'
#' @return
#' A character matrix of dimension \eqn{K \times S}, where \eqn{K} is the total
#' number of valid allele configurations based on the mutation state and
#' observed leaves. Each row is an allele state composed from leaf-observed
#' alleles or the wildcard symbol \code{"WC"}.
all_allele_states_leaves <- function(m_state, leaves) {
  s <- length(m_state)

  alleles <- c()
  s <- 1
  while (s < s + 1) {
    if (m_state[s] == 0) {
      site_possible <- c(0)
      news <- s + 1
      end_site <- s
    } else if (m_state[s] == 1) {
      leaves_site <- leaves[, s]
      possible <- c()
      for (j in seq_along(leaves_site)) { # surely a better way to do this
        if (leaves_site[j] != "0") {
          temp <- strsplit(leaves_site[j], ":")[[1]]
          a <- temp[1]
          if (as.numeric(a) <= s) { # make sure it's not an overlap
            possible <- c(possible, leaves_site[j])
          }
        }
      }
      site_possible <- c(unique(possible), paste(s, "WC", sep = ""))
      news <- s + 1
      end_site <- s
    } else if (m_state[s] > 1) {
      end_site <- max(which(m_state == m_state[s]))
      site_possible <- site_alleles(s * 10 + end_site, m_state[s, end_site])

      possible <- c()
      for (j in seq_len(nrow(leaves))) { # surely a better way to do this
        row <- leaves[j, ]
        if (row[s] == row[end_site] && row[s] != "0") {
          possible <- c(possible, row[s])
        }
      }
      m <- 10 * s + end_site
      site_possible <- c(unique(possible), paste(m, "WC", sep = ""))
      news <- end_site + 1
    }

    p <- length(site_possible)
    if (s == 1) { # first site
      new_matrix <- matrix(0, nrow = p, ncol = (news - 1))
      for (j in 1:p) {
        add_on <- rep(site_possible[j], end_site - s + 1)
        new_matrix[j, ] <- add_on
      }
    } else {
      new_matrix <- matrix(0, nrow = nrow(alleles) * p, ncol = news - 1)
      for (r in seq_along(alleles)) {
        row <- alleles[r, ]
        for (j in 1:p) {
          add_on <- rep(site_possible[j], end_site - s + 1)
          new_state <- c(row, add_on)
          new_matrix[p * (r - 1) + j, ] <- new_state
        }
      }
    }
    alleles <- new_matrix
    s <- news
  }
  alleles
}


#' Compute Transition Probabilities Between Allele States Under Finite Alleles
#' Model
#'
#' @description
#' Computes the probability of transitioning from a given parental allele state
#' to a child allele state along an edge of a phylogenetic tree, under a finite
#' alleles mutation model. The calculation accounts for both mutation state
#' transitions (lumped representation) and allele-level resolutions, assuming
#' either uniform mutation probabilities or user-provided mutation schemes.
#'
#' @details
#' The transition probability between two allele states is calculated in two
#' stages:
#' \enumerate{
#'   \item The probability of transitioning between their corresponding lumped
#'    mutation states, obtained using the lumped transition probability matrix
#'    \code{Pt}.
#'   \item The probability of resolving the actual alleles consistent with the
#'    mutation states, based on either:
#'     \itemize{
#'       \item a uniform distribution of possible alleles drawn from matrix
#'        \code{M}, or
#'       \item a user-provided mutation probability structure \code{mu}.
#'     }
#' }
#'
#' If any allele-level mismatch cannot be reconciled by the transition, the
#' function returns zero probability.
#'
#' @param parent A character vector of length \eqn{S} representing the parental
#' allele state across all sites.
#' @param child A character vector of length \eqn{S} representing the child
#' allele state across all sites.
#' @param p_t A numeric matrix representing the transition probabilities between
#' lumped mutation states (typically computed as
#' \eqn{V \cdot \exp(\Lambda t) \cdot V^{-1}} from eigendecomposition of the
#' rate matrix).
#' @param states A named list mapping lumped allele states (as character
#' strings) to integer indices used to access rows/columns of \code{p_t}.
#' @param m An integer matrix where \code{m[i, i]} gives the number of possible
#' alleles for a single-site mutation at site \eqn{i}, and \code{m[i, j]} (for
#' \eqn{i < j}) gives the number of possible alleles for an overlapping
#' segment from site \eqn{i} to \eqn{j}.
#' @param parent_m Optional numeric vector of length \eqn{S}, giving the lumped
#' mutation state of the parent; computed using \code{\link{lumped_state}} if
#' not provided.
#' @param child_m Optional numeric vector of length \eqn{S}, giving the lumped
#' mutation state of the child, computed if missing.
#' @param k Optional integer index for the row of \code{p_t} corresponding to
#' \code{parent_m}; computed from \code{states} if missing.
#' @param j Optional integer index for the column of \code{p_t} corresponding to
#' \code{child_m}; computed if missing.
#' @param mu Optional list of custom mutation probability distributions. If
#' provided, \code{mu[[s]][[1]]} or \code{mu[[s]][[j - s + 1]]} should give a
#' named vector of allele probabilities for single-site or multi-site overlaps,
#' respectively.
#'
#' @return A numeric scalar giving the transition probability from \code{parent}
#' to \code{child} along an edge, based on the model and mutation parameters.
transition_prob_finite_alleles <- function(
  parent,
  child,
  p_t,
  states,
  m,
  parent_m = NULL,
  child_m = NULL,
  k = NULL,
  j = NULL,
  mu = NULL
) {
  s <- nchar(names(states)[1])
  # First compute probability of mutation state transition
  if (is.null(parent_m)) {
    parent_m <- lumped_state(parent)
  }
  if (is.null(child_m)) {
    child_m <- lumped_state(child)
  }
  if (is.null(k)) {
    k <- states[[
      stringi::stri_join(parent_m, collapse = "")
    ]]
  }
  if (is.null(j)) {
    j <- states[[stringi::stri_join(child_m, collapse = "")]]
  }

  prob <- p_t[k, j]


  # If the mutation states are not compatible, zero probability
  if (prob == 0) {
    return(0)
  }
  # still have to check the allele states are consistent and multiply by allele
  # mutation probabilities
  s <- 1
  while (s < s + 1) {
    # If mutation states are the same, allele states should be the same
    if (parent_m[s] == child_m[s]) {
      if (parent[s] != child[s]) {
        return(0)
      }
      s <- s + 1
    } else {
      # If parent has inactive site, but child site is active, not possible
      # Actually, we don't need to include this because it should already be
      # checked by the transition probability of the mutation state

      if (child_m[s] == 1) { # a new single mutation at site s
        if (is.null(mu)) { # Assume all mutations equally likely
          prob <- prob * (1 / m[s, s])
        } else {
          mu_site <- mu[[s]][[1]]
          prob <- prob * mu_site[[child[s]]]
        }
        s <- s + 1
      } else { # in this case child_m[s] > 1 means a simultaneous cut
        # the end site of the mutation
        end_site <- max(which(child_m == child_m[s]))

        if (is.null(mu)) {
          prob <- prob * (1 / m[s, end_site])
        } else {
          mu_site <- mu[[s]][[end_site - s + 1]]
          prob <- prob * mu_site[[child[s]]]
        }
        s <- end_site + 1
      }
    }
  }
  prob
}


#' Piecewise Exponential Trajectory
#'
#' Computes a piecewise function over time \code{t}, useful for simulations
#' along a fixed tree. The function is:
#' \itemize{
#'   \item \code{5} when \code{t <= 0.1}
#'   \item \code{5 * exp(0.5 - 5t)} when \code{0.1 < t < 0.5}
#'   \item \code{0.6766764} when \code{t >= 0.5}
#' }
#'
#' @param time Numeric vector of time points.
#'
#' @return Numeric vector of trajectory values, same length as \code{t}.
#'
#' @examples
#' t <- seq(0, 1, by = 0.01)
#' plot(t, exp_traj(t), type = "l")
#'
#' @export
exp_traj <- function(time) {
  result <- rep(0, length(time))
  result[time <= 0.1] <- 5
  result[time > 0.1 & time < 0.5] <- 5 *
    exp(0.5 - 5 * time[time > 0.1 & time < 0.5]) # 50 times smaller
  result[time >= 0.5] <- 0.6766764
  result
}
