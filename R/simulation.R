#'State Space Hash Table
#'
#'This generates the state space hash table for a given number of sites
#'
#'@param(s) Number of sites
#'
#'@return A hash table mapping state strings to their indices
#'@examples
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
state_space <- function(S) {
  states = state_space_list(S)
  states = relabel(states) #matrix
  h = hash::hash()
  for (i in 1:nrow(states)) {
    state = paste(states[i, ], collapse = "")
    h[[state]] = i
  }
  return(h)
}

#' State Space List
#'
#' The state space contains all possible vectors of length S,
#' where 0 indicates no mutation, 1 indicates mutation just
#' at that site, and larger numbers i are mutations that are
#' shared between multiple adjacent sites
#'
#' @param(S)
#'
#' @return A list of all possible states
#'
#' @examples
#' state_space_list(3)
#' # [[111],[011],[221],[101],[001],[122],[022],[222],[110],[010],[220],[100],[000]]
#' @export
state_space_list <- function(S) {
  #recursive function
  if (S == 1) {
    return(list(c(1), c(0)))
  }
  else {
    SS_1 = state_space_list(S-1)

    new = list()
    for (state in SS_1) {
      new[[length(new) + 1]] = c(1, state)
      new[[length(new) + 1]] = c(0, state)

      if (state[1] == 1) { #can merge
        m = max(state) + 1
        new[[length(new) + 1]] = c(m, m, state[-1])
      }
      if (state[1] > 1) { #1 is already part of an overlap, merge with that
        new[[length(new) + 1]] = c(state[1], state)
      }
    }
    return(new)
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
#' @param(states) A list of states
#'
#' @return A matrix of relabeled states
#'
#' @examples
#' relabel([[111],[011],[221],[101],[001],[122],[022],[222],[110],[010],[220],[100],[000]])
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
relabel <- function(states) {
  L = length(states) #states is a list
  new = matrix(0, nrow = L, ncol = length(states[[1]]))

  for (i in 1:L) {
    s = states[[i]]
    newS = s
    if (max(s) > 1) {
      for (k in 2:max(s)) {
        overlap = which(s == k)
        first = overlap[1]
        newS[overlap] <- first + 1
      }
    }
    new[i, ]  = newS
  }
  return(new)
}

#' State Space Matrix
#'
#' Turns the state space list into a matrix and relabels
#'
#' @param(S) Number of sites
#'
#' @returns An nxS matrix where each row represents a state
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
state_space_matrix <- function(S) {
  states = state_space_list(S)
  return(relabel(states))
}


#' Q Matrix
#'
#' Q is the generator matrix for a continuous time Markov chain over the list
#' of states. Theta is an SxS matrix with theta[i, j] is the rate of a cut from
#' site i to site j. Theta[i, i] is the rate of a single cut at site i.
#'
#'
#' @param(states) A matrix of states (output of state_space_matrix) with length n
#' @param(theta) An SxS matrix of cut rates
#'
#' @returns A nxn generator matarix
#'
#' @export
Q_matrix <- function(states, theta) {
  L = nrow(states)
  S = ncol(states)
  M = matrix(0, nrow = L, ncol = L)

  for (i in 1:L) {
    s1 = states[i, ]
    active = which(s1 == 0)
    if (length(active) > 0) {
      #First, all possible single mutations
      for (k in 1:length(active)) {
        site = active[k]
        s2 = s1
        s2[site] = 1
        #This tells me which index the new state corresponds to,
        #so I can update the appropriate entry in the rate matrix
        #(Like a cheat version of a dictionary)
        j = which(apply(states, 1, function(x) all.equal(x[1:S], s2)) == "TRUE")
        M[i, j] = theta[site, site]

        if (length(active) - k > 0) {
          #now all pairs of mutations
          for (l in (k+1):length(active)) {
            site2 = active[l]
            s2 = s1
            #This mutation assignment is consistent with what we did before
            s2[site:site2] = site + 1 #we should always have site2 > site
            j = which(apply(states, 1, function(x) all.equal(x[1:S], s2)) == "TRUE")
            M[i, j] = theta[site, site2]
          }
        }
      }
    }
  }

  #Add negatives along diagonal so rows sum to 1
  rows = rowSums(M)
  diagonal = diag(rows, nrow = nrow(M), ncol = ncol(M))
  M = M - diagonal
  return(M)
}
