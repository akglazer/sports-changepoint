library(tidyverse)
library(data.table)
library(permutest)

# -----------------------------------------------------------------------------
# Single Changepoint Detection Algorithm for Normally Distributed Data
#
# Description:
#   Implements an algorithm to detect a single changepoint 
#   in a sequence of continuous data using likelihood ratio tests. Optionally 
#   performs sample splitting for inference to reduce overfitting and control 
#   false positives.
#
#   The algorithm searches over possible changepoints in the odd-indexed data
#   (or the full data if sample splitting is disabled), identifies the 
#   location that maximizes the likelihood ratio statistic, and then uses a 
#   two-sample permutation test on the even-indexed data to assess statistical 
#   significance.
#
# Inputs:
#   - data: A numeric vector of values representing a performance 
#           metric over time.
#   - threshold: A numeric value representing the shift parameter in the 
#                hypothesis test.
#   - sample_splitting: Logical flag (default TRUE). If TRUE, uses odd-indexed 
#           observations for changepoint detection and even-indexed observations 
#           for hypothesis testing. If FALSE, uses the entire dataset for both.
#
# Returns:
#   A list with the following elements:
#     - changepoint: Corresponding index in the full data sequence (adjusted).
#     - p_value: P-value from permutation test with shift parameter.
#     - max_diff: Difference in segment means at the detected changepoint.
#
# Notes:
#   - The function excludes the first and last 20 points in the candidate 
#     changepoint search to reduce spurious edge effects.
#   - Log-likelihoods are calculated under a Normal model for continuous data.
# -----------------------------------------------------------------------------
normal_cp <- function(data, threshold = 0,
                      sample_splitting = TRUE){
  # Initialize variables
  max_lr <- -Inf
  changepoint <- NULL
  p_value <- 1
  
  # if sample splitting is true, split data be even/odd indices
  if(sample_splitting){
    # Get odd and even indices
    odd_indices  <- seq(1, length(data), by = 2)
    even_indices <- seq(2, length(data), by = 2)
    
    # Split the data vector by odd and even indices
    data_odd  <- data[odd_indices] # use to find changepoint
    n_odd <- length(data_odd)
    data_even <- data[even_indices] # use to test changepoint
    n_even <- length(data_even)
  } else {
    # if sample splitting is false, use full dataset to find and test cp
    data_odd <- data
    data_even <- data
    n_odd <- length(data)
    n_even <- n_odd
  }
  
  # Iterate over possible changepoints - exclude first 20 observations
  for (k in 20:(n_odd-20)) {
    # Calculate mean/variance of each segment
    # before changepoint
    mu1 <- mean(data_odd[1:k])
    var1 <- sum((data_odd[1:k] - mu1)^2)/k
    # after changepoint
    mu2 <- mean(data_odd[(k+1):n_odd])
    var2 <- sum((data_odd[(k+1):n_odd] - mu2)^2)/(n_odd-k)
    # whole dataset
    mu <- mean(data_odd)
    var <- sum((data_odd - mu)^2)/n_odd
    
    # Compute likelihood ratio
    L1 <- (-k/2)*log(2*pi*var1)
    L2 <- (-(n_odd - k)/2)*log(2*pi*var2)
    L <- (n_odd/2)*log(2*pi*var)
    LR <- 2 * (L1 + L2 + L)
    
    # Update max likelihood ratio and changepoint if LR is greater
    if (is.na(LR) == F & LR > max_lr) {
      max_lr <- LR
      changepoint <- k
      max_diff <- mu1-mu2
    }
  }
  
  # get shift value
  if(max_diff >= 0){
    alt_value <- "greater"
    mu_value <- abs(threshold)
  } else {
    alt_value <- "less"
    mu_value <- -1*abs(threshold)
  }
  
  # Get p-value for changepoint
  p_value <- two_sample(data_even[1:changepoint], 
                        data_even[changepoint:n_even],
                        shift = mu_value,
                        alternative = alt_value,
                        reps = 10^4)
  
  # calculate changepoint in global index
  if(sample_splitting){
    c <- 2*changepoint-1
  } else {
    c <- changepoint
  }
  
  return(list(changepoint = c,
              p_value = p_value, 
              max_diff = max_diff))
}

# -----------------------------------------------------------------------------
# Multiple Changepoint Detection Algorithm
#
# Description:
#   Implements binary segmentation to detect multiple changepoints 
#   in a sequence of continuous data using likelihood ratio tests. Optionally 
#   performs sample splitting for inference to reduce overfitting and control 
#   false positives.
#
# Inputs:
#   - data: A numeric vector of values representing a performance 
#           metric over time.
#   - alpha: significance level (default 0.05) for rejecting the null hypothesis 
#            of no changepoint in a segment.
#   - threshold: A numeric value representing the shift parameter in the 
#                hypothesis test.
#   - sample_splitting: Logical flag (default TRUE). If TRUE, uses odd-indexed 
#           observations for changepoint detection and even-indexed observations 
#           for hypothesis testing. If FALSE, uses the entire dataset for both.
#   - min_seg_len: Minimum segment length to consider for further splitting. 
#                  Defaults to 80 if sample splitting is TRUE, otherwise 40.
#
# Returns:
#   A list with the following elements:
#     - changepoints: A sorted vector of changepoint locations (global indices).
#     - segments: A data frame with columns `segment`, `start`, and `end` 
#                 indicating segment boundaries.
# Example usage:
#   x <- c(rnorm(500, 0, 1), rnorm(500, 1, 1))
#   binary_segment_normal(x)
# -----------------------------------------------------------------------------
binary_segment_normal <- function(data,
                           alpha = 0.05,
                           threshold = 0,
                           sample_splitting = TRUE,
                           min_seg_len = if (sample_splitting) 80L else 40L) {
  
  # recursive helper
  recurse <- function(x, s, e, cps_accum) {
    seg_len <- e - s + 1L
    if (seg_len < min_seg_len) return(cps_accum)
    
    seg <- data[s:e]
    
    res <- normal_cp(seg, threshold = threshold, sample_splitting = sample_splitting)
    
    # stop if no cp found or not significant
    if (is.null(res$changepoint) || is.na(res$p_value) || res$p_value >= alpha) {
      return(cps_accum)
    }
    
    # map segment cp back to global index
    cp_global <- s + res$changepoint - 1L
    
    # guard against degenerate splits
    if (cp_global <= s || cp_global >= e) {
      return(cps_accum)
    }
    
    # accumulate and recurse left/right
    cps_accum <- c(cps_accum, cp_global)
    cps_accum <- recurse(x = data, s = s, e = cp_global, cps_accum = cps_accum)
    cps_accum <- recurse(x = data, s = cp_global + 1L, e = e, cps_accum = cps_accum)
    cps_accum
  }
  
  cps <- sort(unique(recurse(data, 1L, length(data), integer(0))))
  
  # build a segments data.frame for convenience
  if (length(cps) == 0L) {
    segments <- data.frame(segment = 1L, start = 1L, end = length(data))
  } else {
    cuts <- c(1L, cps + 1L, length(data) + 1L)
    segments <- data.frame(
      segment = seq_along(cuts[-1]),
      start   = cuts[-length(cuts)],
      end     = cuts[-1] - 1L
    )
  }
  
  list(changepoints = cps, segments = segments)
}

