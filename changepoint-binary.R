library(tidyverse)
library(changepoint)
library(permutest)

# -----------------------------------------------------------------------------
# Binary Changepoint Detection Algorithm
#
# Description:
#   Implements an algorithm to detect a single changepoint 
#   in a binary (0/1) sequence of data using likelihood ratio tests. Optionally 
#   performs sample splitting for inference to reduce overfitting and control 
#   false positives.
#
#   The algorithm searches over possible changepoints in the odd-indexed data
#   (or the full data if sample splitting is disabled), identifies the 
#   location that maximizes the likelihood ratio statistic, and then uses a 
#   Fisher's exact test (or z-test when threshold > 0) on the even-indexed data
#   to assess statistical significance.
#
# Inputs:
#   - data: A numeric vector of binary values (0 or 1) representing a performance 
#           metric over time.
#   - threshold: minimum change in proportion (e.g., 0.05 = 5 percentage
#           points) to test for. If 0, Fisher's exact test is used. If > 0, a
#           one-sided unpooled two-proportion z-test of H0: change <= threshold.
#   - sample_splitting: Logical flag (default TRUE). If TRUE, uses odd-indexed 
#           observations for changepoint detection and even-indexed observations 
#           for hypothesis testing. If FALSE, uses the entire dataset for both.
#
# Returns:
#   A list with the following elements:
#     - changepoint: Corresponding index in the full data sequence (adjusted).
#     - p_value: P-value from Fisher's exact test.
#     - max_diff: Difference in segment means (p2 - p1) at the detected changepoint.
#
# Notes:
#   - The function excludes the first and last 20 points in the candidate 
#     changepoint search to reduce spurious edge effects.
#   - Log-likelihoods are calculated under a Bernoulli model for binary data.
# -----------------------------------------------------------------------------
binary_cp <- function(data, threshold = 0, sample_splitting = TRUE){
  # Initialize variables
  max_lr <- -Inf
  changepoint <- NULL
  max_diff <- NA_real_ 
  p_value <- 1
  
  # if sample splitting is true, split data be even/odd indices
  if(sample_splitting){
    # Get odd and even indices
    odd_indices  <- seq(1, length(data), by = 2)
    even_indices <- seq(2, length(data), by = 2)
    
    # Split the data vector by odd and even indices
    data_odd  <- data[odd_indices] # use odd indices to find changepoint
    n_odd <- length(data_odd)
    data_even <- data[even_indices] # use even indices to test changepoint
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
    # Calculate segment probabilities
    p1 <- mean(data_odd[1:k])
    p2 <- mean(data_odd[(k+1):n_odd])
    p <- mean(data_odd)
    
    # Calculate log-likelihoods for segments and entire data
    L1 <- k * (p1 * log(p1) + (1 - p1) * log(1 - p1))
    L2 <- (n_odd - k) * (p2 * log(p2) + (1 - p2) * log(1 - p2))
    L <- n_odd * (p * log(p) + (1 - p) * log(1 - p))
    
    # Compute likelihood ratio
    LR <- 2 * (L1 + L2 - L)
    
    # Update max likelihood ratio and changepoint if LR is greater
    if (is.na(LR) == F & LR > max_lr) {
      max_lr <- LR
      changepoint <- k
      max_diff <- p2-p1
    }
  }
  
  # test on even-indexed timepoints
  before <- data_even[1:changepoint]
  after  <- data_even[(changepoint + 1):n_even]
  dir <- if (max_diff >= 0) 1 else -1 
  
  
  if (threshold == 0) {
    # Fisher's exact test
    tab <- matrix(c(sum(before == 1), sum(before == 0),
                    sum(after == 1),  sum(after == 0)), nrow = 2)
    alt <- if (dir == 1) "less" else "greater"
    p_value <- fisher.test(tab, alternative = alt)$p.value
  } else {
    # one-sided unpooled two-proportion z-test of H0: change <= threshold
    pb <- mean(before); pa <- mean(after)
    se <- sqrt(pb * (1 - pb) / length(before) + pa * (1 - pa) / length(after))
    stat <- dir * (pa - pb) - threshold
    p_value <- if (se == 0) as.numeric(stat <= 0) else pnorm(stat / se, lower.tail = FALSE)
  }
  
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
#   in a binary (0/1) sequence of data using likelihood ratio tests. Optionally 
#   performs sample splitting for inference to reduce overfitting and control 
#   false positives.
#
# Inputs:
#   - data: A numeric vector of binary values (0 or 1) representing a performance 
#           metric over time.
#   - alpha: significance level (default 0.05) for rejecting the null hypothesis 
#            of no changepoint in a segment.
#   - threshold: minimum change in proportion (e.g., 0.05 = 5 percentage
#           points) to test for.
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
#   x <- c(rbinom(500, 1, 0.3), rbinom(500, 1, 0.4))
#   binary_segment(x)
# -----------------------------------------------------------------------------
binary_segment <- function(data,
                           alpha = 0.05,
                           threshold = 0,
                           sample_splitting = TRUE,
                           min_seg_len = if (sample_splitting) 80L else 40L) {
  
  # recursive helper for binary segmentation
  recurse <- function(x, s, e, cps_accum) {
    seg_len <- e - s + 1L
    if (seg_len < min_seg_len) return(cps_accum)
    
    seg <- data[s:e]
    
    res <- binary_cp(seg, threshold = threshold,
                     sample_splitting = sample_splitting)
    
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

