# utils_peaks.R
# Robust peak detection helper for 0, 1, or multiple peaks.
# This version prevents "incorrect number of dimensions" when only one peak is detected.

peaks <- function(data, nups, ndowns, minpeakheight, minpeakdistance) {
  data <- as.data.frame(data)
  if (ncol(data) < 2) {
    stop("peaks(): data must have at least two columns: Time and signal.")
  }
  
  time_vals <- suppressWarnings(as.numeric(data[[1]]))
  signal_vals <- suppressWarnings(as.numeric(data[[2]]))
  
  empty_peak_table <- data.frame(
    absolute_amplitude = numeric(0),
    posision_peaks = numeric(0),
    l_inf = numeric(0),
    l_sup = numeric(0)
  )
  
  empty_peak_positions <- data.frame(
    peak_height = numeric(0),
    peak_index = integer(0),
    left_index = integer(0),
    right_index = integer(0)
  )
  
  if (length(signal_vals) < 3 || all(!is.finite(signal_vals))) {
    return(list(p_eak = empty_peak_table, peak = empty_peak_positions))
  }
  
  if (any(!is.finite(signal_vals))) {
    idx <- seq_along(signal_vals)
    good <- is.finite(signal_vals)
    if (sum(good) >= 2) {
      signal_vals[!good] <- stats::approx(
        x = idx[good],
        y = signal_vals[good],
        xout = idx[!good],
        rule = 2
      )$y
    } else {
      return(list(p_eak = empty_peak_table, peak = empty_peak_positions))
    }
  }
  
  peak_mat <- tryCatch(
    pracma::findpeaks(
      x = signal_vals,
      nups = nups,
      ndowns = ndowns,
      minpeakheight = minpeakheight,
      minpeakdistance = minpeakdistance,
      sortstr = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(peak_mat) || length(peak_mat) == 0) {
    return(list(p_eak = empty_peak_table, peak = empty_peak_positions))
  }
  
  peak_mat <- as.matrix(peak_mat)
  if (is.null(dim(peak_mat))) {
    peak_mat <- matrix(peak_mat, nrow = 1)
  }
  if (ncol(peak_mat) < 4) {
    return(list(p_eak = empty_peak_table, peak = empty_peak_positions))
  }
  
  peak_mat <- peak_mat[, 1:4, drop = FALSE]
  peak_mat <- peak_mat[order(peak_mat[, 2]), , drop = FALSE]
  
  peak_positions <- data.frame(
    peak_height = suppressWarnings(as.numeric(peak_mat[, 1])),
    peak_index = suppressWarnings(as.integer(peak_mat[, 2])),
    left_index = suppressWarnings(as.integer(peak_mat[, 3])),
    right_index = suppressWarnings(as.integer(peak_mat[, 4]))
  )
  
  valid <- is.finite(peak_positions$peak_index) &
    peak_positions$peak_index >= 1 &
    peak_positions$peak_index <= length(time_vals) &
    is.finite(peak_positions$left_index) &
    peak_positions$left_index >= 1 &
    peak_positions$left_index <= length(time_vals) &
    is.finite(peak_positions$right_index) &
    peak_positions$right_index >= 1 &
    peak_positions$right_index <= length(time_vals)
  
  peak_positions <- peak_positions[valid, , drop = FALSE]
  
  if (nrow(peak_positions) == 0) {
    return(list(p_eak = empty_peak_table, peak = empty_peak_positions))
  }
  
  peak_table <- data.frame(
    absolute_amplitude = round(peak_positions$peak_height, 3),
    posision_peaks = time_vals[peak_positions$peak_index],
    l_inf = time_vals[peak_positions$left_index],
    l_sup = time_vals[peak_positions$right_index]
  )
  
  list(
    p_eak = peak_table,
    peak = peak_positions
  )
}
