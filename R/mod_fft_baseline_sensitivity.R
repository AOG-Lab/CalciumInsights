# Dependencies are declared in DESCRIPTION/NAMESPACE; helpers are loaded through DESCRIPTION Collate.
# R/mod_fft_baseline_sensitivity.R
# ============================================================
# FFT-based Denoising Data App
# - Preserves the original denoising workflow
# - Replaces LOESS with FFT low-pass smoothing
# - Adds Fourier explicit reconstruction formula using kept frequencies
# - Uses helper functions saved in the same folder
# - FWHP and FWHM are frozen from the original FFT/local-minimum workflow
#   and are not recalculated or redefined by single or double sigmoid fits.
# - Event-specific baseline left/right minima remain the original calcium-transient
#   FFT minima even when single and double sigmoid fits are activated.
# ============================================================

# ============================================================
# Robust peak wrapper for 0, 1, or multiple detected peaks
# ============================================================
# Why this is needed:
# The original helper in utils_peaks.R can lose its two-dimensional structure
# when only one peak is detected. Later code expects peak[, 2], peak[, 3], etc.
# If the object has become a vector, R throws:
#   Error: incorrect number of dimensions
# This replacement keeps the output as a data frame with the same 4-column
# order for all cases: 0 peaks, 1 peak, or multiple peaks.
# ============================================================
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
  
  # Force matrix/data-frame shape. This is the key fix for single-transient data.
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
  
  # Keep only valid indices.
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
  
  # Return names and column order compatible with the original helper.
  list(
    p_eak = peak_table,
    peak = peak_positions
  )
}

# Defensive normalizer for older cached peak objects or external helper outputs.
normalize_fft_peak_positions <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(data.frame(
      peak_height = numeric(0),
      peak_index = integer(0),
      left_index = integer(0),
      right_index = integer(0)
    ))
  }
  x <- as.data.frame(x)
  if (ncol(x) < 4) {
    return(data.frame(
      peak_height = numeric(0),
      peak_index = integer(0),
      left_index = integer(0),
      right_index = integer(0)
    ))
  }
  x <- x[, 1:4, drop = FALSE]
  names(x)[1:4] <- c("peak_height", "peak_index", "left_index", "right_index")
  x$peak_height <- suppressWarnings(as.numeric(x$peak_height))
  x$peak_index <- suppressWarnings(as.integer(x$peak_index))
  x$left_index <- suppressWarnings(as.integer(x$left_index))
  x$right_index <- suppressWarnings(as.integer(x$right_index))
  x <- x[is.finite(x$peak_index) & is.finite(x$left_index) & is.finite(x$right_index), , drop = FALSE]
  x
}

normalize_fft_peak_table <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(data.frame(
      absolute_amplitude = numeric(0),
      posision_peaks = numeric(0),
      l_inf = numeric(0),
      l_sup = numeric(0)
    ))
  }
  x <- as.data.frame(x)
  if (ncol(x) < 4) {
    return(data.frame(
      absolute_amplitude = numeric(0),
      posision_peaks = numeric(0),
      l_inf = numeric(0),
      l_sup = numeric(0)
    ))
  }
  x <- x[, 1:4, drop = FALSE]
  names(x)[1:4] <- c("absolute_amplitude", "posision_peaks", "l_inf", "l_sup")
  x$absolute_amplitude <- suppressWarnings(as.numeric(x$absolute_amplitude))
  x$posision_peaks <- suppressWarnings(as.numeric(x$posision_peaks))
  x$l_inf <- suppressWarnings(as.numeric(x$l_inf))
  x$l_sup <- suppressWarnings(as.numeric(x$l_sup))
  x
}

# =========================
# FFT low-pass smoothing
# =========================
fft_lowpass_smooth_trace <- function(time, signal, f = 0.20, keep_mean = TRUE) {
  signal <- as.numeric(signal)
  n <- length(signal)
  stopifnot(length(time) == n)
  stopifnot(n >= 4)
  
  X <- fft(signal)
  n_pos <- floor(n / 2) + 1
  k_keep <- max(1L, floor(f * n_pos))
  neg_count <- max(0L, k_keep - 1L)
  neg_idx <- if (neg_count > 0L) seq.int(from = n - neg_count + 1L, to = n) else integer(0)
  keep_idx <- unique(c(seq.int(1L, k_keep), neg_idx))
  
  Xf <- complex(length = n)
  Xf[keep_idx] <- X[keep_idx]
  if (!isTRUE(keep_mean)) Xf[1] <- 0
  
  smooth <- Re(fft(Xf, inverse = TRUE) / n)
  
  list(
    data = data.frame(Time = time, signal = smooth),
    fft = X,
    fft_filt = Xf,
    k_keep = k_keep,
    n_pos = n_pos,
    keep_idx = keep_idx
  )
}

safe_fft_fraction <- function(f, fallback = 0.20) {
  f <- suppressWarnings(as.numeric(f))
  fallback <- suppressWarnings(as.numeric(fallback))
  if (length(fallback) == 0 || !is.finite(fallback[1])) fallback <- 0.20
  if (length(f) == 0 || !is.finite(f[1])) f <- fallback[1]
  f <- f[1]
  max(0.01, min(1, f))
}


# =========================
# Polynomial detrending for calcium traces
# =========================
remove_polynomial_trend_trace <- function(time,
                                          signal,
                                          degree = 2,
                                          center_scale_t = TRUE) {
  signal <- suppressWarnings(as.numeric(signal))
  time <- suppressWarnings(as.numeric(time))
  n <- length(signal)
  
  if (length(time) != n) {
    stop("time and signal must have the same length.")
  }
  
  if (n < degree + 2) {
    stop("Not enough observations to fit the selected polynomial degree.")
  }
  
  df <- data.frame(
    Time = time,
    signal_original = signal,
    t = seq_len(n)
  )
  
  if (isTRUE(center_scale_t)) {
    df$t_used <- as.numeric(scale(df$t))
  } else {
    df$t_used <- df$t
  }
  
  fit <- stats::lm(signal_original ~ poly(t_used, degree, raw = TRUE), data = df)
  df$trend <- as.numeric(stats::predict(fit, newdata = df))
  df$signal_detrended <- df$signal_original - df$trend
  
  list(
    data = data.frame(
      Time = df$Time,
      signal = df$signal_detrended
    ),
    trend_data = data.frame(
      Time = df$Time,
      signal_original = df$signal_original,
      trend = df$trend,
      signal_detrended = df$signal_detrended
    ),
    fit = fit,
    degree = degree,
    center_scale_t = center_scale_t
  )
}



# =========================
# Polynomial trend model information
# =========================
build_polynomial_trend_info <- function(trend_res, digits = 6) {
  if (is.null(trend_res) || is.null(trend_res$fit)) {
    return(list(
      info_text = "Trend removal is not selected.",
      coef_table = data.frame()
    ))
  }
  
  fit <- trend_res$fit
  degree <- trend_res$degree
  center_scale_t <- isTRUE(trend_res$center_scale_t)
  n_obs <- length(stats::residuals(fit))
  coefs <- stats::coef(fit)
  fit_sum <- summary(fit)
  
  fmt <- function(x) {
    format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE, scientific = FALSE)
  }
  
  # Build equation using the internal time variable used by the model.
  t_name <- if (center_scale_t) "t_scaled" else "t"
  equation_terms <- c(fmt(coefs[1]))
  
  if (degree >= 1) {
    for (i in seq_len(degree)) {
      beta <- coefs[i + 1]
      sign_txt <- ifelse(beta >= 0, " + ", " - ")
      power_txt <- if (i == 1) t_name else paste0(t_name, "^", i)
      equation_terms <- c(equation_terms, paste0(sign_txt, fmt(abs(beta)), "*", power_txt))
    }
  }
  
  equation_txt <- paste0("m_hat(t) = ", paste0(equation_terms, collapse = ""))
  residual_equation <- "detrended_trace(t) = original_trace(t) - m_hat(t)"
  
  coef_mat <- as.data.frame(fit_sum$coefficients)
  coef_mat$Term <- rownames(coef_mat)
  rownames(coef_mat) <- NULL
  coef_mat <- coef_mat[, c("Term", "Estimate", "Std. Error", "t value", "Pr(>|t|)")]
  names(coef_mat) <- c("Term", "Estimate", "Std_Error", "t_value", "p_value")
  
  coef_mat$Estimate <- round(coef_mat$Estimate, digits)
  coef_mat$Std_Error <- round(coef_mat$Std_Error, digits)
  coef_mat$t_value <- round(coef_mat$t_value, digits)
  coef_mat$p_value <- signif(coef_mat$p_value, 4)
  
  info_text <- paste0(
    "Fitted polynomial trend model\n",
    "----------------------------------------\n",
    "Polynomial degree: ", degree, "\n",
    "Internal time variable: ", t_name, "\n",
    "Center and scale t: ", ifelse(center_scale_t, "Yes", "No"), "\n",
    "Number of observations: ", n_obs, "\n\n",
    "Trend equation:\n",
    equation_txt, "\n\n",
    "Corrected calcium trace:\n",
    residual_equation, "\n\n",
    "Model fit:\n",
    "R-squared: ", fmt(fit_sum$r.squared), "\n",
    "Adjusted R-squared: ", fmt(fit_sum$adj.r.squared), "\n",
    "Residual standard error: ", fmt(fit_sum$sigma), "\n"
  )
  
  list(
    info_text = info_text,
    coef_table = coef_mat
  )
}

# =========================
# Explicit sine-cosine formula
# =========================
build_fft_formula <- function(X, n, k_keep, keep_mean = TRUE, digits = 3) {
  fmt_num <- function(z, digits = 3) format(round(z, digits), nsmall = digits, trim = TRUE, scientific = FALSE)
  
  terms <- character(0)
  
  # Constant term (DC)
  a0 <- if (isTRUE(keep_mean)) Re(X[1]) / n else 0
  terms <- c(terms, fmt_num(a0, digits))
  
  k_max <- k_keep - 1L
  if (k_max >= 1L) {
    for (k in 1:k_max) {
      # Special Nyquist case when n is even and k = n/2
      if (n %% 2 == 0 && k == n / 2) {
        ak <- Re(X[k + 1]) / n
        sign_cos <- if (ak >= 0) " + " else " - "
        cos_term <- paste0(sign_cos, fmt_num(abs(ak), digits),
                           "*cos(2*pi*", k, "*(t-1)/", n, ")")
        terms <- c(terms, cos_term)
      } else {
        ak <- 2 * Re(X[k + 1]) / n
        bk <- -2 * Im(X[k + 1]) / n
        
        sign_cos <- if (ak >= 0) " + " else " - "
        sign_sin <- if (bk >= 0) " + " else " - "
        
        cos_term <- paste0(sign_cos, fmt_num(abs(ak), digits),
                           "*cos(2*pi*", k, "*(t-1)/", n, ")")
        sin_term <- paste0(sign_sin, fmt_num(abs(bk), digits),
                           "*sin(2*pi*", k, "*(t-1)/", n, ")")
        
        terms <- c(terms, cos_term, sin_term)
      }
    }
  }
  
  formula_txt <- paste0("x_hat[t] = ", paste0(terms, collapse = ""))
  formula_txt
}

# =========================
# Example dataset fallback
# =========================
example_calcium_data <- function() {
  set.seed(123)
  t <- seq(0, 100, by = 0.5)
  mk_signal <- function(shift = 0, noise = 0.15) {
    base <- 0.2 * sin(2 * pi * t / 35) + 0.1 * sin(2 * pi * t / 9)
    peaks <- exp(-0.5 * ((t - (20 + shift)) / 1.8)^2) +
      1.2 * exp(-0.5 * ((t - (48 + shift)) / 2.3)^2) +
      0.9 * exp(-0.5 * ((t - (75 + shift)) / 2.1)^2)
    base + peaks + rnorm(length(t), 0, noise) + 0.5
  }
  data.frame(
    Time = t,
    ROI_1 = mk_signal(0),
    ROI_2 = mk_signal(3),
    ROI_3 = mk_signal(-2),
    ROI_4 = mk_signal(5)
  )
}



# =========================
# Percentile baselines for long-term calcium transients
# =========================

# Option 6: rolling percentile baseline
# For each time point, a window centered at that time point is used.
# Therefore, the baseline can change at every time point.
compute_rolling_percentile_baseline <- function(data_smoothed,
                                                window_size = 600,
                                                percentile = 2) {
  time_vals <- suppressWarnings(as.numeric(data_smoothed$Time))
  signal_vals <- suppressWarnings(as.numeric(data_smoothed$signal))
  
  window_size <- suppressWarnings(as.numeric(window_size))
  percentile <- suppressWarnings(as.numeric(percentile))
  
  if (length(time_vals) == 0 || length(signal_vals) == 0 || all(!is.finite(time_vals))) {
    return(data.frame(Time = numeric(0), Baseline = numeric(0)))
  }
  
  if (!is.finite(window_size) || window_size <= 0) {
    window_size <- diff(range(time_vals, na.rm = TRUE))
  }
  if (!is.finite(window_size) || window_size <= 0) {
    window_size <- 1
  }
  
  if (!is.finite(percentile)) percentile <- 2
  percentile <- max(0, min(100, percentile))
  
  time_min <- min(time_vals, na.rm = TRUE)
  time_max <- max(time_vals, na.rm = TRUE)
  half_window <- window_size / 2
  
  baseline_vals <- rep(NA_real_, length(signal_vals))
  window_start_vals <- rep(NA_real_, length(signal_vals))
  window_end_vals <- rep(NA_real_, length(signal_vals))
  
  for (i in seq_along(time_vals)) {
    t0 <- time_vals[i]
    xmin <- max(t0 - half_window, time_min)
    xmax <- min(t0 + half_window, time_max)
    
    idx <- which(time_vals >= xmin & time_vals <= xmax)
    
    if (length(idx) > 0) {
      baseline_vals[i] <- as.numeric(stats::quantile(
        signal_vals[idx],
        probs = percentile / 100,
        na.rm = TRUE,
        type = 7
      ))
      window_start_vals[i] <- xmin
      window_end_vals[i] <- xmax
    }
  }
  
  missing_idx <- which(!is.finite(baseline_vals))
  if (length(missing_idx) > 0) {
    fallback <- as.numeric(stats::quantile(
      signal_vals,
      probs = percentile / 100,
      na.rm = TRUE,
      type = 7
    ))
    baseline_vals[missing_idx] <- fallback
  }
  
  data.frame(
    Time = time_vals,
    Baseline = baseline_vals,
    Window_Start = window_start_vals,
    Window_End = window_end_vals
  )
}

# Option 7: fixed-window percentile baseline
# For each fixed non-overlapping window, the percentile is calculated once.
# Therefore, all points inside the same window share one baseline value.
compute_windowed_percentile_baseline <- function(data_smoothed,
                                                 window_size = 600,
                                                 percentile = 2) {
  time_vals <- suppressWarnings(as.numeric(data_smoothed$Time))
  signal_vals <- suppressWarnings(as.numeric(data_smoothed$signal))
  
  window_size <- suppressWarnings(as.numeric(window_size))
  percentile <- suppressWarnings(as.numeric(percentile))
  
  if (length(time_vals) == 0 || length(signal_vals) == 0 || all(!is.finite(time_vals))) {
    return(data.frame(Time = numeric(0), Baseline = numeric(0)))
  }
  
  if (!is.finite(window_size) || window_size <= 0) {
    window_size <- diff(range(time_vals, na.rm = TRUE))
  }
  if (!is.finite(window_size) || window_size <= 0) {
    window_size <- 1
  }
  
  if (!is.finite(percentile)) percentile <- 2
  percentile <- max(0, min(100, percentile))
  
  time_min <- min(time_vals, na.rm = TRUE)
  time_max <- max(time_vals, na.rm = TRUE)
  
  baseline_vals <- rep(NA_real_, length(signal_vals))
  window_id <- rep(NA_integer_, length(signal_vals))
  window_start_vals <- rep(NA_real_, length(signal_vals))
  window_end_vals <- rep(NA_real_, length(signal_vals))
  
  window_starts <- seq(from = time_min, to = time_max, by = window_size)
  window_starts <- window_starts[window_starts < time_max]
  
  if (length(window_starts) == 0) {
    window_starts <- time_min
  }
  
  for (i in seq_along(window_starts)) {
    xmin <- window_starts[i]
    xmax <- min(xmin + window_size, time_max)
    
    if (i < length(window_starts)) {
      idx <- which(time_vals >= xmin & time_vals < xmax)
    } else {
      idx <- which(time_vals >= xmin & time_vals <= xmax)
    }
    
    if (length(idx) > 0) {
      b <- as.numeric(stats::quantile(
        signal_vals[idx],
        probs = percentile / 100,
        na.rm = TRUE,
        type = 7
      ))
      
      baseline_vals[idx] <- b
      window_id[idx] <- i
      window_start_vals[idx] <- xmin
      window_end_vals[idx] <- xmax
    }
  }
  
  missing_idx <- which(!is.finite(baseline_vals))
  if (length(missing_idx) > 0) {
    fallback <- as.numeric(stats::quantile(
      signal_vals,
      probs = percentile / 100,
      na.rm = TRUE,
      type = 7
    ))
    baseline_vals[missing_idx] <- fallback
  }
  
  data.frame(
    Time = time_vals,
    Baseline = baseline_vals,
    Window_ID = window_id,
    Window_Start = window_start_vals,
    Window_End = window_end_vals
  )
}

get_baseline_info <- function(data_smoothed,
                              baseline_mode = 2,
                              lim_inf = 0,
                              lim_sup = 20,
                              own_baseline = 0,
                              initial_fraction = 0.10,
                              minimal_fraction = 0.10,
                              time_start_increasin_peak = NULL,
                              rolling_window = 600,
                              rolling_percentile = 2) {
  signal_vals <- suppressWarnings(as.numeric(data_smoothed$signal))
  time_vals <- suppressWarnings(as.numeric(data_smoothed$Time))
  baseline_trace <- NULL
  
  # Shiny can evaluate reactives before all inputs exist. These checks prevent
  # errors such as: Error in if: argument is of length zero.
  if (is.null(baseline_mode) || length(baseline_mode) == 0 || is.na(baseline_mode[1])) {
    baseline_mode <- "1"
  } else {
    baseline_mode <- as.character(baseline_mode[1])
  }
  
  safe_mean <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (length(x) == 0 || all(!is.finite(x))) return(0)
    mean(x, na.rm = TRUE)
  }
  
  safe_scalar <- function(x, fallback = safe_mean(signal_vals)) {
    x <- suppressWarnings(as.numeric(x))
    if (length(x) == 0 || !is.finite(x[1])) return(fallback)
    x[1]
  }
  
  scalar <- switch(
    baseline_mode,
    "1" = 0,
    "2" = {
      # Standard definition: baseline estimated from the initial fraction of the trace.
      frac <- safe_scalar(initial_fraction, fallback = 0.10)
      frac <- max(0.01, min(1, frac))
      n_initial <- max(1L, floor(length(signal_vals) * frac))
      safe_mean(signal_vals[seq_len(n_initial)])
    },
    "3" = {
      # Low-fluorescence region baseline: mean of the lowest signal values.
      frac <- safe_scalar(minimal_fraction, fallback = 0.10)
      frac <- max(0.01, min(1, frac))
      y_sorted <- sort(signal_vals[is.finite(signal_vals)])
      if (length(y_sorted) == 0) {
        safe_mean(signal_vals)
      } else {
        n_low <- max(1L, floor(length(y_sorted) * frac))
        safe_mean(y_sorted[seq_len(n_low)])
      }
    },
    "4" = safe_scalar(own_baseline, fallback = safe_mean(signal_vals)),
    "5" = {
      if (length(signal_vals) == 0 || all(!is.finite(signal_vals))) 0 else min(signal_vals, na.rm = TRUE)
    },
    "6" = {
      baseline_trace <- compute_rolling_percentile_baseline(
        data_smoothed = data_smoothed,
        window_size = safe_scalar(rolling_window, fallback = 600),
        percentile = safe_scalar(rolling_percentile, fallback = 2)
      )
      safe_scalar(stats::median(baseline_trace$Baseline, na.rm = TRUE), fallback = safe_mean(signal_vals))
    },
    "7" = {
      baseline_trace <- compute_windowed_percentile_baseline(
        data_smoothed = data_smoothed,
        window_size = safe_scalar(rolling_window, fallback = 600),
        percentile = safe_scalar(rolling_percentile, fallback = 2)
      )
      safe_scalar(stats::median(baseline_trace$Baseline, na.rm = TRUE), fallback = safe_mean(signal_vals))
    },
    "8" = {
      # Event-specific baseline is finalized after peak detection in peaks_plot(),
      # because it needs peak positions.
      safe_mean(signal_vals)
    },
    "9" = {
      if (length(signal_vals) == 0 || all(!is.finite(signal_vals))) 0 else median(signal_vals, na.rm = TRUE)
    },
    "10" = safe_mean(signal_vals),
    "11" = {
      # Constant baseline: mean of the calcium trace from trace start to Time_Onset.
      onset_times <- numeric(0)
      if (!is.null(time_start_increasin_peak) &&
          is.data.frame(time_start_increasin_peak) &&
          "Time" %in% names(time_start_increasin_peak)) {
        onset_times <- suppressWarnings(as.numeric(time_start_increasin_peak$Time))
        onset_times <- onset_times[is.finite(onset_times)]
      }
      if (length(onset_times) > 0) {
        time_onset <- min(onset_times, na.rm = TRUE)
        idx_onset <- which(time_vals <= time_onset)
        if (length(idx_onset) > 0) {
          baseline_trace <- data.frame(Time = time_vals, Baseline = rep(safe_mean(signal_vals[idx_onset]), length(time_vals)))
          safe_mean(signal_vals[idx_onset])
        } else {
          safe_mean(signal_vals)
        }
      } else {
        safe_mean(signal_vals)
      }
    },
    safe_mean(signal_vals)
  )
  
  scalar <- safe_scalar(scalar, fallback = safe_mean(signal_vals))
  list(scalar = scalar, trace = baseline_trace, event_baselines = NULL, mode = baseline_mode)
}

baseline_at_times <- function(time_points, baseline_info) {
  time_points <- as.numeric(time_points)
  
  # Event-specific baseline support:
  # When each detected peak has its own local baseline, return the Baseline of the nearest peak.
  if (!is.null(baseline_info$event_baselines) &&
      is.data.frame(baseline_info$event_baselines) &&
      nrow(baseline_info$event_baselines) > 0 &&
      "Peak_Time" %in% names(baseline_info$event_baselines) &&
      "Baseline_at_peak" %in% names(baseline_info$event_baselines)) {
    
    peak_times <- suppressWarnings(as.numeric(baseline_info$event_baselines$Peak_Time))
    peak_baselines <- suppressWarnings(as.numeric(baseline_info$event_baselines$Baseline_at_peak))
    
    return(vapply(time_points, function(tt) {
      if (!is.finite(tt) || length(peak_times) == 0 || all(!is.finite(peak_times))) {
        return(baseline_info$scalar)
      }
      j <- which.min(abs(peak_times - tt))
      if (length(j) == 0 || !is.finite(peak_baselines[j])) baseline_info$scalar else peak_baselines[j]
    }, numeric(1)))
  }
  
  if (!is.null(baseline_info$trace) && nrow(baseline_info$trace) > 0) {
    return(approx(
      x = baseline_info$trace$Time,
      y = baseline_info$trace$Baseline,
      xout = time_points,
      rule = 2
    )$y)
  }
  rep(baseline_info$scalar, length(time_points))
}

find_left_minimum_smoothed <- function(signal_vals, peak_idx) {
  n <- length(signal_vals)
  peak_idx <- suppressWarnings(as.integer(peak_idx))
  if (!is.finite(peak_idx) || peak_idx <= 1 || peak_idx > n) return(NA_integer_)
  
  y <- signal_vals[seq_len(peak_idx)]
  if (length(y) < 3) return(which.min(y))
  
  local_min_idx <- which(
    y[2:(length(y) - 1)] <= y[1:(length(y) - 2)] &
      y[2:(length(y) - 1)] <= y[3:length(y)]
  ) + 1
  
  if (length(local_min_idx) > 0) tail(local_min_idx, 1) else which.min(y)
}

find_right_minimum_smoothed <- function(signal_vals, peak_idx) {
  n <- length(signal_vals)
  peak_idx <- suppressWarnings(as.integer(peak_idx))
  if (!is.finite(peak_idx) || peak_idx < 1 || peak_idx >= n) return(NA_integer_)
  
  y <- signal_vals[peak_idx:n]
  if (length(y) < 3) {
    local_rel_idx <- which.min(y)
  } else {
    local_min_idx <- which(
      y[2:(length(y) - 1)] <= y[1:(length(y) - 2)] &
        y[2:(length(y) - 1)] <= y[3:length(y)]
    ) + 1
    
    local_rel_idx <- if (length(local_min_idx) > 0) local_min_idx[1] else which.min(y)
  }
  
  peak_idx + local_rel_idx - 1
}

compute_event_specific_minima_baseline <- function(data_smoothed,
                                                   table_positions_peaks,
                                                   time_start_increasin_peak = NULL,
                                                   baseline_smoothed = NULL,
                                                   baseline_fft_fraction = NA_real_) {
  # data_smoothed is the main FFT-smoothed calcium trace used for peak detection,
  # AUC, FWHM and plotting.
  time_vals <- suppressWarnings(as.numeric(data_smoothed[[1]]))
  signal_vals <- suppressWarnings(as.numeric(data_smoothed[[2]]))
  n_signal <- length(signal_vals)
  
  # baseline_smoothed is a second FFT-smoothed trace used only to find the
  # left/right minima that define the event-specific baseline.
  # This lets the user tune the event-specific baseline smoothing independently
  # from the main FFT smoothing selected in section 4.
  if (is.null(baseline_smoothed)) {
    baseline_smoothed <- data_smoothed
  }
  
  baseline_time_vals <- suppressWarnings(as.numeric(baseline_smoothed[[1]]))
  baseline_signal_vals <- suppressWarnings(as.numeric(baseline_smoothed[[2]]))
  
  baseline_fft_fraction <- suppressWarnings(as.numeric(baseline_fft_fraction))
  baseline_fft_fraction <- if (length(baseline_fft_fraction) == 0 ||
                               !is.finite(baseline_fft_fraction[1])) {
    NA_real_
  } else {
    baseline_fft_fraction[1]
  }
  
  # Safety: if the baseline-specific smoothed trace is not aligned with the main
  # trace, interpolate it onto the main Time grid. Peak indices are then valid.
  if (length(baseline_time_vals) != n_signal ||
      length(baseline_signal_vals) != n_signal ||
      all(!is.finite(baseline_time_vals)) ||
      all(!is.finite(baseline_signal_vals))) {
    if (length(baseline_time_vals) > 1 &&
        length(baseline_signal_vals) > 1 &&
        any(is.finite(baseline_time_vals)) &&
        any(is.finite(baseline_signal_vals))) {
      m_align <- min(length(baseline_time_vals), length(baseline_signal_vals))
      baseline_time_tmp <- baseline_time_vals[seq_len(m_align)]
      baseline_signal_tmp <- baseline_signal_vals[seq_len(m_align)]
      ok <- is.finite(baseline_time_tmp) & is.finite(baseline_signal_tmp)
      if (sum(ok) >= 2) {
        baseline_signal_vals <- stats::approx(
          x = baseline_time_tmp[ok],
          y = baseline_signal_tmp[ok],
          xout = time_vals,
          rule = 2
        )$y
      } else {
        baseline_signal_vals <- signal_vals
      }
      baseline_time_vals <- time_vals
    } else {
      baseline_time_vals <- time_vals
      baseline_signal_vals <- signal_vals
    }
  }
  
  # Fill missing values in the baseline-specific smoothed trace.
  missing_baseline <- which(!is.finite(baseline_signal_vals))
  if (length(missing_baseline) > 0) {
    ok <- which(is.finite(baseline_signal_vals))
    if (length(ok) >= 2) {
      baseline_signal_vals[missing_baseline] <- stats::approx(
        x = ok,
        y = baseline_signal_vals[ok],
        xout = missing_baseline,
        rule = 2
      )$y
    } else {
      baseline_signal_vals[missing_baseline] <- signal_vals[missing_baseline]
    }
  }
  
  fallback_scalar <- if (length(signal_vals) > 0 && any(is.finite(signal_vals))) {
    mean(signal_vals, na.rm = TRUE)
  } else {
    0
  }
  
  empty <- list(
    scalar = fallback_scalar,
    trace = data.frame(Time = numeric(0), Baseline = numeric(0)),
    event_baselines = data.frame(),
    baseline_smoothed = data.frame(Time = baseline_time_vals, signal = baseline_signal_vals),
    baseline_fft_fraction = baseline_fft_fraction,
    mode = "8"
  )
  
  if (is.null(table_positions_peaks) || NROW(table_positions_peaks) == 0) {
    return(empty)
  }
  
  peak_idx <- suppressWarnings(as.integer(table_positions_peaks[, 2]))
  valid_idx <- is.finite(peak_idx) & peak_idx >= 1 & peak_idx <= length(signal_vals)
  peak_idx <- peak_idx[valid_idx]
  if (length(peak_idx) == 0) return(empty)
  
  peak_idx <- sort(unique(peak_idx))
  peak_times <- time_vals[peak_idx]
  
  event_list <- lapply(seq_along(peak_idx), function(i) {
    pidx <- peak_idx[i]
    
    # Minima are found on the baseline-specific FFT-smoothed trace.
    left_idx <- find_left_minimum_smoothed(baseline_signal_vals, pidx)
    right_idx <- find_right_minimum_smoothed(baseline_signal_vals, pidx)
    
    if (!is.finite(left_idx) || !is.finite(right_idx) || right_idx <= left_idx) {
      return(NULL)
    }
    
    left_time <- time_vals[left_idx]
    right_time <- time_vals[right_idx]
    left_value <- baseline_signal_vals[left_idx]
    right_value <- baseline_signal_vals[right_idx]
    peak_time <- time_vals[pidx]
    peak_value <- signal_vals[pidx]
    
    # Local baseline is defined by the left and right minima from the
    # baseline-specific FFT-smoothed trace.
    baseline_at_peak <- approx(
      x = c(left_time, right_time),
      y = c(left_value, right_value),
      xout = peak_time,
      rule = 2
    )$y
    
    # ------------------------------------------------------------
    # Corrected plotting/AUC interval:
    # Instead of filling everything from left minimum to right minimum,
    # calculate the baseline across that interval and keep only the
    # continuous positive component that contains the peak.
    #
    # This prevents green area from appearing before/after the transient
    # and prevents the baseline line from visually connecting transients.
    # ------------------------------------------------------------
    idx_full <- seq(left_idx, right_idx)
    
    baseline_full <- approx(
      x = c(left_time, right_time),
      y = c(left_value, right_value),
      xout = time_vals[idx_full],
      rule = 2
    )$y
    
    # The event-specific baseline is estimated from baseline_signal_vals,
    # but the biological signal/AUC is evaluated on the main smoothed trace.
    above <- signal_vals[idx_full] > baseline_full
    
    peak_rel <- which(idx_full == pidx)
    if (length(peak_rel) == 0) return(NULL)
    
    if (!isTRUE(above[peak_rel])) {
      # If the peak is not above the baseline, this event should not produce
      # event-specific AUC.
      return(NULL)
    }
    
    left_rel <- peak_rel
    while (left_rel > 1 && isTRUE(above[left_rel - 1])) {
      left_rel <- left_rel - 1
    }
    
    right_rel <- peak_rel
    while (right_rel < length(above) && isTRUE(above[right_rel + 1])) {
      right_rel <- right_rel + 1
    }
    
    plot_left_idx <- idx_full[left_rel]
    plot_right_idx <- idx_full[right_rel]
    
    # Additional safety clipping by midpoints between neighboring peaks.
    if (i > 1 && is.finite(peak_times[i - 1])) {
      prev_mid <- (peak_times[i - 1] + peak_time) / 2
      plot_left_idx <- max(plot_left_idx, which.min(abs(time_vals - prev_mid)))
    }
    
    if (i < length(peak_times) && is.finite(peak_times[i + 1])) {
      next_mid <- (peak_time + peak_times[i + 1]) / 2
      plot_right_idx <- min(plot_right_idx, which.min(abs(time_vals - next_mid)))
    }
    
    if (!is.finite(plot_left_idx) ||
        !is.finite(plot_right_idx) ||
        plot_right_idx <= plot_left_idx) {
      return(NULL)
    }
    
    plot_left_time <- time_vals[plot_left_idx]
    plot_right_time <- time_vals[plot_right_idx]
    
    plot_left_value <- approx(
      x = c(left_time, right_time),
      y = c(left_value, right_value),
      xout = plot_left_time,
      rule = 2
    )$y
    
    plot_right_value <- approx(
      x = c(left_time, right_time),
      y = c(left_value, right_value),
      xout = plot_right_time,
      rule = 2
    )$y
    
    data.frame(
      Event_ID = i,
      Peak_Index = pidx,
      Peak_Time = peak_time,
      Peak_Value = peak_value,
      Left_Min_Time = left_time,
      Left_Min_Value = left_value,
      Right_Min_Time = right_time,
      Right_Min_Value = right_value,
      Plot_Left_Time = plot_left_time,
      Plot_Left_Value = plot_left_value,
      Plot_Right_Time = plot_right_time,
      Plot_Right_Value = plot_right_value,
      Baseline_at_peak = baseline_at_peak,
      Baseline_FFT_Fraction = baseline_fft_fraction,
      stringsAsFactors = FALSE
    )
  })
  
  event_list <- event_list[!sapply(event_list, is.null)]
  if (length(event_list) == 0) return(empty)
  
  event_baselines <- do.call(rbind, event_list)
  
  # Build AUC/ribbon data only for the restricted positive component
  # around each peak.
  trace_list <- lapply(seq_len(nrow(event_baselines)), function(i) {
    
    idx <- which(
      time_vals >= event_baselines$Plot_Left_Time[i] &
        time_vals <= event_baselines$Plot_Right_Time[i]
    )
    
    if (length(idx) < 2) return(NULL)
    
    baseline_vec <- approx(
      x = c(event_baselines$Left_Min_Time[i], event_baselines$Right_Min_Time[i]),
      y = c(event_baselines$Left_Min_Value[i], event_baselines$Right_Min_Value[i]),
      xout = time_vals[idx],
      rule = 2
    )$y
    
    data.frame(
      Event_ID = event_baselines$Event_ID[i],
      Time = time_vals[idx],
      Baseline = baseline_vec,
      Signal = signal_vals[idx],
      Baseline_Smoothed_Signal = baseline_signal_vals[idx],
      AUC_Upper = pmax(signal_vals[idx], baseline_vec),
      stringsAsFactors = FALSE
    )
  })
  
  trace_list <- trace_list[!sapply(trace_list, is.null)]
  
  baseline_trace <- if (length(trace_list) > 0) {
    do.call(rbind, trace_list)
  } else {
    data.frame(Time = numeric(0), Baseline = numeric(0))
  }
  
  list(
    scalar = mean(event_baselines$Baseline_at_peak, na.rm = TRUE),
    trace = baseline_trace,
    event_baselines = event_baselines,
    baseline_smoothed = data.frame(Time = baseline_time_vals, signal = baseline_signal_vals),
    baseline_fft_fraction = baseline_fft_fraction,
    mode = "8"
  )
  
}

compute_trace_start_to_onset_mean_baseline <- function(data_smoothed,
                                                       time_start_increasin_peak = NULL) {
  # Constant baseline requested by the user:
  # take the first part of the calcium trace from the beginning of the trace
  # up to Time_Onset, calculate its mean, and use that single value as the
  # baseline for the complete trace and for all detected transients.
  time_vals <- suppressWarnings(as.numeric(data_smoothed[[1]]))
  signal_vals <- suppressWarnings(as.numeric(data_smoothed[[2]]))
  
  fallback_scalar <- if (length(signal_vals) > 0 && any(is.finite(signal_vals))) {
    mean(signal_vals, na.rm = TRUE)
  } else {
    0
  }
  
  empty_trace <- data.frame(Time = numeric(0), Baseline = numeric(0))
  
  if (length(time_vals) == 0 || length(signal_vals) == 0 ||
      all(!is.finite(time_vals)) || all(!is.finite(signal_vals))) {
    return(list(
      scalar = fallback_scalar,
      trace = empty_trace,
      event_baselines = NULL,
      mode = "11",
      Time_Onset_For_Baseline = NA_real_,
      Baseline_Window_Start = NA_real_,
      Baseline_Window_End = NA_real_,
      Baseline_Window_N = 0L
    ))
  }
  
  onset_times <- numeric(0)
  if (!is.null(time_start_increasin_peak) &&
      is.data.frame(time_start_increasin_peak) &&
      "Time" %in% names(time_start_increasin_peak)) {
    onset_times <- suppressWarnings(as.numeric(time_start_increasin_peak$Time))
    onset_times <- onset_times[is.finite(onset_times)]
  }
  
  # Time_Onset is the first transient onset in the trace.
  time_onset <- if (length(onset_times) > 0) {
    min(onset_times, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  if (!is.finite(time_onset)) {
    idx_baseline <- seq_len(max(1L, floor(0.10 * length(signal_vals))))
    time_onset <- max(time_vals[idx_baseline], na.rm = TRUE)
  } else {
    idx_baseline <- which(time_vals <= time_onset)
  }
  
  idx_baseline <- idx_baseline[
    idx_baseline >= 1 &
      idx_baseline <= length(signal_vals) &
      is.finite(signal_vals[idx_baseline]) &
      is.finite(time_vals[idx_baseline])
  ]
  
  if (length(idx_baseline) == 0) {
    idx_baseline <- which(is.finite(signal_vals) & is.finite(time_vals))
    if (length(idx_baseline) > 0) idx_baseline <- idx_baseline[1]
  }
  
  baseline_value <- if (length(idx_baseline) > 0) {
    mean(signal_vals[idx_baseline], na.rm = TRUE)
  } else {
    fallback_scalar
  }
  
  if (!is.finite(baseline_value)) baseline_value <- fallback_scalar
  
  # Keep a full constant trace so the plot can draw the selected baseline as
  # one horizontal line while baseline_at_times() still returns a constant value.
  baseline_trace <- data.frame(
    Time = time_vals,
    Baseline = rep(baseline_value, length(time_vals))
  )
  
  list(
    scalar = baseline_value,
    trace = baseline_trace,
    event_baselines = NULL,
    mode = "11",
    Time_Onset_For_Baseline = time_onset,
    Baseline_Window_Start = if (length(idx_baseline) > 0) min(time_vals[idx_baseline], na.rm = TRUE) else NA_real_,
    Baseline_Window_End = if (length(idx_baseline) > 0) max(time_vals[idx_baseline], na.rm = TRUE) else NA_real_,
    Baseline_Window_N = length(idx_baseline)
  )
}

# Backward-compatible alias. Previous drafts used this name for an event-specific
# mean-rise baseline. It now intentionally returns the corrected constant baseline.
compute_event_rise_mean_baseline <- function(data_smoothed,
                                             table_positions_peaks = NULL,
                                             time_start_increasin_peak = NULL) {
  compute_trace_start_to_onset_mean_baseline(
    data_smoothed = data_smoothed,
    time_start_increasin_peak = time_start_increasin_peak
  )
}

AUC_event_specific_minima <- function(data_smoothed, baseline_info) {
  if (is.null(baseline_info$trace) || nrow(baseline_info$trace) == 0) {
    return(list(area = NA_real_, with_absolute_error = NA_real_, P_min = NA_real_, P_max = NA_real_))
  }
  
  auc_by_event <- lapply(split(baseline_info$trace, baseline_info$trace$Event_ID), function(df) {
    if (nrow(df) < 2) return(NA_real_)
    positive_signal <- pmax(df$Signal - df$Baseline, 0)
    sum(diff(df$Time) * (head(positive_signal, -1) + tail(positive_signal, -1)) / 2, na.rm = TRUE)
  })
  
  area <- sum(unlist(auc_by_event), na.rm = TRUE)
  above_df <- baseline_info$trace[baseline_info$trace$Signal > baseline_info$trace$Baseline, , drop = FALSE]
  
  list(
    area = area,
    with_absolute_error = NA_real_,
    P_min = if (nrow(above_df) > 0) min(above_df$Time, na.rm = TRUE) else NA_real_,
    P_max = if (nrow(above_df) > 0) max(above_df$Time, na.rm = TRUE) else NA_real_
  )
}

AUC2_dynamic <- function(datos, baseline_trace) {
  colnames(datos) <- c("Time", "sing")
  baseline_interp <- approx(
    x = baseline_trace$Time,
    y = baseline_trace$Baseline,
    xout = datos$Time,
    rule = 2
  )$y
  positive_signal <- pmax(datos$sing - baseline_interp, 0)
  if (length(positive_signal) < 2) {
    return(list(area = NA_real_, with_absolute_error = NA_real_, P_min = NA_real_, P_max = NA_real_))
  }
  area <- sum(diff(datos$Time) * (head(positive_signal, -1) + tail(positive_signal, -1)) / 2, na.rm = TRUE)
  above <- which(positive_signal > 0)
  list(
    area = area,
    with_absolute_error = NA_real_,
    P_min = if (length(above) > 0) min(datos$Time[above], na.rm = TRUE) else NA_real_,
    P_max = if (length(above) > 0) max(datos$Time[above], na.rm = TRUE) else NA_real_
  )
}

# =========================
# Helper functions for FFT grid assessment
# =========================
build_fft_frequency_table <- function(X, n, k_keep) {
  n_pos <- floor(n / 2) + 1
  k <- 0:(n_pos - 1)
  omega <- 2 * pi * k / n
  freq_cyc <- k / n
  period <- ifelse(k == 0, Inf, n / k)
  X_pos <- X[1:n_pos]
  data.frame(
    k = k,
    omega = omega,
    freq = freq_cyc,
    period = period,
    X_re = Re(X_pos),
    X_im = Im(X_pos),
    X_k = sprintf("%.6f%+.6fi", Re(X_pos), Im(X_pos)),
    Mod = Mod(X_pos),
    Power = Mod(X_pos)^2,
    I = (Mod(X_pos)^2) / n,
    kept = k <= (k_keep - 1L)
  )
}

compute_baseline_value_from_mode <- function(
    data_smoothed,
    baseline_mode = 2,
    lim_inf = 0,
    lim_sup = 20,
    own_baseline = 0,
    initial_fraction = 0.10,
    minimal_fraction = 0.10,
    time_start_increasin_peak = NULL,
    rolling_window = 600,
    rolling_percentile = 2
) {
  baseline_mode <- as.character(baseline_mode)
  signal_vals <- data_smoothed$signal
  time_vals <- data_smoothed$Time
  
  switch(
    baseline_mode,
    "1" = 0,
    "2" = {
      frac <- suppressWarnings(as.numeric(initial_fraction))
      if (length(frac) == 0 || !is.finite(frac[1])) frac <- 0.10
      frac <- max(0.01, min(1, frac[1]))
      n_initial <- max(1L, floor(length(signal_vals) * frac))
      mean(signal_vals[seq_len(n_initial)], na.rm = TRUE)
    },
    "3" = {
      frac <- suppressWarnings(as.numeric(minimal_fraction))
      if (length(frac) == 0 || !is.finite(frac[1])) frac <- 0.10
      frac <- max(0.01, min(1, frac[1]))
      y_sorted <- sort(signal_vals[is.finite(signal_vals)])
      if (length(y_sorted) == 0) {
        mean(signal_vals, na.rm = TRUE)
      } else {
        n_low <- max(1L, floor(length(y_sorted) * frac))
        mean(y_sorted[seq_len(n_low)], na.rm = TRUE)
      }
    },
    "4" = own_baseline,
    "5" = min(signal_vals, na.rm = TRUE),
    "6" = {
      bt <- compute_rolling_percentile_baseline(
        data_smoothed = data_smoothed,
        window_size = rolling_window,
        percentile = rolling_percentile
      )
      stats::median(bt$Baseline, na.rm = TRUE)
    },
    "7" = {
      bt <- compute_windowed_percentile_baseline(
        data_smoothed = data_smoothed,
        window_size = rolling_window,
        percentile = rolling_percentile
      )
      stats::median(bt$Baseline, na.rm = TRUE)
    },
    "8" = mean(signal_vals, na.rm = TRUE),
    "9" = median(signal_vals, na.rm = TRUE),
    "10" = mean(signal_vals, na.rm = TRUE),
    "11" = {
      onset_times <- numeric(0)
      if (!is.null(time_start_increasin_peak) &&
          is.data.frame(time_start_increasin_peak) &&
          "Time" %in% names(time_start_increasin_peak)) {
        onset_times <- suppressWarnings(as.numeric(time_start_increasin_peak$Time))
        onset_times <- onset_times[is.finite(onset_times)]
      }
      if (length(onset_times) > 0) {
        time_onset <- min(onset_times, na.rm = TRUE)
        idx_onset <- which(time_vals <= time_onset)
        if (length(idx_onset) > 0) mean(signal_vals[idx_onset], na.rm = TRUE) else mean(signal_vals, na.rm = TRUE)
      } else {
        mean(signal_vals, na.rm = TRUE)
      }
    },
    mean(signal_vals, na.rm = TRUE)
  )
}

extract_peak_metrics_from_smoothed <- function(
    data_smoothed,
    nups,
    ndowns,
    minpeakheight,
    minpeakdistance,
    baseline_mode = 2,
    lim_inf = 0,
    lim_sup = 20,
    own_baseline = 0,
    initial_fraction = 0.10,
    minimal_fraction = 0.10,
    rolling_window = 600,
    rolling_percentile = 2,
    min_FWHP = 0,
    min_prominence = 0,
    baseline_info_override = NULL
) {
  empty_summary <- data.frame(
    n_peaks = 0,
    mean_peak = NA_real_,
    mean_peak_occurrence = NA_real_,
    mean_peak_rise_time = NA_real_,
    mean_prominence = NA_real_,
    mean_FWHM = NA_real_
  )
  empty_result <- list(
    df_p = data.frame(),
    baseline = NA_real_,
    summary = empty_summary,
    error = NULL
  )
  
  tryCatch({
    peaks_found <- peaks(
      data = data_smoothed,
      nups = nups,
      ndowns = ndowns,
      minpeakheight = minpeakheight,
      minpeakdistance = minpeakdistance
    )
    
    if (is.null(peaks_found$p_eak) || NROW(peaks_found$p_eak) == 0 ||
        is.null(peaks_found$peak) || NROW(peaks_found$peak) == 0) {
      if (!is.null(baseline_info_override) &&
          is.list(baseline_info_override) &&
          !is.null(baseline_info_override$scalar)) {
        empty_result$baseline <- suppressWarnings(as.numeric(baseline_info_override$scalar))[1]
      } else {
        empty_result$baseline <- compute_baseline_value_from_mode(
          data_smoothed = data_smoothed,
          baseline_mode = baseline_mode,
          lim_inf = lim_inf,
          lim_sup = lim_sup,
          own_baseline = own_baseline,
          initial_fraction = initial_fraction,
          minimal_fraction = minimal_fraction,
          time_start_increasin_peak = NULL,
          rolling_window = rolling_window,
          rolling_percentile = rolling_percentile
        )
      }
      return(empty_result)
    }
    
    table_peak <- normalize_fft_peak_table(peaks_found$p_eak)
    table_positions_peaks <- normalize_fft_peak_positions(peaks_found$peak)
    if (nrow(table_positions_peaks) == 0 || nrow(table_peak) == 0) return(empty_result)
    peaks_idx <- table_positions_peaks$peak_index
    
    MSCPFP <- Time_of_the_first_peak(
      data1 = data_smoothed,
      peak = table_positions_peaks
    )$cambios_menor_que_pfp
    
    prom <- prominens2(
      data = data_smoothed,
      peak = table_positions_peaks,
      MSCPFP = MSCPFP
    )
    
    df_peaks_parcia <- prom$df_peaks_parcia
    time_start_increasin_peak <- prom$time_start_increasin_peak
    
    Puntos_medios <- FWHP2(
      peaks = data_smoothed[, 1][peaks_idx],
      df_peaks_parcia = df_peaks_parcia
    )$Puntos_medios
    
    table_peak$prominence <- prom$prominens_amplitud
    table_peak$Prominence_Midpoint <- Puntos_medios$p_eak_mediun
    
    right_left <- right_left_FWHP(
      data1 = data_smoothed,
      peak = table_positions_peaks,
      P_M = Puntos_medios
    )
    
    table_peak$Time_left_FWHP <- right_left$df$Time_left_FWHP
    table_peak$Time_right_FWHP <- right_left$df2$Time_right_FWHP
    table_peak$FWHP <- right_left$df2$Time_right_FWHP - right_left$df$Time_left_FWHP
    table_peak$Time_to_peak <- table_peak$posision_peaks - time_start_increasin_peak$Time
    table_peak$puntominimo_y <- prom$df_peaks_parcia$p_fin1
    table_peak$Transient_Ocurrence_Time <- time_start_increasin_peak$Time
    
    deri1_vals <- prospectr::savitzkyGolay(
      X = data_smoothed$signal,
      m = 1,
      p = 2,
      w = 5
    )
    
    n_der <- length(deri1_vals)
    n_time <- nrow(data_smoothed)
    start_idx <- floor((n_time - n_der) / 2) + 1
    end_idx <- start_idx + n_der - 1
    
    primera_derivada1 <- data.frame(
      Time = data_smoothed$Time[start_idx:end_idx],
      deri1 = deri1_vals
    )
    
    data_min <- prom$data_min
    data_minimos_crecientes <- data.frame(
      x1 = time_start_increasin_peak$Time,
      y1 = data_min$y,
      x2 = table_peak$posision_peaks,
      y2 = table_peak$absolute_amplitude
    )
    
    slope <- numeric(nrow(data_minimos_crecientes))
    for (i in seq_len(nrow(data_minimos_crecientes))) {
      resultados_filtrados <- primera_derivada1[
        primera_derivada1$Time >= data_minimos_crecientes$x1[i] &
          primera_derivada1$Time <= data_minimos_crecientes$x2[i],
        ,
        drop = FALSE
      ]
      slope[i] <- if (nrow(resultados_filtrados) > 0 && any(is.finite(resultados_filtrados$deri1))) {
        max(resultados_filtrados$deri1, na.rm = TRUE)
      } else {
        NA_real_
      }
    }
    table_peak$slope <- slope
    
    if (!is.null(baseline_info_override) &&
        is.list(baseline_info_override) &&
        !is.null(baseline_info_override$scalar)) {
      baseline1 <- suppressWarnings(as.numeric(baseline_info_override$scalar))[1]
    } else {
      baseline1 <- compute_baseline_value_from_mode(
        data_smoothed = data_smoothed,
        baseline_mode = baseline_mode,
        lim_inf = lim_inf,
        lim_sup = lim_sup,
        own_baseline = own_baseline,
        initial_fraction = initial_fraction,
        minimal_fraction = minimal_fraction,
        time_start_increasin_peak = time_start_increasin_peak,
        rolling_window = rolling_window,
        rolling_percentile = rolling_percentile
      )
    }
    
    baseline_for_fwhm <- rep(baseline1, length(peaks_idx))
    if (!is.null(baseline_info_override) && is.list(baseline_info_override)) {
      baseline_for_fwhm <- baseline_at_times(
        time_points = data_smoothed[, 1][peaks_idx],
        baseline_info = baseline_info_override
      )
    }
    
    p_eak_mediun <- c((table_positions_peaks[, 1] + baseline_for_fwhm) / 2)
    Puntos_medios_FWHM <- data.frame(
      posiscion_medio = data_smoothed[, 1][peaks_idx],
      p_eak_mediun = p_eak_mediun
    )
    
    right_left_FWHM <- right_left_FWHP(
      data1 = data_smoothed,
      peak = table_positions_peaks,
      P_M = Puntos_medios_FWHM
    )
    
    df_p <- table_peak
    df_p$FWHM <- right_left_FWHM$df2$Time_right_FWHP - right_left_FWHM$df$Time_left_FWHP
    colnames(df_p) <- c(
      "Amplitude", "Peak_Occurence_Time", "L_inf", "L_sup",
      "Prominence", "Prominence_Midpoint", "Time_left_FWHP",
      "Time_right_FWHP", "FWHP", "Peak_Rise_Time",
      "puntominimo_y", "Transient_Ocurrence_Time",
      "Rise_Rate", "FWHM"
    )
    
    df_FWHM1 <- data.frame(
      Time_left_FWHM = right_left_FWHM$df$Time_left_FWHP,
      Time_right_FWHM = right_left_FWHM$df2$Time_right_FWHP,
      Amplitude_Midpoint = p_eak_mediun
    )
    
    df_p <- cbind(df_p, df_FWHM1)
    df_p <- df_p[df_p$FWHP > min_FWHP, , drop = FALSE]
    df_p <- df_p[df_p$Prominence > min_prominence, , drop = FALSE]
    
    if (nrow(df_p) > 0) {
      if (!is.null(baseline_info_override) && is.list(baseline_info_override)) {
        df_p$Baseline_at_peak <- baseline_at_times(
          time_points = df_p$Peak_Occurence_Time,
          baseline_info = baseline_info_override
        )
      } else {
        df_p$Baseline_at_peak <- rep(baseline1, nrow(df_p))
      }
      
      valid_baseline_filter <- is.finite(df_p$Baseline_at_peak)
      if (any(valid_baseline_filter)) {
        df_p <- df_p[valid_baseline_filter & df_p$Amplitude > df_p$Baseline_at_peak, , drop = FALSE]
      }
    }
    
    if (nrow(df_p) > 0) {
      df_p <- df_p[order(df_p$Peak_Occurence_Time), , drop = FALSE]
      mean_peak_vals <- df_p$Amplitude
      if ("Baseline_at_peak" %in% names(df_p) && any(is.finite(df_p$Baseline_at_peak))) {
        mean_peak_vals <- mean_peak_vals - df_p$Baseline_at_peak
      } else if (is.finite(baseline1)) {
        mean_peak_vals <- mean_peak_vals - baseline1
      }
      summary_df <- data.frame(
        n_peaks = nrow(table_peak),
        mean_peak = mean(table_peak$absolute_amplitude, na.rm = TRUE),
        mean_peak_occurrence = mean(table_peak$posision_peaks, na.rm = TRUE),
        mean_peak_rise_time = mean(df_p$Peak_Rise_Time, na.rm = TRUE),
        mean_prominence = mean(df_p$Prominence, na.rm = TRUE),
        mean_FWHM = mean(df_p$FWHM, na.rm = TRUE)
      )
    } else {
      summary_df <- empty_summary
    }
    
    list(
      df_p = df_p,
      baseline = baseline1,
      summary = summary_df,
      error = NULL
    )
  }, error = function(e) {
    empty_result$error <- conditionMessage(e)
    empty_result
  })
}

compare_peak_metrics_to_reference <- function(candidate_df, reference_df, dt) {
  out <- list(
    same_n_peaks = FALSE,
    peak_occurrence_shift = NA_real_,
    rise_time_rel_change = NA_real_,
    prominence_rel_change = NA_real_,
    fwhm_rel_change = NA_real_,
    pass_metrics = FALSE
  )
  
  if (nrow(candidate_df) == 0 && nrow(reference_df) == 0) {
    out$same_n_peaks <- TRUE
    out$pass_metrics <- TRUE
    return(out)
  }
  
  if (nrow(candidate_df) == 0 || nrow(reference_df) == 0) return(out)
  
  candidate_df <- candidate_df[order(candidate_df$Peak_Occurence_Time), , drop = FALSE]
  reference_df <- reference_df[order(reference_df$Peak_Occurence_Time), , drop = FALSE]
  m <- min(nrow(candidate_df), nrow(reference_df))
  if (m <= 0) return(out)
  
  eps <- 1e-8
  out$same_n_peaks <- nrow(candidate_df) == nrow(reference_df)
  out$peak_occurrence_shift <- mean(
    abs(candidate_df$Peak_Occurence_Time[1:m] - reference_df$Peak_Occurence_Time[1:m]),
    na.rm = TRUE
  )
  out$rise_time_rel_change <- median(
    abs(candidate_df$Peak_Rise_Time[1:m] - reference_df$Peak_Rise_Time[1:m]) /
      (abs(reference_df$Peak_Rise_Time[1:m]) + eps),
    na.rm = TRUE
  )
  out$prominence_rel_change <- median(
    abs(candidate_df$Prominence[1:m] - reference_df$Prominence[1:m]) /
      (abs(reference_df$Prominence[1:m]) + eps),
    na.rm = TRUE
  )
  out$fwhm_rel_change <- median(
    abs(candidate_df$FWHM[1:m] - reference_df$FWHM[1:m]) /
      (abs(reference_df$FWHM[1:m]) + eps),
    na.rm = TRUE
  )
  
  out$pass_metrics <- isTRUE(out$same_n_peaks) &&
    is.finite(out$peak_occurrence_shift) && out$peak_occurrence_shift <= (2 * dt) &&
    is.finite(out$rise_time_rel_change) && out$rise_time_rel_change <= 0.20 &&
    is.finite(out$prominence_rel_change) && out$prominence_rel_change <= 0.20 &&
    is.finite(out$fwhm_rel_change) && out$fwhm_rel_change <= 0.20
  
  out
}



# ============================================================
# Optional sigmoid modeling for FFT module
# ============================================================
# These functions are intentionally optional. They do not modify the
# original FFT smoothing, peak detection, baseline, AUC, FWHM, or FWHP
# workflow. They only add extra columns when the user activates the
# corresponding checkbox in the UI.
#
# Single sigmoid: onset and rise-time refinement.
# Double sigmoid: full-transient model, including offset, duration,
# rise slope, decay slope, and model-based metrics.
# ============================================================

fft_safe_num <- function(x, fallback = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || !is.finite(x[1])) return(fallback)
  x[1]
}

fft_sigmoid_simple_fun <- function(t, Baseline, A, k, t0) {
  Baseline + A / (1 + exp(-k * (t - t0)))
}

fft_double_sigmoid_fun <- function(t, Baseline, A, k1, t1, k2, t2) {
  s1 <- 1 / (1 + exp(-k1 * (t - t1)))
  s2 <- 1 / (1 + exp(-k2 * (t - t2)))
  Baseline + A * (s1 - s2)
}

fft_safe_r2 <- function(y, yhat) {
  y <- suppressWarnings(as.numeric(y))
  yhat <- suppressWarnings(as.numeric(yhat))
  ok <- is.finite(y) & is.finite(yhat)
  if (sum(ok) < 3) return(NA_real_)
  ss_res <- sum((y[ok] - yhat[ok])^2, na.rm = TRUE)
  ss_tot <- sum((y[ok] - mean(y[ok], na.rm = TRUE))^2, na.rm = TRUE)
  if (!is.finite(ss_tot) || ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}

fft_compute_event_baseline_rms <- function(time,
                                           signal,
                                           reference_time,
                                           pre_event_window = 5,
                                           fallback_baseline = NA_real_) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  reference_time <- fft_safe_num(reference_time, NA_real_)
  pre_event_window <- fft_safe_num(pre_event_window, 5)
  fallback_baseline <- fft_safe_num(fallback_baseline, NA_real_)
  
  if (!is.finite(pre_event_window) || pre_event_window <= 0) pre_event_window <- 5
  
  idx <- which(time >= reference_time - pre_event_window & time < reference_time)
  if (length(idx) < 3) {
    idx <- which(time < reference_time)
    if (length(idx) > 10) idx <- tail(idx, 10)
  }
  
  if (length(idx) >= 2) {
    Baseline <- stats::median(signal[idx], na.rm = TRUE)
    rms <- sqrt(mean((signal[idx] - Baseline)^2, na.rm = TRUE))
  } else {
    Baseline <- fallback_baseline
    rms <- stats::sd(signal, na.rm = TRUE)
  }
  
  if (!is.finite(Baseline)) Baseline <- stats::median(signal, na.rm = TRUE)
  if (!is.finite(rms) || rms <= 0) rms <- stats::sd(signal, na.rm = TRUE)
  if (!is.finite(rms) || rms <= 0) rms <- 0
  
  list(
    Baseline = Baseline,
    RMS = rms,
    Threshold = Baseline + rms,
    Pre_Event_N = length(idx),
    Pre_Event_Start = ifelse(length(idx) > 0, min(time[idx], na.rm = TRUE), NA_real_),
    Pre_Event_End = ifelse(length(idx) > 0, max(time[idx], na.rm = TRUE), NA_real_)
  )
}

fit_fft_single_sigmoid_event <- function(time,
                                         signal,
                                         left_reference_time,
                                         peak_time,
                                         peak_value,
                                         fallback_baseline = NA_real_,
                                         pre_event_window = 5,
                                         fit_left_pad = 0,
                                         fit_right_pad = 0,
                                         dense_n = 300) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  left_reference_time <- fft_safe_num(left_reference_time, NA_real_)
  peak_time <- fft_safe_num(peak_time, NA_real_)
  peak_value <- fft_safe_num(peak_value, NA_real_)
  fallback_baseline <- fft_safe_num(fallback_baseline, NA_real_)
  
  empty <- list(
    SingleSigmoid_Onset_Time = NA_real_,
    SingleSigmoid_Rise_Time = NA_real_,
    SingleSigmoid_Rise_Rate = NA_real_,
    SingleSigmoid_Baseline = NA_real_,
    SingleSigmoid_RMS = NA_real_,
    SingleSigmoid_Threshold = NA_real_,
    SingleSigmoid_A = NA_real_,
    SingleSigmoid_k = NA_real_,
    SingleSigmoid_t0 = NA_real_,
    SingleSigmoid_Slope_at_t0 = NA_real_,
    SingleSigmoid_R2 = NA_real_,
    SingleSigmoid_Fit_Status = "Not fitted",
    fit_df = data.frame()
  )
  
  if (!is.finite(left_reference_time) || !is.finite(peak_time) || peak_time <= left_reference_time) {
    empty$SingleSigmoid_Fit_Status <- "Failed: invalid left reference/peak times"
    return(empty)
  }
  
  baseline_info <- fft_compute_event_baseline_rms(
    time = time,
    signal = signal,
    reference_time = left_reference_time,
    pre_event_window = pre_event_window,
    fallback_baseline = fallback_baseline
  )
  
  Baseline_fixed <- baseline_info$Baseline
  threshold <- baseline_info$Threshold
  fit_left_pad <- fft_safe_num(fit_left_pad, 0)
  fit_right_pad <- fft_safe_num(fit_right_pad, 0)
  if (!is.finite(fit_left_pad) || fit_left_pad < 0) fit_left_pad <- 0
  if (!is.finite(fit_right_pad) || fit_right_pad < 0) fit_right_pad <- 0
  
  fit_start <- max(min(time, na.rm = TRUE), left_reference_time - fit_left_pad)
  fit_end <- min(max(time, na.rm = TRUE), peak_time + fit_right_pad)
  idx <- which(time >= fit_start & time <= fit_end)
  
  if (length(idx) < 5) {
    empty$SingleSigmoid_Baseline <- Baseline_fixed
    empty$SingleSigmoid_RMS <- baseline_info$RMS
    empty$SingleSigmoid_Threshold <- threshold
    empty$SingleSigmoid_Fit_Status <- "Failed: too few points"
    return(empty)
  }
  
  fit_data <- data.frame(t = time[idx], y = signal[idx])
  fit_data <- fit_data[is.finite(fit_data$t) & is.finite(fit_data$y), , drop = FALSE]
  if (nrow(fit_data) < 5) {
    empty$SingleSigmoid_Baseline <- Baseline_fixed
    empty$SingleSigmoid_RMS <- baseline_info$RMS
    empty$SingleSigmoid_Threshold <- threshold
    empty$SingleSigmoid_Fit_Status <- "Failed: insufficient finite points"
    return(empty)
  }
  
  if (!is.finite(peak_value)) {
    peak_value <- signal[which.min(abs(time - peak_time))]
  }
  
  A_start <- max(peak_value - Baseline_fixed, diff(range(fit_data$y, na.rm = TRUE)), 1e-6)
  rise_span <- peak_time - left_reference_time
  if (!is.finite(rise_span) || rise_span <= 0) rise_span <- stats::median(diff(time), na.rm = TRUE)
  if (!is.finite(rise_span) || rise_span <= 0) rise_span <- 1
  
  k_start <- 4 / rise_span
  t0_start <- left_reference_time + 0.5 * (peak_time - left_reference_time)
  
  fit_obj <- tryCatch(
    stats::nls(
      y ~ Baseline_fixed + A / (1 + exp(-k * (t - t0))),
      data = fit_data,
      start = list(A = A_start, k = k_start, t0 = t0_start),
      algorithm = "port",
      lower = c(A = 0, k = 1e-6, t0 = min(fit_data$t, na.rm = TRUE)),
      upper = c(A = 10 * max(abs(A_start), 1e-6), k = 1e4, t0 = peak_time),
      control = stats::nls.control(maxiter = 300, warnOnly = TRUE)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_obj)) {
    empty$SingleSigmoid_Baseline <- Baseline_fixed
    empty$SingleSigmoid_RMS <- baseline_info$RMS
    empty$SingleSigmoid_Threshold <- threshold
    empty$SingleSigmoid_Fit_Status <- "Failed: nls did not converge"
    return(empty)
  }
  
  pars <- stats::coef(fit_obj)
  yhat_train <- Baseline_fixed + pars[["A"]] / (1 + exp(-pars[["k"]] * (fit_data$t - pars[["t0"]])))
  r2 <- fft_safe_r2(fit_data$y, yhat_train)
  
  dense_t <- seq(fit_start, fit_end, length.out = dense_n)
  dense_y <- Baseline_fixed + pars[["A"]] / (1 + exp(-pars[["k"]] * (dense_t - pars[["t0"]])))
  onset_idx <- which(dense_t <= peak_time & dense_y >= threshold)
  onset_time <- if (length(onset_idx) > 0) dense_t[min(onset_idx)] else NA_real_
  rise_time <- ifelse(is.finite(onset_time), peak_time - onset_time, NA_real_)
  rise_rate <- ifelse(is.finite(rise_time) && rise_time > 0, (peak_value - Baseline_fixed) / rise_time, NA_real_)
  slope_t0 <- pars[["A"]] * pars[["k"]] / 4
  
  list(
    SingleSigmoid_Onset_Time = onset_time,
    SingleSigmoid_Rise_Time = rise_time,
    SingleSigmoid_Rise_Rate = rise_rate,
    SingleSigmoid_Baseline = Baseline_fixed,
    SingleSigmoid_RMS = baseline_info$RMS,
    SingleSigmoid_Threshold = threshold,
    SingleSigmoid_A = pars[["A"]],
    SingleSigmoid_k = pars[["k"]],
    SingleSigmoid_t0 = pars[["t0"]],
    SingleSigmoid_Slope_at_t0 = slope_t0,
    SingleSigmoid_R2 = r2,
    SingleSigmoid_Fit_Status = "OK",
    fit_df = data.frame(
      Event_ID = NA_integer_,
      Time = dense_t,
      Fitted = dense_y,
      Model = "Single sigmoid",
      stringsAsFactors = FALSE
    )
  )
}


# ------------------------------------------------------------
# Automatic decay-end finder for the double-sigmoid fit
# ------------------------------------------------------------
# The previous double-sigmoid fit often used Time_right_FWHM/FWHP as the
# right reference when an event-specific right minimum was not available.
# That makes the fitted purple curve stop around half-height instead of
# following the transient down to the end of the visible calcium activity.
# This helper searches after each detected peak until the smoothed trace
# returns close to the local baseline and then uses the first local minimum
# after that return as the right reference for the double-sigmoid fit.
# The search is clipped by the midpoint to the next detected peak so that
# one transient does not borrow the decay interval of the following event.
fft_find_double_sigmoid_decay_end <- function(time,
                                              signal,
                                              peak_time,
                                              left_reference_time = NA_real_,
                                              fallback_baseline = NA_real_,
                                              pre_event_window = 5,
                                              next_peak_time = NA_real_,
                                              min_right_time = NA_real_) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  peak_time <- fft_safe_num(peak_time, NA_real_)
  left_reference_time <- fft_safe_num(left_reference_time, NA_real_)
  fallback_baseline <- fft_safe_num(fallback_baseline, NA_real_)
  pre_event_window <- fft_safe_num(pre_event_window, 5)
  next_peak_time <- fft_safe_num(next_peak_time, NA_real_)
  min_right_time <- fft_safe_num(min_right_time, NA_real_)
  
  if (length(time) != length(signal) || length(time) < 4 || !is.finite(peak_time)) {
    return(NA_real_)
  }
  if (all(!is.finite(time)) || all(!is.finite(signal))) {
    return(NA_real_)
  }
  
  # Fill missing signal values by linear interpolation so local-minimum search
  # is not broken by isolated NA/Inf values.
  bad_signal <- which(!is.finite(signal))
  if (length(bad_signal) > 0) {
    good_signal <- which(is.finite(signal))
    if (length(good_signal) >= 2) {
      signal[bad_signal] <- stats::approx(
        x = time[good_signal],
        y = signal[good_signal],
        xout = time[bad_signal],
        rule = 2
      )$y
    } else {
      return(NA_real_)
    }
  }
  
  pidx <- which.min(abs(time - peak_time))
  if (length(pidx) == 0 || !is.finite(pidx) || pidx >= length(time)) {
    return(NA_real_)
  }
  
  time_min <- min(time, na.rm = TRUE)
  time_max <- max(time, na.rm = TRUE)
  if (!is.finite(time_min) || !is.finite(time_max) || time_max <= time_min) {
    return(NA_real_)
  }
  
  # Do not let the decay of one transient cross into the next transient.
  right_bound_time <- if (is.finite(next_peak_time) && next_peak_time > peak_time) {
    min(time_max, (peak_time + next_peak_time) / 2)
  } else {
    time_max
  }
  if (!is.finite(right_bound_time) || right_bound_time <= peak_time) {
    right_bound_time <- time_max
  }
  
  right_bound_idx <- max(which(time <= right_bound_time))
  if (length(right_bound_idx) == 0 || !is.finite(right_bound_idx) || right_bound_idx <= pidx) {
    right_bound_idx <- length(time)
  }
  
  if (!is.finite(min_right_time) || min_right_time <= peak_time) {
    min_right_time <- peak_time
  }
  min_right_time <- min(min_right_time, time[right_bound_idx])
  
  baseline_reference_time <- if (is.finite(left_reference_time)) {
    left_reference_time
  } else {
    peak_time
  }
  
  local_baseline <- fft_compute_event_baseline_rms(
    time = time,
    signal = signal,
    reference_time = baseline_reference_time,
    pre_event_window = pre_event_window,
    fallback_baseline = fallback_baseline
  )
  
  baseline_level <- local_baseline$Baseline
  rms_level <- local_baseline$RMS
  peak_value <- signal[pidx]
  
  if (!is.finite(baseline_level)) {
    baseline_level <- stats::median(signal, na.rm = TRUE)
  }
  if (!is.finite(rms_level) || rms_level < 0) {
    rms_level <- 0
  }
  
  amplitude_above_baseline <- peak_value - baseline_level
  if (!is.finite(amplitude_above_baseline)) {
    amplitude_above_baseline <- diff(range(signal, na.rm = TRUE))
  }
  if (!is.finite(amplitude_above_baseline) || amplitude_above_baseline < 0) {
    amplitude_above_baseline <- 0
  }
  
  # A small amplitude fraction is used as a floor so that noiseless or very low
  # RMS traces still have a meaningful return-to-baseline criterion.
  end_threshold <- baseline_level + max(rms_level, 0.05 * amplitude_above_baseline, 1e-8)
  
  search_idx <- seq.int(pidx + 1L, right_bound_idx)
  search_idx <- search_idx[time[search_idx] >= min_right_time]
  if (length(search_idx) == 0) {
    return(time[right_bound_idx])
  }
  
  local_minima_in <- function(idx_vec) {
    idx_vec <- idx_vec[is.finite(idx_vec)]
    idx_vec <- idx_vec[idx_vec > 1 & idx_vec < length(signal)]
    if (length(idx_vec) == 0) return(integer(0))
    idx_vec[
      signal[idx_vec] <= signal[idx_vec - 1L] &
        signal[idx_vec] <= signal[idx_vec + 1L]
    ]
  }
  
  below_threshold <- search_idx[signal[search_idx] <= end_threshold]
  candidate_idx <- NA_integer_
  
  if (length(below_threshold) > 0) {
    # Prefer the first local minimum after the trace has returned close to the
    # local baseline. This visually follows the full decay and avoids stopping
    # at FWHM/FWHP.
    after_return <- search_idx[search_idx >= min(below_threshold)]
    local_after_return <- local_minima_in(after_return)
    if (length(local_after_return) > 0) {
      candidate_idx <- local_after_return[1]
    } else {
      candidate_idx <- below_threshold[1]
    }
  } else {
    # If the trace never reaches the baseline + RMS threshold before the next
    # peak/boundary, use the deepest point in the allowed decay interval.
    local_candidates <- local_minima_in(search_idx)
    if (length(local_candidates) > 0) {
      candidate_idx <- local_candidates[which.min(signal[local_candidates])]
    } else {
      candidate_idx <- search_idx[which.min(signal[search_idx])]
    }
  }
  
  if (length(candidate_idx) == 0 || !is.finite(candidate_idx)) {
    return(time[right_bound_idx])
  }
  
  candidate_time <- time[candidate_idx]
  if (!is.finite(candidate_time) || candidate_time <= peak_time) {
    candidate_time <- time[right_bound_idx]
  }
  
  candidate_time
}

fit_fft_double_sigmoid_event <- function(time,
                                         signal,
                                         left_reference_time,
                                         peak_time,
                                         right_reference_time,
                                         peak_value,
                                         fallback_baseline = NA_real_,
                                         pre_event_window = 5,
                                         fit_left_pad = 0,
                                         fit_right_pad = 0,
                                         dense_n = 400) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  left_reference_time <- fft_safe_num(left_reference_time, NA_real_)
  peak_time <- fft_safe_num(peak_time, NA_real_)
  right_reference_time <- fft_safe_num(right_reference_time, NA_real_)
  peak_value <- fft_safe_num(peak_value, NA_real_)
  fallback_baseline <- fft_safe_num(fallback_baseline, NA_real_)
  
  empty <- list(
    DoubleSigmoid_Onset_Time = NA_real_,
    DoubleSigmoid_Offset_Time = NA_real_,
    DoubleSigmoid_Duration = NA_real_,
    DoubleSigmoid_Rise_Time = NA_real_,
    DoubleSigmoid_Decay_Time = NA_real_,
    DoubleSigmoid_Baseline = NA_real_,
    DoubleSigmoid_RMS = NA_real_,
    DoubleSigmoid_Threshold = NA_real_,
    DoubleSigmoid_A = NA_real_,
    DoubleSigmoid_k1 = NA_real_,
    DoubleSigmoid_t1 = NA_real_,
    DoubleSigmoid_k2 = NA_real_,
    DoubleSigmoid_t2 = NA_real_,
    DoubleSigmoid_Rise_Slope = NA_real_,
    DoubleSigmoid_Decay_Slope = NA_real_,
    DoubleSigmoid_Model_Peak_Time = NA_real_,
    DoubleSigmoid_Model_Peak_Value = NA_real_,
    DoubleSigmoid_Left_Min_Time = NA_real_,
    DoubleSigmoid_Left_Min_Value = NA_real_,
    DoubleSigmoid_Right_Min_Time = NA_real_,
    DoubleSigmoid_Right_Min_Value = NA_real_,
    DoubleSigmoid_Baseline_at_peak = NA_real_,
    DoubleSigmoid_Prominence = NA_real_,
    DoubleSigmoid_R2 = NA_real_,
    DoubleSigmoid_Fit_Status = "Not fitted",
    fit_df = data.frame()
  )
  
  if (!is.finite(left_reference_time) || !is.finite(peak_time) || !is.finite(right_reference_time) ||
      peak_time <= left_reference_time || right_reference_time <= peak_time) {
    empty$DoubleSigmoid_Fit_Status <- "Failed: invalid left/peak/right times"
    return(empty)
  }
  
  baseline_info <- fft_compute_event_baseline_rms(
    time = time,
    signal = signal,
    reference_time = left_reference_time,
    pre_event_window = pre_event_window,
    fallback_baseline = fallback_baseline
  )
  
  Baseline_fixed <- baseline_info$Baseline
  threshold <- baseline_info$Threshold
  fit_left_pad <- fft_safe_num(fit_left_pad, 0)
  fit_right_pad <- fft_safe_num(fit_right_pad, 0)
  if (!is.finite(fit_left_pad) || fit_left_pad < 0) fit_left_pad <- 0
  if (!is.finite(fit_right_pad) || fit_right_pad < 0) fit_right_pad <- 0
  
  fit_start <- max(min(time, na.rm = TRUE), left_reference_time - fit_left_pad)
  fit_end <- min(max(time, na.rm = TRUE), right_reference_time + fit_right_pad)
  idx <- which(time >= fit_start & time <= fit_end)
  
  if (length(idx) < 8) {
    empty$DoubleSigmoid_Baseline <- Baseline_fixed
    empty$DoubleSigmoid_RMS <- baseline_info$RMS
    empty$DoubleSigmoid_Threshold <- threshold
    empty$DoubleSigmoid_Fit_Status <- "Failed: too few points"
    return(empty)
  }
  
  fit_data <- data.frame(t = time[idx], y = signal[idx])
  fit_data <- fit_data[is.finite(fit_data$t) & is.finite(fit_data$y), , drop = FALSE]
  if (nrow(fit_data) < 8) {
    empty$DoubleSigmoid_Baseline <- Baseline_fixed
    empty$DoubleSigmoid_RMS <- baseline_info$RMS
    empty$DoubleSigmoid_Threshold <- threshold
    empty$DoubleSigmoid_Fit_Status <- "Failed: insufficient finite points"
    return(empty)
  }
  
  if (!is.finite(peak_value)) peak_value <- signal[which.min(abs(time - peak_time))]
  A_start <- max(peak_value - Baseline_fixed, diff(range(fit_data$y, na.rm = TRUE)), 1e-6)
  rise_span <- peak_time - left_reference_time
  decay_span <- right_reference_time - peak_time
  if (!is.finite(rise_span) || rise_span <= 0) rise_span <- stats::median(diff(time), na.rm = TRUE)
  if (!is.finite(decay_span) || decay_span <= 0) decay_span <- stats::median(diff(time), na.rm = TRUE)
  if (!is.finite(rise_span) || rise_span <= 0) rise_span <- 1
  if (!is.finite(decay_span) || decay_span <= 0) decay_span <- 1
  
  k1_start <- 4 / rise_span
  k2_start <- 4 / decay_span
  t1_start <- left_reference_time + 0.5 * (peak_time - left_reference_time)
  t2_start <- peak_time + 0.5 * (right_reference_time - peak_time)
  
  fit_obj <- tryCatch(
    stats::nls(
      y ~ Baseline_fixed + A * (
        1 / (1 + exp(-k1 * (t - t1))) -
          1 / (1 + exp(-k2 * (t - t2)))
      ),
      data = fit_data,
      start = list(A = A_start, k1 = k1_start, t1 = t1_start, k2 = k2_start, t2 = t2_start),
      algorithm = "port",
      lower = c(A = 0, k1 = 1e-6, t1 = min(fit_data$t, na.rm = TRUE), k2 = 1e-6, t2 = peak_time),
      upper = c(A = 10 * max(abs(A_start), 1e-6), k1 = 1e4, t1 = peak_time, k2 = 1e4, t2 = max(fit_data$t, na.rm = TRUE)),
      control = stats::nls.control(maxiter = 400, warnOnly = TRUE)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_obj)) {
    empty$DoubleSigmoid_Baseline <- Baseline_fixed
    empty$DoubleSigmoid_RMS <- baseline_info$RMS
    empty$DoubleSigmoid_Threshold <- threshold
    empty$DoubleSigmoid_Fit_Status <- "Failed: nls did not converge"
    return(empty)
  }
  
  pars <- stats::coef(fit_obj)
  yhat_train <- Baseline_fixed + pars[["A"]] * (
    1 / (1 + exp(-pars[["k1"]] * (fit_data$t - pars[["t1"]]))) -
      1 / (1 + exp(-pars[["k2"]] * (fit_data$t - pars[["t2"]])))
  )
  r2 <- fft_safe_r2(fit_data$y, yhat_train)
  
  dense_t <- seq(fit_start, fit_end, length.out = dense_n)
  dense_y <- fft_double_sigmoid_fun(dense_t, Baseline_fixed, pars[["A"]], pars[["k1"]], pars[["t1"]], pars[["k2"]], pars[["t2"]])
  model_peak_idx <- which.max(dense_y)
  model_peak_time <- dense_t[model_peak_idx]
  model_peak_value <- dense_y[model_peak_idx]
  
  # User-requested double-sigmoid minima.
  # These values are used only when the user activates the double-sigmoid option.
  # Left minimum: minimum of the fitted double-sigmoid curve on the left side of the detected peak.
  # Right minimum: minimum of the fitted double-sigmoid curve on the right side of the detected peak.
  left_min_time <- NA_real_
  left_min_value <- NA_real_
  right_min_time <- NA_real_
  right_min_value <- NA_real_
  ds_baseline_at_peak <- NA_real_
  ds_prominence <- NA_real_
  
  left_side_idx <- which(dense_t <= peak_time)
  if (length(left_side_idx) > 0) {
    j_left <- left_side_idx[which.min(dense_y[left_side_idx])]
    left_min_time <- dense_t[j_left]
    left_min_value <- dense_y[j_left]
  }
  
  right_side_idx <- which(dense_t >= peak_time)
  if (length(right_side_idx) > 0) {
    j_right <- right_side_idx[which.min(dense_y[right_side_idx])]
    right_min_time <- dense_t[j_right]
    right_min_value <- dense_y[j_right]
  }
  
  if (is.finite(left_min_time) && is.finite(right_min_time) &&
      is.finite(left_min_value) && is.finite(right_min_value) &&
      right_min_time > left_min_time) {
    ds_baseline_at_peak <- stats::approx(
      x = c(left_min_time, right_min_time),
      y = c(left_min_value, right_min_value),
      xout = peak_time,
      rule = 2
    )$y
  }
  
  if (is.finite(peak_value) && is.finite(left_min_value)) {
    ds_prominence <- peak_value - left_min_value
  }
  
  onset_idx <- which(dense_t <= model_peak_time & dense_y >= threshold)
  offset_idx <- which(dense_t >= model_peak_time & dense_y <= threshold)
  
  ds_onset <- if (length(onset_idx) > 0) {
    dense_t[min(onset_idx)]
  } else {
    NA_real_
  }
  
  ds_offset <- if (length(offset_idx) > 0) {
    dense_t[min(offset_idx)]
  } else {
    NA_real_
  }
  
  # ------------------------------------------------------------
  # Robust double-sigmoid decay metrics
  # ------------------------------------------------------------
  # Primary definition:
  #   DoubleSigmoid_Decay_Time = DoubleSigmoid_Offset_Time - detected peak time
  #
  # If the fitted curve does not cross the Baseline + RMS threshold after the peak,
  # DoubleSigmoid_Offset_Time can be NA. In that case, use the fitted right
  # minimum of the double-sigmoid curve as a biologically meaningful fallback:
  #   DoubleSigmoid_Decay_Time = DoubleSigmoid_Right_Min_Time - detected peak time
  #
  # This keeps the decay-time column from appearing empty when the model fit is
  # otherwise successful but the threshold-crossing offset is not found.
  # ------------------------------------------------------------
  
  ds_duration <- ifelse(
    is.finite(ds_onset) && is.finite(ds_offset),
    ds_offset - ds_onset,
    NA_real_
  )
  
  ds_rise_time <- ifelse(
    is.finite(ds_onset) && is.finite(peak_time),
    peak_time - ds_onset,
    NA_real_
  )
  
  ds_decay_time <- ifelse(
    is.finite(ds_offset) && is.finite(peak_time),
    ds_offset - peak_time,
    NA_real_
  )
  
  # Fallback decay time based on fitted right minimum
  if (!is.finite(ds_decay_time) &&
      is.finite(right_min_time) &&
      is.finite(peak_time) &&
      right_min_time > peak_time) {
    ds_decay_time <- right_min_time - peak_time
  }
  
  # Fallback duration based on fitted right minimum
  if (!is.finite(ds_duration) &&
      is.finite(ds_onset) &&
      is.finite(right_min_time) &&
      right_min_time > ds_onset) {
    ds_duration <- right_min_time - ds_onset
  }
  
  # Parameter-based slopes. Decay slope is negative by convention.
  rise_slope <- ifelse(
    is.finite(pars[["A"]]) && is.finite(pars[["k1"]]),
    pars[["A"]] * pars[["k1"]] / 4,
    NA_real_
  )
  
  decay_slope <- ifelse(
    is.finite(pars[["A"]]) && is.finite(pars[["k2"]]),
    -pars[["A"]] * pars[["k2"]] / 4,
    NA_real_
  )
  
  # Numerical fallback for decay slope if the parameter-based value is missing.
  # It uses the steepest negative slope of the fitted double-sigmoid curve
  # after the fitted peak.
  if (!is.finite(decay_slope)) {
    dy <- diff(dense_y) / diff(dense_t)
    t_mid <- (head(dense_t, -1) + tail(dense_t, -1)) / 2
    right_deriv_idx <- which(t_mid >= model_peak_time & is.finite(dy))
    if (length(right_deriv_idx) > 0) {
      decay_slope <- min(dy[right_deriv_idx], na.rm = TRUE)
    }
  }
  
  list(
    DoubleSigmoid_Onset_Time = ds_onset,
    DoubleSigmoid_Offset_Time = ds_offset,
    DoubleSigmoid_Duration = ds_duration,
    DoubleSigmoid_Rise_Time = ds_rise_time,
    DoubleSigmoid_Decay_Time = ds_decay_time,
    DoubleSigmoid_Baseline = Baseline_fixed,
    DoubleSigmoid_RMS = baseline_info$RMS,
    DoubleSigmoid_Threshold = threshold,
    DoubleSigmoid_A = pars[["A"]],
    DoubleSigmoid_k1 = pars[["k1"]],
    DoubleSigmoid_t1 = pars[["t1"]],
    DoubleSigmoid_k2 = pars[["k2"]],
    DoubleSigmoid_t2 = pars[["t2"]],
    DoubleSigmoid_Rise_Slope = rise_slope,
    DoubleSigmoid_Decay_Slope = decay_slope,
    DoubleSigmoid_Model_Peak_Time = model_peak_time,
    DoubleSigmoid_Model_Peak_Value = model_peak_value,
    DoubleSigmoid_Left_Min_Time = left_min_time,
    DoubleSigmoid_Left_Min_Value = left_min_value,
    DoubleSigmoid_Right_Min_Time = right_min_time,
    DoubleSigmoid_Right_Min_Value = right_min_value,
    DoubleSigmoid_Baseline_at_peak = ds_baseline_at_peak,
    DoubleSigmoid_Prominence = ds_prominence,
    DoubleSigmoid_R2 = r2,
    DoubleSigmoid_Fit_Status = "OK",
    fit_df = data.frame(
      Event_ID = NA_integer_,
      Time = dense_t,
      Fitted = dense_y,
      Model = "Double sigmoid",
      stringsAsFactors = FALSE
    )
  )
}

append_fft_sigmoid_metrics <- function(metrics,
                                       time,
                                       signal,
                                       use_single_sigmoid = FALSE,
                                       use_double_sigmoid = FALSE,
                                       pre_event_window = 5,
                                       single_left_pad = 0,
                                       single_right_pad = 0,
                                       double_left_pad = 0,
                                       double_right_pad = 0) {
  if (is.null(metrics) || nrow(metrics) == 0) {
    return(list(metrics = metrics, fitted_curves = data.frame()))
  }
  
  # If no sigmoid option is selected, preserve the original workflow exactly.
  if (!isTRUE(use_single_sigmoid) && !isTRUE(use_double_sigmoid)) {
    return(list(metrics = metrics, fitted_curves = data.frame()))
  }
  
  curve_list <- list()
  single_rows <- vector("list", nrow(metrics))
  double_rows <- vector("list", nrow(metrics))
  
  # Empty/status rows used when a detected candidate is not a valid
  # biological transient for sigmoid modeling. In particular, if the
  # detected local maximum is at or below the selected baseline, the app
  # must not draw a single- or double-sigmoid curve for that candidate.
  empty_single_sigmoid_row <- function(status = "Not requested") {
    data.frame(
      SingleSigmoid_Onset_Time = NA_real_,
      SingleSigmoid_Rise_Time = NA_real_,
      SingleSigmoid_Rise_Rate = NA_real_,
      SingleSigmoid_Baseline = NA_real_,
      SingleSigmoid_RMS = NA_real_,
      SingleSigmoid_Threshold = NA_real_,
      SingleSigmoid_A = NA_real_,
      SingleSigmoid_k = NA_real_,
      SingleSigmoid_t0 = NA_real_,
      SingleSigmoid_Slope_at_t0 = NA_real_,
      SingleSigmoid_R2 = NA_real_,
      SingleSigmoid_Fit_Status = status,
      stringsAsFactors = FALSE
    )
  }
  
  empty_double_sigmoid_row <- function(status = "Not requested") {
    data.frame(
      DoubleSigmoid_Onset_Time = NA_real_,
      DoubleSigmoid_Offset_Time = NA_real_,
      DoubleSigmoid_Duration = NA_real_,
      DoubleSigmoid_Rise_Time = NA_real_,
      DoubleSigmoid_Decay_Time = NA_real_,
      DoubleSigmoid_Baseline = NA_real_,
      DoubleSigmoid_RMS = NA_real_,
      DoubleSigmoid_Threshold = NA_real_,
      DoubleSigmoid_A = NA_real_,
      DoubleSigmoid_k1 = NA_real_,
      DoubleSigmoid_t1 = NA_real_,
      DoubleSigmoid_k2 = NA_real_,
      DoubleSigmoid_t2 = NA_real_,
      DoubleSigmoid_Rise_Slope = NA_real_,
      DoubleSigmoid_Decay_Slope = NA_real_,
      DoubleSigmoid_Model_Peak_Time = NA_real_,
      DoubleSigmoid_Model_Peak_Value = NA_real_,
      DoubleSigmoid_Left_Min_Time = NA_real_,
      DoubleSigmoid_Left_Min_Value = NA_real_,
      DoubleSigmoid_Right_Min_Time = NA_real_,
      DoubleSigmoid_Right_Min_Value = NA_real_,
      DoubleSigmoid_Baseline_at_peak = NA_real_,
      DoubleSigmoid_Prominence = NA_real_,
      DoubleSigmoid_R2 = NA_real_,
      DoubleSigmoid_Fit_Status = status,
      stringsAsFactors = FALSE
    )
  }
  
  for (i in seq_len(nrow(metrics))) {
    
    peak_time <- metrics$Peak_Occurence_Time[i]
    peak_value <- metrics$Amplitude[i]
    
    fallback_baseline <- if ("Baseline_at_peak" %in% names(metrics)) {
      metrics$Baseline_at_peak[i]
    } else {
      NA_real_
    }
    
    # Stable event ID used for plotting fitted curves.
    # This ID must match the original peak-candidate/Event_ID used by
    # the event-specific baseline object. Do not use the row number after
    # filtering, because removed peaks would shift the IDs and could make
    # a baseline or sigmoid curve appear on a rejected transient.
    event_id_for_output <- i
    if (".Peak_Candidate_ID" %in% names(metrics)) {
      tmp_event_id <- suppressWarnings(as.integer(metrics$.Peak_Candidate_ID[i]))
      if (length(tmp_event_id) > 0 && is.finite(tmp_event_id[1])) {
        event_id_for_output <- tmp_event_id[1]
      }
    }
    
    # ------------------------------------------------------------
    # Do not fit sigmoid models to candidates that are not above baseline.
    # This prevents blue/purple fitted curves from being drawn on small
    # local bumps located below the selected baseline. Those candidates are
    # not treated as valid calcium transients for single/double sigmoid
    # modeling.
    # ------------------------------------------------------------
    if (is.finite(peak_value) && is.finite(fallback_baseline) &&
        peak_value <= fallback_baseline) {
      single_rows[[i]] <- empty_single_sigmoid_row(
        if (isTRUE(use_single_sigmoid)) "Skipped: peak at or below baseline" else "Not requested"
      )
      double_rows[[i]] <- empty_double_sigmoid_row(
        if (isTRUE(use_double_sigmoid)) "Skipped: peak at or below baseline" else "Not requested"
      )
      next
    }
    
    # ------------------------------------------------------------
    # Left reference for sigmoid fitting
    # ------------------------------------------------------------
    # Prefer event-specific left minimum if available.
    # Otherwise use the original Transient_Ocurrence_Time.
    # ------------------------------------------------------------
    left_ref <- if ("Baseline_Left_Min_Time" %in% names(metrics) &&
                    is.finite(metrics$Baseline_Left_Min_Time[i])) {
      metrics$Baseline_Left_Min_Time[i]
    } else if ("Transient_Ocurrence_Time" %in% names(metrics) &&
               is.finite(metrics$Transient_Ocurrence_Time[i])) {
      metrics$Transient_Ocurrence_Time[i]
    } else {
      NA_real_
    }
    
    # ------------------------------------------------------------
    # Right reference for double sigmoid fitting
    # ------------------------------------------------------------
    # The double sigmoid should cover the full transient decay. Therefore,
    # do not stop the fitted purple curve at Time_right_FWHM/FWHP unless no
    # better decay endpoint is available. First estimate an automatic
    # decay-end time from the smoothed trace, clipped by the next peak.
    peak_times_all <- suppressWarnings(as.numeric(metrics$Peak_Occurence_Time))
    future_peaks <- peak_times_all[is.finite(peak_times_all) & peak_times_all > peak_time]
    next_peak_time <- if (length(future_peaks) > 0) min(future_peaks, na.rm = TRUE) else NA_real_
    
    min_right_candidates <- c(peak_time)
    if ("Time_right_FWHM" %in% names(metrics) && is.finite(metrics$Time_right_FWHM[i])) {
      min_right_candidates <- c(min_right_candidates, metrics$Time_right_FWHM[i])
    }
    if ("Time_right_FWHP" %in% names(metrics) && is.finite(metrics$Time_right_FWHP[i])) {
      min_right_candidates <- c(min_right_candidates, metrics$Time_right_FWHP[i])
    }
    min_right_time <- max(min_right_candidates, na.rm = TRUE)
    
    right_ref_auto <- fft_find_double_sigmoid_decay_end(
      time = time,
      signal = signal,
      peak_time = peak_time,
      left_reference_time = left_ref,
      fallback_baseline = fallback_baseline,
      pre_event_window = pre_event_window,
      next_peak_time = next_peak_time,
      min_right_time = min_right_time
    )
    
    right_ref_fallback <- if ("Baseline_AUC_Right_Time" %in% names(metrics) &&
                              is.finite(metrics$Baseline_AUC_Right_Time[i])) {
      metrics$Baseline_AUC_Right_Time[i]
    } else if ("Baseline_Right_Min_Time" %in% names(metrics) &&
               is.finite(metrics$Baseline_Right_Min_Time[i])) {
      metrics$Baseline_Right_Min_Time[i]
    } else if ("Time_right_FWHM" %in% names(metrics) &&
               is.finite(metrics$Time_right_FWHM[i])) {
      metrics$Time_right_FWHM[i]
    } else if ("Time_right_FWHP" %in% names(metrics) &&
               is.finite(metrics$Time_right_FWHP[i])) {
      metrics$Time_right_FWHP[i]
    } else {
      NA_real_
    }
    
    right_ref <- if (is.finite(right_ref_auto) && right_ref_auto > peak_time) {
      right_ref_auto
    } else {
      right_ref_fallback
    }
    
    # ============================================================
    # Single sigmoid: optional onset/rise-time replacement
    # ============================================================
    if (isTRUE(use_single_sigmoid)) {
      
      single_fit <- fit_fft_single_sigmoid_event(
        time = time,
        signal = signal,
        left_reference_time = left_ref,
        peak_time = peak_time,
        peak_value = peak_value,
        fallback_baseline = fallback_baseline,
        pre_event_window = pre_event_window,
        fit_left_pad = single_left_pad,
        fit_right_pad = single_right_pad
      )
      
      single_rows[[i]] <- data.frame(
        SingleSigmoid_Onset_Time = single_fit$SingleSigmoid_Onset_Time,
        SingleSigmoid_Rise_Time = single_fit$SingleSigmoid_Rise_Time,
        SingleSigmoid_Rise_Rate = single_fit$SingleSigmoid_Rise_Rate,
        SingleSigmoid_Baseline = single_fit$SingleSigmoid_Baseline,
        SingleSigmoid_RMS = single_fit$SingleSigmoid_RMS,
        SingleSigmoid_Threshold = single_fit$SingleSigmoid_Threshold,
        SingleSigmoid_A = single_fit$SingleSigmoid_A,
        SingleSigmoid_k = single_fit$SingleSigmoid_k,
        SingleSigmoid_t0 = single_fit$SingleSigmoid_t0,
        SingleSigmoid_Slope_at_t0 = single_fit$SingleSigmoid_Slope_at_t0,
        SingleSigmoid_R2 = single_fit$SingleSigmoid_R2,
        SingleSigmoid_Fit_Status = single_fit$SingleSigmoid_Fit_Status,
        stringsAsFactors = FALSE
      )
      
      if (!is.null(single_fit$fit_df) && nrow(single_fit$fit_df) > 0) {
        single_fit$fit_df$Event_ID <- event_id_for_output
        curve_list[[length(curve_list) + 1]] <- single_fit$fit_df
      }
      
    } else {
      
      single_rows[[i]] <- data.frame(
        SingleSigmoid_Onset_Time = NA_real_,
        SingleSigmoid_Rise_Time = NA_real_,
        SingleSigmoid_Rise_Rate = NA_real_,
        SingleSigmoid_Baseline = NA_real_,
        SingleSigmoid_RMS = NA_real_,
        SingleSigmoid_Threshold = NA_real_,
        SingleSigmoid_A = NA_real_,
        SingleSigmoid_k = NA_real_,
        SingleSigmoid_t0 = NA_real_,
        SingleSigmoid_Slope_at_t0 = NA_real_,
        SingleSigmoid_R2 = NA_real_,
        SingleSigmoid_Fit_Status = "Not requested",
        stringsAsFactors = FALSE
      )
    }
    
    # ============================================================
    # Double sigmoid: optional full-transient metrics
    # ============================================================
    if (isTRUE(use_double_sigmoid)) {
      
      double_fit <- fit_fft_double_sigmoid_event(
        time = time,
        signal = signal,
        left_reference_time = left_ref,
        peak_time = peak_time,
        right_reference_time = right_ref,
        peak_value = peak_value,
        fallback_baseline = fallback_baseline,
        pre_event_window = pre_event_window,
        fit_left_pad = double_left_pad,
        fit_right_pad = double_right_pad
      )
      
      double_rows[[i]] <- data.frame(
        DoubleSigmoid_Onset_Time = double_fit$DoubleSigmoid_Onset_Time,
        DoubleSigmoid_Offset_Time = double_fit$DoubleSigmoid_Offset_Time,
        DoubleSigmoid_Duration = double_fit$DoubleSigmoid_Duration,
        DoubleSigmoid_Rise_Time = double_fit$DoubleSigmoid_Rise_Time,
        DoubleSigmoid_Decay_Time = double_fit$DoubleSigmoid_Decay_Time,
        DoubleSigmoid_Baseline = double_fit$DoubleSigmoid_Baseline,
        DoubleSigmoid_RMS = double_fit$DoubleSigmoid_RMS,
        DoubleSigmoid_Threshold = double_fit$DoubleSigmoid_Threshold,
        DoubleSigmoid_A = double_fit$DoubleSigmoid_A,
        DoubleSigmoid_k1 = double_fit$DoubleSigmoid_k1,
        DoubleSigmoid_t1 = double_fit$DoubleSigmoid_t1,
        DoubleSigmoid_k2 = double_fit$DoubleSigmoid_k2,
        DoubleSigmoid_t2 = double_fit$DoubleSigmoid_t2,
        DoubleSigmoid_Rise_Slope = double_fit$DoubleSigmoid_Rise_Slope,
        DoubleSigmoid_Decay_Slope = double_fit$DoubleSigmoid_Decay_Slope,
        DoubleSigmoid_Model_Peak_Time = double_fit$DoubleSigmoid_Model_Peak_Time,
        DoubleSigmoid_Model_Peak_Value = double_fit$DoubleSigmoid_Model_Peak_Value,
        DoubleSigmoid_Left_Min_Time = double_fit$DoubleSigmoid_Left_Min_Time,
        DoubleSigmoid_Left_Min_Value = double_fit$DoubleSigmoid_Left_Min_Value,
        DoubleSigmoid_Right_Min_Time = double_fit$DoubleSigmoid_Right_Min_Time,
        DoubleSigmoid_Right_Min_Value = double_fit$DoubleSigmoid_Right_Min_Value,
        DoubleSigmoid_Baseline_at_peak = double_fit$DoubleSigmoid_Baseline_at_peak,
        DoubleSigmoid_Prominence = double_fit$DoubleSigmoid_Prominence,
        DoubleSigmoid_R2 = double_fit$DoubleSigmoid_R2,
        DoubleSigmoid_Fit_Status = double_fit$DoubleSigmoid_Fit_Status,
        stringsAsFactors = FALSE
      )
      
      if (!is.null(double_fit$fit_df) && nrow(double_fit$fit_df) > 0) {
        double_fit$fit_df$Event_ID <- event_id_for_output
        curve_list[[length(curve_list) + 1]] <- double_fit$fit_df
      }
      
    } else {
      
      double_rows[[i]] <- data.frame(
        DoubleSigmoid_Onset_Time = NA_real_,
        DoubleSigmoid_Offset_Time = NA_real_,
        DoubleSigmoid_Duration = NA_real_,
        DoubleSigmoid_Rise_Time = NA_real_,
        DoubleSigmoid_Decay_Time = NA_real_,
        DoubleSigmoid_Baseline = NA_real_,
        DoubleSigmoid_RMS = NA_real_,
        DoubleSigmoid_Threshold = NA_real_,
        DoubleSigmoid_A = NA_real_,
        DoubleSigmoid_k1 = NA_real_,
        DoubleSigmoid_t1 = NA_real_,
        DoubleSigmoid_k2 = NA_real_,
        DoubleSigmoid_t2 = NA_real_,
        DoubleSigmoid_Rise_Slope = NA_real_,
        DoubleSigmoid_Decay_Slope = NA_real_,
        DoubleSigmoid_Model_Peak_Time = NA_real_,
        DoubleSigmoid_Model_Peak_Value = NA_real_,
        DoubleSigmoid_Left_Min_Time = NA_real_,
        DoubleSigmoid_Left_Min_Value = NA_real_,
        DoubleSigmoid_Right_Min_Time = NA_real_,
        DoubleSigmoid_Right_Min_Value = NA_real_,
        DoubleSigmoid_Baseline_at_peak = NA_real_,
        DoubleSigmoid_Prominence = NA_real_,
        DoubleSigmoid_R2 = NA_real_,
        DoubleSigmoid_Fit_Status = "Not requested",
        stringsAsFactors = FALSE
      )
    }
  }
  
  single_df <- do.call(rbind, single_rows)
  double_df <- do.call(rbind, double_rows)
  
  fitted_curves <- if (length(curve_list) > 0) {
    do.call(rbind, curve_list)
  } else {
    data.frame()
  }
  
  metrics_out <- cbind(metrics, single_df, double_df)
  
  # ============================================================
  # User-requested behavior
  # ============================================================
  # When "10.1. Apply Single Sigmoid Fit for Onset and Rise-Time Metrics"
  # is selected, replace the main onset-related columns with the
  # single-sigmoid values:
  #   Transient_Ocurrence_Time <- SingleSigmoid_Onset_Time
  #   Peak_Rise_Time <- Peak_Occurence_Time - SingleSigmoid_Onset_Time
  #   Rise_Rate <- amplitude above baseline / Peak_Rise_Time
  #
  # When the checkbox is not selected, the function returns earlier or
  # keeps the original columns unchanged.
  # ============================================================
  
  if (isTRUE(use_single_sigmoid) &&
      "SingleSigmoid_Onset_Time" %in% names(metrics_out)) {
    
    valid_single <- is.finite(metrics_out$SingleSigmoid_Onset_Time) &
      metrics_out$SingleSigmoid_Fit_Status == "OK"
    
    # Keep original values for comparison and transparency.
    if ("Transient_Ocurrence_Time" %in% names(metrics_out)) {
      metrics_out$Original_Transient_Ocurrence_Time <-
        metrics_out$Transient_Ocurrence_Time
    }
    
    if ("Peak_Rise_Time" %in% names(metrics_out)) {
      metrics_out$Original_Peak_Rise_Time <-
        metrics_out$Peak_Rise_Time
    }
    
    if ("Rise_Rate" %in% names(metrics_out)) {
      metrics_out$Original_Rise_Rate <-
        metrics_out$Rise_Rate
    }
    
    metrics_out$Onset_Method <- ifelse(
      valid_single,
      "Single sigmoid + RMS",
      "Original local-minimum onset"
    )
    
    # Replace the main onset column only for successful fits.
    metrics_out$Transient_Ocurrence_Time[valid_single] <-
      metrics_out$SingleSigmoid_Onset_Time[valid_single]
    
    # Recalculate rise time using the new sigmoid-based onset.
    metrics_out$Peak_Rise_Time[valid_single] <-
      metrics_out$Peak_Occurence_Time[valid_single] -
      metrics_out$Transient_Ocurrence_Time[valid_single]
    
    # Recalculate rise rate.
    amp_above_baseline <- metrics_out$Amplitude
    
    if ("Baseline_at_peak" %in% names(metrics_out)) {
      amp_above_baseline <- metrics_out$Amplitude - metrics_out$Baseline_at_peak
    }
    
    metrics_out$Rise_Rate[valid_single] <- ifelse(
      is.finite(metrics_out$Peak_Rise_Time[valid_single]) &
        metrics_out$Peak_Rise_Time[valid_single] > 0,
      amp_above_baseline[valid_single] /
        metrics_out$Peak_Rise_Time[valid_single],
      NA_real_
    )
    
    # Keep the single-sigmoid-specific columns consistent with the main metrics.
    metrics_out$SingleSigmoid_Rise_Time[valid_single] <-
      metrics_out$Peak_Rise_Time[valid_single]
    
    metrics_out$SingleSigmoid_Rise_Rate[valid_single] <-
      metrics_out$Rise_Rate[valid_single]
    
  } else {
    # If only double sigmoid is requested, keep original onset metrics.
    if (!"Onset_Method" %in% names(metrics_out)) {
      metrics_out$Onset_Method <- "Original local-minimum onset"
    }
  }
  
  list(
    metrics = metrics_out,
    fitted_curves = fitted_curves
  )
}


# ------------------------------------------------------------
# Deprecated helper: build a baseline object from double-sigmoid minima
# ------------------------------------------------------------
# Kept for backward compatibility, but the corrected FFT workflow no longer
# uses double-sigmoid minima to redefine the event-specific baseline.
# When baseline mode is "Event-Specific Baseline: Left-Right FFT Minima",
# the plotted baseline line and AUC baseline must come from the original
# calcium-transient left/right FFT minima computed by
# compute_event_specific_minima_baseline().
# ------------------------------------------------------------
build_double_sigmoid_event_baseline_info <- function(metrics, time, signal) {
  empty <- list(
    scalar = NA_real_,
    trace = data.frame(Time = numeric(0), Baseline = numeric(0)),
    event_baselines = data.frame(),
    mode = "8_double_sigmoid"
  )
  
  if (is.null(metrics) || nrow(metrics) == 0) return(empty)
  required_cols <- c(
    "Peak_Occurence_Time",
    "DoubleSigmoid_Left_Min_Time", "DoubleSigmoid_Left_Min_Value",
    "DoubleSigmoid_Right_Min_Time", "DoubleSigmoid_Right_Min_Value",
    "DoubleSigmoid_Baseline_at_peak"
  )
  if (!all(required_cols %in% names(metrics))) return(empty)
  
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  
  eb_list <- lapply(seq_len(nrow(metrics)), function(i) {
    lt <- suppressWarnings(as.numeric(metrics$DoubleSigmoid_Left_Min_Time[i]))
    lv <- suppressWarnings(as.numeric(metrics$DoubleSigmoid_Left_Min_Value[i]))
    rt <- suppressWarnings(as.numeric(metrics$DoubleSigmoid_Right_Min_Time[i]))
    rv <- suppressWarnings(as.numeric(metrics$DoubleSigmoid_Right_Min_Value[i]))
    pt <- suppressWarnings(as.numeric(metrics$Peak_Occurence_Time[i]))
    bp <- suppressWarnings(as.numeric(metrics$DoubleSigmoid_Baseline_at_peak[i]))
    
    if (!is.finite(lt) || !is.finite(lv) ||
        !is.finite(rt) || !is.finite(rv) ||
        !is.finite(pt) || !is.finite(bp) ||
        rt <= lt) {
      return(NULL)
    }
    
    data.frame(
      Event_ID = i,
      Peak_Time = pt,
      Left_Min_Time = lt,
      Left_Min_Value = lv,
      Right_Min_Time = rt,
      Right_Min_Value = rv,
      Plot_Left_Time = lt,
      Plot_Left_Value = lv,
      Plot_Right_Time = rt,
      Plot_Right_Value = rv,
      Baseline_at_peak = bp,
      stringsAsFactors = FALSE
    )
  })
  
  eb_list <- eb_list[!sapply(eb_list, is.null)]
  if (length(eb_list) == 0) return(empty)
  event_baselines <- do.call(rbind, eb_list)
  
  trace_list <- lapply(seq_len(nrow(event_baselines)), function(i) {
    idx <- which(time >= event_baselines$Plot_Left_Time[i] &
                   time <= event_baselines$Plot_Right_Time[i])
    if (length(idx) < 2) return(NULL)
    
    baseline_vec <- stats::approx(
      x = c(event_baselines$Left_Min_Time[i], event_baselines$Right_Min_Time[i]),
      y = c(event_baselines$Left_Min_Value[i], event_baselines$Right_Min_Value[i]),
      xout = time[idx],
      rule = 2
    )$y
    
    data.frame(
      Event_ID = event_baselines$Event_ID[i],
      Time = time[idx],
      Baseline = baseline_vec,
      Signal = signal[idx],
      AUC_Upper = pmax(signal[idx], baseline_vec),
      stringsAsFactors = FALSE
    )
  })
  
  trace_list <- trace_list[!sapply(trace_list, is.null)]
  baseline_trace <- if (length(trace_list) > 0) do.call(rbind, trace_list) else empty$trace
  
  list(
    scalar = mean(event_baselines$Baseline_at_peak, na.rm = TRUE),
    trace = baseline_trace,
    event_baselines = event_baselines,
    mode = "8_double_sigmoid"
  )
}


# ------------------------------------------------------------
# Freeze FWHP and FWHM before optional sigmoid modeling
# ------------------------------------------------------------
# The single sigmoid is allowed to refine onset/rise-time metrics.
# The double sigmoid is allowed to add full-transient/decay metrics and,
# when requested, to update prominence. It must not redefine the event-specific baseline.
# However, FWHP and FWHM must remain the values obtained from the original
# FFT-smoothed trace + local-minimum/baseline workflow. This helper restores
# those width columns after optional sigmoid modeling so they never depend on
# either the single sigmoid or the double sigmoid.
restore_sigmoid_independent_width_metrics <- function(metrics) {
  if (is.null(metrics) || nrow(metrics) == 0) return(metrics)
  
  restore_col <- function(df, target_col, frozen_col) {
    if (frozen_col %in% names(df)) {
      df[[target_col]] <- df[[frozen_col]]
    }
    df
  }
  
  # FWHP geometry and FWHP numeric value.
  metrics <- restore_col(metrics, "FWHP", "FWHP_NoSigmoid")
  metrics <- restore_col(metrics, "Time_left_FWHP", "Time_left_FWHP_NoSigmoid")
  metrics <- restore_col(metrics, "Time_right_FWHP", "Time_right_FWHP_NoSigmoid")
  
  # FWHM geometry and FWHM numeric value.
  metrics <- restore_col(metrics, "FWHM", "FWHM_NoSigmoid")
  metrics <- restore_col(metrics, "Time_left_FWHM", "Time_left_FWHM_NoSigmoid")
  metrics <- restore_col(metrics, "Time_right_FWHM", "Time_right_FWHM_NoSigmoid")
  metrics <- restore_col(metrics, "Amplitude_Midpoint", "FWHM_Level_NoSigmoid")
  
  metrics$Width_Metrics_Method <- "Original FFT/local-minimum workflow; independent of single and double sigmoid"
  metrics
}

# =========================
# UI Module
# =========================
mod_Denoising_data_fft_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$head(
      tags$style(HTML("
        .small-button { font-size: 10px; padding: 2px 4px; }
        .helper-card {
          background: #f7f7f7;
          border: 1px solid #e5e5e5;
          border-radius: 10px;
          padding: 12px;
          margin-top: 8px;
        }
        .mono {
          font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace;
        }
        .formula-box {
          background: #fbfbfb;
          border: 1px solid #d9d9d9;
          border-radius: 10px;
          padding: 12px;
          white-space: pre-wrap;
          word-break: break-word;
          font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace;
          font-size: 13px;
          line-height: 1.5;
          max-height: 320px;
          overflow-y: auto;
        }
      "))
    ),
    
    sidebarLayout(
      sidebarPanel(
        width = 4,
        shinyjs::useShinyjs(),
        
        # actionButton(
        #   ns("param_info_button11"), "Help",
        #   class = "btn-sm",
        #   style = "position:absolute; top:0; right:15px; margin:5px;"
        # ),
        
        radioButtons(
          ns("data_simulate"), "1. Example Data",
          choices = c("Yes" = 1, "No" = 0),
          selected = 1
        ),
        
        fileInput(
          ns("fileBcsv2"),
          accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".tsv"),
          label = h5("2. Dataset")
        ),
        
        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),
        
        numericInput(ns("Cell2"), "3. Region of Interest (ROI):", value = 1, min = 1),
        
        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),
        
        # tags$h4("Trend Removal", style = "color: gray; margin-top: 10px;"),
        # checkboxInput(
        #   ns("remove_trend"),
        #   "Remove polynomial trend from calcium trace",
        #   value = FALSE
        # ),
        conditionalPanel(
          condition = "input.remove_trend == true",
          ns = ns,
          sliderInput(
            ns("trend_degree"),
            "Polynomial degree:",
            min = 1,
            max = 6,
            value = 2,
            step = 1
          ),
          # checkboxInput(
          #   ns("trend_center_scale_t"),
          #   "Center and scale t before fitting trend",
          #   value = TRUE
          # ),
          helpText("When selected, the app fits a polynomial trend and uses the residuals as the new calcium trace.")
        ),
        
        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),
        
        sliderInput(ns("fft_fraction"), "4. FFT Low-Frequency Fraction:", min = 0.01, max = 1,
                    value = 0.20, step = 0.01),
        checkboxInput(ns("keep_mean"), "4.1. Keep Mean (DC Component)", TRUE),
        
        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),
        
        tags$h4("5. Find Peaks Function Arguments", style = "color: gray; margin-top: 10px;"),
        
        fluidRow(
          column(
            width = 6,
            numericInput(ns("minpeakheight2"), "5.1. Peak Height (min)", value = 0, min = 0, max = 100, step = 0.1),
            numericInput(ns("ndowns2"), "5.3. Peak Descent", value = 1, min = 0, max = 100)
          ),
          column(
            width = 6,
            numericInput(ns("nups2"), "5.2. Peak Ascent", value = 1, min = 0, max = 100),
            numericInput(ns("minpeakdistance2"), "5.4. Min Peak Distance", value = 1, min = 0, max = 100)
          ),
          column(
            width = 6,
            numericInput(ns("min_FWHP"), "5.5. FWHP (min)", value = 0, min = 0, step = 0.1)
          ),
          column(
            width = 6,
            numericInput(ns("min_prominence"), "5.6. Prominence (min)", value = 0, min = 0, step = 0.1)
          )
        ),
        
        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),
        
        selectInput(
          ns("Baseline"), "6. Baseline:",
          choices = c(
            "6.1. Standard Definition" = 2,
            "6.2. Low-Fluorescence Region Baseline" = 3,
            "6.3. Rolling Percentile Baseline" = 6,
            "6.4. Event-Specific Baseline: Left-Right FFT Minima" = 8,
            "6.5. Constant Baseline: Mean From Trace Start to Time_Onset" = 11,
            "6.6. Median Signal (Legacy)" = 9,
            "6.7. Mean Signal (Legacy)" = 10,
            "6.8. Min" = 5
          ),
          selected = 2
        ),
        
        conditionalPanel(
          condition = "input.Baseline==2", ns = ns,
          numericInput(
            ns("initial_trace_fraction"),
            "6.1.1. Initial Trace Fraction:",
            value = 0.10,
            min = 0.01,
            max = 1,
            step = 0.01
          ),
          helpText("Uses the initial fraction of the trace to estimate a global baseline.")
        ),
        
        conditionalPanel(
          condition = "input.Baseline==3", ns = ns,
          numericInput(
            ns("minimal_lower_fraction"),
            "6.2.1. Minimal Lower Fraction:",
            value = 0.10,
            min = 0.01,
            max = 1,
            step = 0.01
          ),
          helpText("Uses the lowest fraction of fluorescence values to estimate the baseline.")
        ),
        
        conditionalPanel(
          condition = "input.Baseline==6 || input.Baseline==7", ns = ns,
          fluidRow(
            column(
              6,
              numericInput(
                ns("rolling_window"),
                "6.3.1. Moving Window Size:",
                value = 20,
                min = 0.0001,
                step = 1
              )
            ),
            column(
              6,
              numericInput(
                ns("rolling_percentile"),
                "6.3.2. Percentile:",
                value = 2,
                min = 0,
                max = 100,
                step = 0.1
              )
            )
          ),
          helpText("Use the same time unit as the Time column. For example, if Time is in seconds, 10 minutes = 600 seconds.")
        ),
        
        conditionalPanel(
          condition = "input.Baseline==7", ns = ns,
          checkboxInput(
            ns("show_fixed_windows"),
            "6.3.3. Show Fixed-Window Intervals on Calcium Trace",
            value = TRUE
          )
          
        ),
        
        conditionalPanel(
          condition = "input.Baseline==8", ns = ns,
          sliderInput(
            ns("event_baseline_fft_fraction"),
            "6.4.1. FFT Low-Frequency Fraction for Event-Specific Baseline:",
            min = 0.01,
            max = 1,
            value = 0.20,
            step = 0.01
          ),
          helpText(
            "This controls only the FFT smoothing used to find the left and right minima for the event-specific baseline. Peak detection and the main calcium trace smoothing still use section 4."
          )
        ),
        
        conditionalPanel(
          condition = "input.Baseline==11", ns = ns,
          helpText(
            "This option finds the first Time_Onset among the detected transients, takes the calcium trace from the beginning up to that Time_Onset, calculates the mean of that portion, and uses that single constant value as the baseline for the complete trace and all detected transients."
          )
        ),
        
        radioButtons(ns("auc2"), "7. Area Under the Curve (AUC):", choices = c("No" = 1, "Yes" = 2), selected = 1),
        radioButtons(ns("raw_data"), "8. Raw Data", choices = c("No" = 1, "Yes" = 2), selected = 2),
        radioButtons(ns("FWHM"), "9. Full Width at Half Maximum:", choices = c("No" = 1, "Yes" = 2), selected = 1),
        
        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),
        tags$h4("10. Optional Sigmoid Modeling", style = "color: gray; margin-top: 10px;"),
        tags$div(
          class = "helper-card",
          tags$strong("10. Optional Sigmoid Fits"),
          tags$p(
            "These models do not change the FFT smoothing, peak detection, baseline, AUC, FWHM, or FWHP calculations. ",
            "They only add extra metrics when activated. The single sigmoid estimates onset and rise-time metrics; ",
            "the double sigmoid models the full transient and adds decay/offset metrics."
          )
        ),
        checkboxInput(
          ns("use_single_sigmoid_fft"),
          "10.1. Apply Single Sigmoid Fit for Onset and Rise-Time Metrics",
          value = FALSE
        ),
        checkboxInput(
          ns("use_double_sigmoid_fft"),
          "10.2. Apply Double Sigmoid Fit for Full-Transient Metrics",
          value = FALSE
        ),
        conditionalPanel(
          condition = "input.use_single_sigmoid_fft == true || input.use_double_sigmoid_fft == true",
          ns = ns,
          numericInput(
            ns("sigmoid_pre_event_window_fft"),
            "10.3. Pre-Event Window for Baseline RMS:",
            value = 5,
            min = 0.1,
            step = 0.5
          )
        ),
        conditionalPanel(
          condition = "input.use_single_sigmoid_fft == true",
          ns = ns,
          fluidRow(
            column(
              6,
              numericInput(
                ns("single_sigmoid_left_pad_fft"),
                "10.4. Single Sigmoid Left Padding:",
                value = 0,
                min = 0,
                step = 0.25
              )
            ),
            column(
              6,
              numericInput(
                ns("single_sigmoid_right_pad_fft"),
                "10.5. Single Sigmoid Right Padding:",
                value = 0,
                min = 0,
                step = 0.25
              )
            )
          )
        ),
        conditionalPanel(
          condition = "input.use_double_sigmoid_fft == true",
          ns = ns,
          fluidRow(
            column(
              6,
              numericInput(
                ns("double_sigmoid_left_pad_fft"),
                "10.6. Double Sigmoid Left Padding:",
                value = 0,
                min = 0,
                step = 0.25
              )
            ),
            column(
              6,
              numericInput(
                ns("double_sigmoid_right_pad_fft"),
                "10.7. Double Sigmoid Right Padding:",
                value = 0,
                min = 0,
                step = 0.25
              )
            )
          )
        ),
        
        tags$h4("11. Download Options", style = "color: gray; margin-top: 10px;"),
        downloadButton(ns("descargarP"), "11.1. Trace Metrics"),
        downloadButton(ns("descargar"), "11.2. Transient Metrics"),
        downloadButton(ns("Calcium_Trance_Graph"), "11.3. Calcium Trace Graph")
      ),
      
      mainPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "12. SummaryData",
            tags$h4("12.1. Dataset Summary", style = "color: gray; margin-top: 10px;"),
            DTOutput(ns("infotable2")),
            tags$br(),
            tags$h4("12.2. Dataset", style = "color: gray; margin-top: 10px;"),
            DTOutput(ns("data2")),
            tags$br(),
            tags$h4("12.3. Selected ROI Trace", style = "color: gray; margin-top: 10px;"),
            plotOutput(ns("selected_trace_plot"), height = "420px"),
            tags$br(),
            tags$h4("12.4. Selected ROI Descriptive Summary", style = "color: gray; margin-top: 10px;"),
            DTOutput(ns("selected_trace_summary"))
          ),
          tabPanel(
            "13. Peaks",
            tabsetPanel(
              type = "tabs",
              tabPanel(
                "13.1. Metrics",
                DTOutput(ns("table_peaks2")),
                DTOutput(ns("table_peaks22"))
              ),
              tabPanel(
                "13.2. Metric Plots",
                plotOutput(ns("plot_peak3")),
                plotOutput(ns("derivative")),
                #plotOutput(ns("plot_raw_smoothed")),  #### grafica de regresion lineal
                plotOutput(ns("plot_trend_removal")),
                conditionalPanel(
                  condition = "input.remove_trend == true",
                  ns = ns,
                  tags$div(
                    class = "helper-card",
                    h4("13.2.1. Fitted Polynomial Trend Information"),
                    verbatimTextOutput(ns("trend_poly_info")),
                    DTOutput(ns("trend_coef_table"))
                  )
                )
              ),
              tabPanel(
                "13.3. Fourier Analysis",
                fluidRow(
                  column(
                    width = 6,
                    h4("13.3.1. FFT Summary"),
                    verbatimTextOutput(ns("fft_infoText")),
                    plotOutput(ns("fft_fitPlot"), height = "320px"),
                    plotOutput(ns("fft_residPlot"), height = "260px")
                  ),
                  column(
                    width = 6,
                    h4("13.3.2. Spectrum and Frequency Table"),
                    plotOutput(ns("fft_specPlot"), height = "320px"),
                    DTOutput(ns("fft_freqTable"))
                  )
                ),
                fluidRow(
                  column(
                    width = 12,
                    tags$div(
                      class = "helper-card",
                      h4("13.3.3. Explicit Reconstruction Formula Using the Kept FFT Frequencies"),
                      tags$p("This formula writes the recovered series explicitly as a finite sum of cosines and sines using only the frequencies retained by the FFT low-pass filter."),
                      div(class = "formula-box", textOutput(ns("fft_formula_text")))
                    )
                  )
                ),
                fluidRow(
                  column(
                    width = 12,
                    h4("13.3.4. Theory and Interpretation"),
                    uiOutput(ns("fft_formula_box"))
                  )
                )
              )
            )
          ),
          
          tabPanel(
            "14. Baseline Sensitivity Analysis",
            tags$div(
              class = "helper-card",
              h4("14.1. Methodological Note: Why Baseline Matters"),
              tags$p(
                "Many calcium transient analyses transform the fluorescence trace into ΔF/Baseline. ",
                "However, Baseline is not uniquely defined. It may be estimated from an initial segment, ",
                "a minimum-intensity region, a user-defined interval, a fixed value, or a moving percentile baseline. ",
                "Because amplitude, AUC, FWHM, and even the number of accepted events depend on the selected reference level, ",
                "different baseline definitions can produce different quantitative conclusions."
              ),
              tags$p(
                "This panel compares the same calcium trace and the same FFT-smoothed signal under several baseline definitions. ",
                "The purpose is not to choose a universally correct baseline, but to evaluate whether the conclusions are stable ",
                "or highly sensitive to the baseline method."
              ),
              tags$p(
                tags$strong("Interpretation: "),
                "If the values are similar across baseline methods, the analysis is relatively robust. ",
                "If the values change strongly, the interpretation should explicitly report that the results are baseline-sensitive."
              ),
              tags$p(
                tags$strong("Dynamic baseline note: "),
                "For static baselines, the table reports the same scalar value as mean, median, minimum and maximum. ",
                "For the rolling percentile baseline, the baseline changes over time; therefore, the table reports summary statistics of the full time-varying baseline curve."
              )
            ),
            br(),
            DTOutput(ns("baseline_sensitivity_table")),
            br(),
            plotOutput(ns("baseline_sensitivity_plot"), height = "900px"),
            br(),
            uiOutput(ns("baseline_sensitivity_interpretation"))
          ),
          tabPanel(
            "15. FFT Selection Criteria",
            tags$div(
              class = "helper-card",
              h4("15.1. Selection of FFT Low-Frequency Fraction"),
              tags$p("This panel evaluates the candidate fractions f = 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, and 0.40 using four descriptive criteria."),
              tags$p("Specifically, the assessment considers: (1) preservation of the original signal morphology after FFT-based smoothing, (2) the structure and magnitude of the residual signal, (3) the fraction of spectral power retained in the periodogram, and (4) the closeness of each candidate mean peak amplitude (mean_peak) to a weighted reference mean computed across all candidate fractions."),
              tags$p("The weighted reference mean of mean_peak is calculated using the number of detected peaks as weights. Based on Criterion 4, the app proposes a candidate fraction f whose mean_peak is closest to this weighted reference mean. However, the final selection of f remains at the user's discretion."),
              #tags$p(tags$strong("Note:"), "In this panel only, peak detection uses fixed internal exploratory settings: nups = 1, ndowns = 1, min peak distance = 1, FWHP filter = 0, prominence filter = 0, and the minimum peak height is automatically set to the minimum intensity of each FFT-smoothed trace. Therefore, the user inputs under 'Find Peaks Function Arguments' do not affect this FFT selection table.")
            ),
            # tags$div(
            #   class = "helper-card",
            #   h4("Fixed parameters used in this panel"),
            #   uiOutput(ns("fft_grid_current_params"))
            # ),
            #uiOutput(ns("fft_grid_recommendation")),
            DTOutput(ns("fft_grid_summary")),
            fluidRow(
              column(
                width = 12,
                plotOutput(ns("fft_grid_signal_plot"), height = "800px")
              )
            ),
            fluidRow(
              column(
                width = 12,
                plotOutput(ns("fft_grid_resid_plot"), height = "700px")
              )
            ),
            fluidRow(
              column(
                width = 12,
                plotOutput(ns("fft_grid_spec_plot"), height = "700px")
              )
            ),
            uiOutput(ns("fft_grid_metrics_note")),
            fluidRow(
              column(
                width = 12,
                plotOutput(ns("fft_grid_metrics_plot"), height = "520px")
              )
            )
          )
        )
      )
    )
  )
}

# =========================
# Server Module
# =========================
mod_Denoising_data_fft_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    observeEvent(input$param_info_button11, {
      showModal(modalDialog(
        title = "Help",
        size = "l",
        HTML("<p style='text-align:justify;'><strong>This module reproduces the original denoising workflow, but the smoothing step is now based on FFT low-pass filtering instead of LOESS.</strong></p>"),
        HTML("<p style='text-align:justify;'><strong>FFT Low-Frequency Fraction:</strong> controls how many low frequencies are preserved. Smaller values produce stronger smoothing; values near 1 preserve more of the original signal.</p>"),
        HTML("<p style='text-align:justify;'><strong>Keep mean (DC component):</strong> keeps the global average level of the signal during FFT reconstruction.</p>"),
        HTML("<p style='text-align:justify;'><strong>Remove polynomial trend:</strong> when selected, the app fits a polynomial trend to the selected calcium trace using least squares. The residuals, defined as original signal minus fitted trend, become the new calcium trace used for FFT smoothing, peak detection, AUC, FWHM, FWHP, and all downstream metrics.</p>"),
        footer = modalButton("Close")
      ))
    })
    
    filedata <- reactive({
      if (input$data_simulate > 0) {
        fileInput <- example_calcium_data()
        fileInput2 <- NULL
      } else {
        req(input$fileBcsv2)
        ext <- tools::file_ext(input$fileBcsv2$name)
        fileInput1 <- load_file(input$fileBcsv2$name, input$fileBcsv2$datapath, ext)
        
        if (ext %in% c("csv", "tsv")) {
          fileInput <- as.data.frame(fileInput1)
          fileInput2 <- NULL
        } else if (ext == "json") {
          fileInput2 <- fileInput1
          comp <- fileInput2$components
          com <- t(comp)
          time <- seq(0, fileInput2$image_data[2] - 1, by = 1) * fileInput2$image_data[1]
          com <- cbind(time, com)
          fileInput <- as.data.frame(com)
          colnames(fileInput)[1] <- "Time"
        } else {
          stop("Unsupported file format.")
        }
      }
      list(fileInput = fileInput, fileInput2 = fileInput2)
    })
    
    data_info <- reactive({
      req(filedata()$fileInput)
      Nobservations <- nrow(filedata()$fileInput)
      Ncells <- ncol(filedata()$fileInput) - 1
      SummaryData <- data.frame(Number = c(Ncells, Nobservations))
      rownames(SummaryData) <- c("Region of Interest (ROI)", "Time observations")
      list(
        SummaryData = SummaryData,
        data = data.frame(filedata()$fileInput, row.names = NULL)
      )
    })
    
    output$data2 <- renderDT({
      datatable(data_info()$data, options = list(pagingType = "simple"),
                caption = tags$caption(tags$strong("Dataset:")))
    })
    
    output$infotable2 <- renderDT({
      datatable(data_info()$SummaryData, options = list(pagingType = "simple", dom = "t"),
                caption = tags$caption(tags$strong("Dataset Summary:")))
    })
    
    selected_trace <- reactive({
      req(filedata()$fileInput)
      dat <- filedata()$fileInput
      shiny::validate(shiny::need(ncol(dat) >= 2, "Dataset must contain Time plus at least one ROI."))
      
      roi_index <- suppressWarnings(as.integer(input$Cell2))
      if (length(roi_index) == 0 || !is.finite(roi_index) || roi_index < 1) roi_index <- 1L
      roi_index <- min(roi_index, ncol(dat) - 1L)
      roi_col <- roi_index + 1L
      
      data.frame(
        Time = suppressWarnings(as.numeric(dat[[1]])),
        signal = suppressWarnings(as.numeric(dat[[roi_col]]))
      )
    })
    
    output$selected_trace_plot <- renderPlot({
      df <- selected_trace()
      ggplot(df, aes(x = Time, y = signal)) +
        geom_line(linewidth = 0.6) +
        labs(
          title = "Selected ROI Calcium Trace",
          subtitle = paste0("ROI: ", input$Cell2),
          x = "Time",
          y = "Signal"
        ) +
        theme_bw(base_size = 13) +
        theme(
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank()
        )
    })
    
    output$selected_trace_summary <- renderDT({
      df <- selected_trace()
      signal_vals <- suppressWarnings(as.numeric(df$signal))
      time_vals <- suppressWarnings(as.numeric(df$Time))
      
      trace_summary <- data.frame(
        Metric = c(
          "Selected ROI",
          "Time observations",
          "Minimum time",
          "Maximum time",
          "Minimum signal",
          "Maximum signal",
          "Mean signal",
          "Median signal",
          "Standard deviation"
        ),
        Value = c(
          as.character(input$Cell2),
          length(signal_vals),
          round(min(time_vals, na.rm = TRUE), 6),
          round(max(time_vals, na.rm = TRUE), 6),
          round(min(signal_vals, na.rm = TRUE), 6),
          round(max(signal_vals, na.rm = TRUE), 6),
          round(mean(signal_vals, na.rm = TRUE), 6),
          round(stats::median(signal_vals, na.rm = TRUE), 6),
          round(stats::sd(signal_vals, na.rm = TRUE), 6)
        )
      )
      
      datatable(
        trace_summary,
        rownames = FALSE,
        options = list(pagingType = "simple", dom = "t"),
        caption = tags$caption(tags$strong("Selected ROI Summary:"))
      )
    })
    
    peaks_df <- reactive({
      req(filedata()$fileInput)
      data <- filedata()$fileInput
      shiny::validate(shiny::need(ncol(data) >= 2, "Dataset must contain a time column and at least one ROI."))
      
      time <- as.numeric(data[[1]])
      roi_index <- as.numeric(input$Cell2) + 1
      shiny::validate(shiny::need(roi_index <= ncol(data), "Selected ROI is out of range."))
      
      data_original <- data.frame(Time = time, signal = as.numeric(data[[roi_index]]))
      
      trend_res <- NULL
      if (isTRUE(input$remove_trend)) {
        shiny::validate(shiny::need(
          nrow(data_original) >= input$trend_degree + 2,
          "Not enough observations to fit the selected polynomial trend."
        ))
        
        trend_res <- remove_polynomial_trend_trace(
          time = data_original$Time,
          signal = data_original$signal,
          degree = input$trend_degree,
          center_scale_t = isTRUE(input$trend_center_scale_t)
        )
        
        # This becomes the new calcium trace used by the rest of the app.
        data_raw <- trend_res$data
      } else {
        data_raw <- data_original
      }
      
      fft_res <- fft_lowpass_smooth_trace(
        time = data_raw$Time,
        signal = data_raw$signal,
        f = input$fft_fraction,
        keep_mean = isTRUE(input$keep_mean)
      )
      df_smoothed <- fft_res$data
      
      # Separate FFT smoothing used only by Baseline option 6.4
      # (Event-Specific Baseline: Left-Right FFT Minima).
      # This allows the user to tune the low-frequency fraction for baseline
      # minima without changing the main FFT smoothing or peak detection.
      event_baseline_fft_fraction <- safe_fft_fraction(
        input$event_baseline_fft_fraction,
        fallback = 0.20
      )
      
      event_baseline_fft_res <- fft_lowpass_smooth_trace(
        time = data_raw$Time,
        signal = data_raw$signal,
        f = event_baseline_fft_fraction,
        keep_mean = isTRUE(input$keep_mean)
      )
      df_event_baseline_smoothed <- event_baseline_fft_res$data
      
      peaks_found <- peaks(
        data = df_smoothed,
        nups = input$nups2,
        ndowns = input$ndowns2,
        minpeakheight = input$minpeakheight2,
        minpeakdistance = input$minpeakdistance2
      )
      
      list(
        table_peak = normalize_fft_peak_table(peaks_found$p_eak),
        table_positions_peaks = normalize_fft_peak_positions(peaks_found$peak),
        data_original = data_original,
        data_raw = data_raw,
        df_smoothed = df_smoothed,
        df_event_baseline_smoothed = df_event_baseline_smoothed,
        data_matrix = as.matrix(data[, -1, drop = FALSE]),
        fft_res = fft_res,
        event_baseline_fft_res = event_baseline_fft_res,
        event_baseline_fft_fraction = event_baseline_fft_fraction,
        trend_res = trend_res,
        trend_removed = isTRUE(input$remove_trend)
      )
    })
    
    
    # ============================================================
    # Baseline sensitivity analysis
    # ============================================================
    baseline_sensitivity_analysis <- reactive({
      req(peaks_df()$df_smoothed)
      
      data_smoothed <- peaks_df()$df_smoothed
      event_baseline_smoothed <- peaks_df()$df_event_baseline_smoothed
      event_baseline_fft_fraction_num <- peaks_df()$event_baseline_fft_fraction
      
      baseline_methods <- data.frame(
        baseline_mode = c("2", "3", "6", "8", "11", "9", "10", "5"),
        Baseline_Method = c(
          "Standard definition",
          "Low-fluorescence region baseline",
          "Rolling percentile baseline",
          "Event-specific baseline: left-right FFT minima",
          "Constant baseline: mean from trace start to Time_Onset",
          "Median signal (legacy)",
          "Mean signal (legacy)",
          "Min"
        ),
        stringsAsFactors = FALSE
      )
      
      lim_inf_num <- suppressWarnings(as.numeric(input$Lim_inf))
      lim_sup_num <- suppressWarnings(as.numeric(input$Lim_sup))
      own_baseline_num <- suppressWarnings(as.numeric(input$own_baseline))
      initial_fraction_num <- suppressWarnings(as.numeric(input$initial_trace_fraction))
      minimal_fraction_num <- suppressWarnings(as.numeric(input$minimal_lower_fraction))
      rolling_window_num <- suppressWarnings(as.numeric(input$rolling_window))
      rolling_percentile_num <- suppressWarnings(as.numeric(input$rolling_percentile))
      
      if (length(lim_inf_num) == 0 || !is.finite(lim_inf_num[1])) lim_inf_num <- min(data_smoothed$Time, na.rm = TRUE) else lim_inf_num <- lim_inf_num[1]
      if (length(lim_sup_num) == 0 || !is.finite(lim_sup_num[1])) lim_sup_num <- max(data_smoothed$Time, na.rm = TRUE) else lim_sup_num <- lim_sup_num[1]
      if (length(own_baseline_num) == 0 || !is.finite(own_baseline_num[1])) own_baseline_num <- 0 else own_baseline_num <- own_baseline_num[1]
      if (length(initial_fraction_num) == 0 || !is.finite(initial_fraction_num[1])) initial_fraction_num <- 0.10 else initial_fraction_num <- initial_fraction_num[1]
      if (length(minimal_fraction_num) == 0 || !is.finite(minimal_fraction_num[1])) minimal_fraction_num <- 0.10 else minimal_fraction_num <- minimal_fraction_num[1]
      if (length(rolling_window_num) == 0 || !is.finite(rolling_window_num[1])) rolling_window_num <- 20 else rolling_window_num <- rolling_window_num[1]
      if (length(rolling_percentile_num) == 0 || !is.finite(rolling_percentile_num[1])) rolling_percentile_num <- 2 else rolling_percentile_num <- rolling_percentile_num[1]
      
      # Detect peaks once to obtain the time at which each transient starts increasing.
      time_start_increasin_peak <- NULL
      
      peaks_found <- tryCatch({
        peaks(
          data = data_smoothed,
          nups = input$nups2,
          ndowns = input$ndowns2,
          minpeakheight = input$minpeakheight2,
          minpeakdistance = input$minpeakdistance2
        )
      }, error = function(e) NULL)
      
      if (!is.null(peaks_found) &&
          !is.null(peaks_found$peak) &&
          NROW(peaks_found$peak) > 0) {
        
        MSCPFP <- Time_of_the_first_peak(
          data1 = data_smoothed,
          peak = peaks_found$peak
        )$cambios_menor_que_pfp
        
        prom <- prominens2(
          data = data_smoothed,
          peak = peaks_found$peak,
          MSCPFP = MSCPFP
        )
        
        time_start_increasin_peak <- prom$time_start_increasin_peak
      }
      
      event_specific_baseline_info <- NULL
      trace_onset_mean_baseline_info <- NULL
      if (!is.null(peaks_found) &&
          !is.null(peaks_found$peak) &&
          NROW(peaks_found$peak) > 0) {
        event_specific_baseline_info <- compute_event_specific_minima_baseline(
          data_smoothed = data_smoothed,
          table_positions_peaks = peaks_found$peak,
          time_start_increasin_peak = time_start_increasin_peak,
          baseline_smoothed = event_baseline_smoothed,
          baseline_fft_fraction = event_baseline_fft_fraction_num
        )
        
        trace_onset_mean_baseline_info <- compute_trace_start_to_onset_mean_baseline(
          data_smoothed = data_smoothed,
          time_start_increasin_peak = time_start_increasin_peak
        )
      }
      
      sensitivity_list <- lapply(seq_len(nrow(baseline_methods)), function(i) {
        
        mode_i <- baseline_methods$baseline_mode[i]
        method_i <- baseline_methods$Baseline_Method[i]
        
        baseline_info_i <- get_baseline_info(
          data_smoothed = data_smoothed,
          baseline_mode = mode_i,
          lim_inf = lim_inf_num,
          lim_sup = lim_sup_num,
          own_baseline = own_baseline_num,
          initial_fraction = initial_fraction_num,
          minimal_fraction = minimal_fraction_num,
          time_start_increasin_peak = time_start_increasin_peak,
          rolling_window = rolling_window_num,
          rolling_percentile = rolling_percentile_num
        )
        
        if (identical(as.character(mode_i), "8") &&
            !is.null(event_specific_baseline_info) &&
            is.list(event_specific_baseline_info)) {
          baseline_info_i <- event_specific_baseline_info
        }
        
        if (identical(as.character(mode_i), "11") &&
            !is.null(trace_onset_mean_baseline_info) &&
            is.list(trace_onset_mean_baseline_info)) {
          baseline_info_i <- trace_onset_mean_baseline_info
        }
        
        # ------------------------------------------------------------
        # Baseline summary
        # ------------------------------------------------------------
        # Static baselines are represented by one scalar value.
        # Rolling percentile baseline is dynamic/time-varying, so the table
        # reports summary statistics of the full baseline curve.
        if (!is.null(baseline_info_i$trace) && nrow(baseline_info_i$trace) > 0) {
          
          baseline_type_i <- if (identical(as.character(mode_i), "8")) {
            "Event-specific"
          } else if (identical(as.character(mode_i), "11")) {
            "Static / scalar"
          } else {
            "Dynamic / time-varying"
          }
          baseline_values_i <- suppressWarnings(as.numeric(baseline_info_i$trace$Baseline))
          baseline_values_i <- baseline_values_i[is.finite(baseline_values_i)]
          
          baseline_mean_i <- if (length(baseline_values_i) > 0) mean(baseline_values_i, na.rm = TRUE) else NA_real_
          baseline_median_i <- if (length(baseline_values_i) > 0) median(baseline_values_i, na.rm = TRUE) else NA_real_
          baseline_min_i <- if (length(baseline_values_i) > 0) min(baseline_values_i, na.rm = TRUE) else NA_real_
          baseline_max_i <- if (length(baseline_values_i) > 0) max(baseline_values_i, na.rm = TRUE) else NA_real_
          baseline_iqr_i <- if (length(baseline_values_i) > 0) IQR(baseline_values_i, na.rm = TRUE) else NA_real_
          
        } else {
          
          baseline_type_i <- "Static / scalar"
          baseline_values_i <- suppressWarnings(as.numeric(baseline_info_i$scalar))
          if (length(baseline_values_i) == 0 || !is.finite(baseline_values_i[1])) {
            baseline_values_i <- NA_real_
          } else {
            baseline_values_i <- baseline_values_i[1]
          }
          
          baseline_mean_i <- baseline_values_i
          baseline_median_i <- baseline_values_i
          baseline_min_i <- baseline_values_i
          baseline_max_i <- baseline_values_i
          baseline_iqr_i <- 0
        }
        
        metrics_i <- extract_peak_metrics_from_smoothed(
          data_smoothed = data_smoothed,
          nups = input$nups2,
          ndowns = input$ndowns2,
          minpeakheight = input$minpeakheight2,
          minpeakdistance = input$minpeakdistance2,
          baseline_mode = mode_i,
          lim_inf = lim_inf_num,
          lim_sup = lim_sup_num,
          own_baseline = own_baseline_num,
          initial_fraction = initial_fraction_num,
          minimal_fraction = minimal_fraction_num,
          rolling_window = rolling_window_num,
          rolling_percentile = rolling_percentile_num,
          min_FWHP = input$min_FWHP,
          min_prominence = input$min_prominence,
          baseline_info_override = if (as.character(mode_i) %in% c("8", "11")) baseline_info_i else NULL
        )
        
        df_i <- metrics_i$df_p
        
        if (nrow(df_i) == 0) {
          return(data.frame(
            Baseline_Method = method_i,
            Baseline_Type = baseline_type_i,
            Baseline_Mean = baseline_mean_i,
            Baseline_Median = baseline_median_i,
            Baseline_Min = baseline_min_i,
            Baseline_Max = baseline_max_i,
            Baseline_IQR = baseline_iqr_i,
            Number_of_Peaks = 0,
            Frequency = 0,
            Mean_Amplitude_Above_Baseline = NA_real_,
            Mean_Prominence = NA_real_,
            Mean_FWHP = NA_real_,
            Mean_FWHM = NA_real_,
            Mean_Peak_Rise_Time = NA_real_,
            Mean_Peak_Occurence_Time = NA_real_,
            Mean_Transient_Occurence_Time = NA_real_,
            Mean_Onset = NA_real_,
            Mean_Rise_Rate = NA_real_,
            AUC = NA_real_,
            stringsAsFactors = FALSE
          ))
        }
        
        baseline_at_peak_i <- baseline_at_times(
          time_points = df_i$Peak_Occurence_Time,
          baseline_info = baseline_info_i
        )
        
        amplitude_above_baseline_i <- df_i$Amplitude - baseline_at_peak_i
        
        time_range <- max(data_smoothed$Time, na.rm = TRUE) -
          min(data_smoothed$Time, na.rm = TRUE)
        
        frequency_i <- ifelse(
          is.finite(time_range) && time_range > 0,
          nrow(df_i) / time_range,
          NA_real_
        )
        
        auc_i <- tryCatch({
          if (identical(as.character(mode_i), "8")) {
            AUC_event_specific_minima(
              data_smoothed = data_smoothed,
              baseline_info = baseline_info_i
            )$area
          } else if (!is.null(baseline_info_i$trace) && nrow(baseline_info_i$trace) > 0) {
            AUC2_dynamic(
              datos = data_smoothed,
              baseline_trace = baseline_info_i$trace
            )$area
          } else {
            AUC2(
              datos = data_smoothed,
              Integration_Reference = baseline_info_i$scalar
            )$area
          }
        }, error = function(e) NA_real_)
        
        {
          safe_mean_col <- function(df, colname) {
            if (!colname %in% names(df)) return(NA_real_)
            mean(df[[colname]], na.rm = TRUE)
          }
          
          mean_peak_time_i <- safe_mean_col(df_i, "Peak_Occurence_Time")
          mean_onset_i <- safe_mean_col(df_i, "Transient_Ocurrence_Time")
          mean_rise_rate_i <- safe_mean_col(df_i, "Rise_Rate")
          
          data.frame(
            Baseline_Method = method_i,
            Baseline_Type = baseline_type_i,
            Baseline_Mean = baseline_mean_i,
            Baseline_Median = baseline_median_i,
            Baseline_Min = baseline_min_i,
            Baseline_Max = baseline_max_i,
            Baseline_IQR = baseline_iqr_i,
            Number_of_Peaks = nrow(df_i),
            Frequency = frequency_i,
            Mean_Amplitude_Above_Baseline = mean(amplitude_above_baseline_i, na.rm = TRUE),
            Mean_Prominence = mean(df_i$Prominence, na.rm = TRUE),
            Mean_FWHP = mean(df_i$FWHP, na.rm = TRUE),
            Mean_FWHM = mean(df_i$FWHM, na.rm = TRUE),
            Mean_Peak_Rise_Time = mean(df_i$Peak_Rise_Time, na.rm = TRUE),
            Mean_Peak_Occurence_Time = mean_peak_time_i,
            Mean_Transient_Occurence_Time = mean_onset_i,
            Mean_Onset = mean_onset_i,
            Mean_Rise_Rate = mean_rise_rate_i,
            AUC = auc_i,
            stringsAsFactors = FALSE
          )
        }
      })
      
      sensitivity_df <- do.call(rbind, sensitivity_list)
      sensitivity_df
    })
    
    output$baseline_sensitivity_table <- renderDT({
      sensitivity_df <- baseline_sensitivity_analysis()
      
      sensitivity_df_display <- sensitivity_df
      numeric_cols <- vapply(sensitivity_df_display, is.numeric, logical(1))
      sensitivity_df_display[numeric_cols] <- lapply(
        sensitivity_df_display[numeric_cols],
        function(x) round(x, 4)
      )
      
      datatable(
        sensitivity_df_display,
        caption = tags$caption(
          tags$strong("Baseline Sensitivity Analysis: comparison of calcium transient metrics across baseline definitions")
        ),
        options = list(
          pageLength = 10,
          scrollX = TRUE
        )
      )
    })
    
    output$baseline_sensitivity_plot <- renderPlot({
      sensitivity_df <- baseline_sensitivity_analysis()
      
      metric_cols <- c(
        "Number_of_Peaks",
        "Frequency",
        "Mean_Amplitude_Above_Baseline",
        "Mean_Prominence",
        "Mean_FWHP",
        "Mean_FWHM",
        "Mean_Peak_Rise_Time",
        "Mean_Peak_Occurence_Time",
        "Mean_Transient_Occurence_Time",
        "Mean_Onset",
        "Mean_Rise_Rate",
        "AUC"
      )
      
      metric_cols <- metric_cols[metric_cols %in% names(sensitivity_df)]
      
      plot_long <- do.call(
        rbind,
        lapply(metric_cols, function(m) {
          data.frame(
            Baseline_Method = sensitivity_df$Baseline_Method,
            Metric = m,
            Value = sensitivity_df[[m]],
            stringsAsFactors = FALSE
          )
        })
      )
      
      plot_long$Baseline_Method <- factor(
        plot_long$Baseline_Method,
        levels = sensitivity_df$Baseline_Method
      )
      
      ggplot(plot_long, aes(x = Baseline_Method, y = Value)) +
        geom_col() +
        facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
        labs(
          title = "Sensitivity of calcium transient metrics to baseline definition",
          subtitle = "Comparison across all transient-analysis metrics",
          x = "Baseline method",
          y = "Metric value"
        ) +
        theme_classic(base_size = 14) +
        theme(
          axis.text.x = element_text(angle = 35, hjust = 1),
          plot.title = element_text(face = "bold", size = 18),
          strip.text = element_text(face = "bold")
        )
    })
    
    output$baseline_sensitivity_interpretation <- renderUI({
      sensitivity_df <- baseline_sensitivity_analysis()
      
      amp_vals <- sensitivity_df$Mean_Amplitude_Above_Baseline
      auc_vals <- sensitivity_df$AUC
      fwhm_vals <- sensitivity_df$Mean_FWHM
      
      cv_fun <- function(x) {
        x <- x[is.finite(x)]
        if (length(x) < 2 || mean(abs(x), na.rm = TRUE) == 0) return(NA_real_)
        100 * sd(x, na.rm = TRUE) / mean(abs(x), na.rm = TRUE)
      }
      
      cv_amp <- cv_fun(amp_vals)
      cv_auc <- cv_fun(auc_vals)
      cv_fwhm <- cv_fun(fwhm_vals)
      
      cv_vector <- c(cv_amp, cv_auc, cv_fwhm)
      max_cv <- if (all(!is.finite(cv_vector))) NA_real_ else max(cv_vector, na.rm = TRUE)
      
      interpretation <- if (!is.finite(max_cv)) {
        "The sensitivity index could not be computed because some metrics are unavailable."
      } else if (max_cv < 10) {
        "The selected calcium trace appears relatively robust to the baseline definition. The main metrics change only slightly across baseline methods."
      } else if (max_cv < 30) {
        "The selected calcium trace shows moderate baseline sensitivity. The results should report which baseline method was used."
      } else {
        "The selected calcium trace shows strong baseline sensitivity. Biological interpretation should be made carefully because amplitude, AUC, or FWHM change substantially depending on the baseline definition."
      }
      
      tags$div(
        class = "helper-card",
        h4("Automatic interpretation"),
        tags$p(interpretation),
        tags$p(
          tags$strong("Important note for dynamic baselines: "),
          "For static baseline methods, the baseline is represented by a single scalar value. For rolling percentile baseline, the baseline is time-varying; therefore, the table reports summary statistics of the baseline curve, including mean, median, minimum, maximum, and IQR. Metrics such as amplitude above baseline and AUC are computed using the time-varying baseline when applicable."
        ),
        tags$p(
          tags$strong("Sensitivity index used here: "),
          "coefficient of variation across baseline methods for amplitude, AUC, and FWHM."
        ),
        tags$ul(
          tags$li(paste0("Amplitude CV: ", ifelse(is.finite(cv_amp), round(cv_amp, 2), "NA"), "%")),
          tags$li(paste0("AUC CV: ", ifelse(is.finite(cv_auc), round(cv_auc, 2), "NA"), "%")),
          tags$li(paste0("FWHM CV: ", ifelse(is.finite(cv_fwhm), round(cv_fwhm, 2), "NA"), "%"))
        )
      )
    })
    
    fft_candidate_fractions <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)
    fft_grid_fixed_params <- list(
      minpeakheight = 1,
      nups = 1,
      ndowns = 1,
      minpeakdistance = 1,
      min_FWHP = 0,
      min_prominence = 0,
      baseline_mode = "2",
      lim_inf = 0,
      lim_sup = 20,
      own_baseline = 0,
      auc = "1",
      fwhm = "1"
    )
    
    
    fft_panel_data <- reactive({
      req(peaks_df()$data_raw, peaks_df()$fft_res)
      
      df_raw <- peaks_df()$data_raw
      df_smoothed <- peaks_df()$df_smoothed
      res <- peaks_df()$fft_res
      
      df_out <- data.frame(
        Time = df_raw$Time,
        t = seq_len(nrow(df_raw)),
        x_t = df_raw$signal,
        s_hat = df_smoothed$signal,
        resid = df_raw$signal - df_smoothed$signal
      )
      
      n <- nrow(df_out)
      n_pos <- res$n_pos
      k_keep <- res$k_keep
      df_freq <- build_fft_frequency_table(X = res$fft, n = n, k_keep = k_keep)
      
      formula_txt <- build_fft_formula(
        X = res$fft,
        n = n,
        k_keep = k_keep,
        keep_mean = isTRUE(input$keep_mean),
        digits = 3
      )
      
      list(
        df_out = df_out,
        df_freq = df_freq,
        n = n,
        n_pos = n_pos,
        k_keep = k_keep,
        formula_txt = formula_txt
      )
    })
    
    fft_fraction_grid <- reactive({
      req(filedata()$fileInput)
      data <- filedata()$fileInput
      shiny::validate(shiny::need(ncol(data) >= 2, "Dataset must contain a time column and at least one ROI."))
      
      time <- as.numeric(data[[1]])
      roi_index <- as.numeric(input$Cell2) + 1
      shiny::validate(shiny::need(roi_index <= ncol(data), "Selected ROI is out of range."))
      
      df_original <- data.frame(Time = time, signal = as.numeric(data[[roi_index]]))
      
      if (isTRUE(input$remove_trend)) {
        shiny::validate(shiny::need(
          nrow(df_original) >= input$trend_degree + 2,
          "Not enough observations to fit the selected polynomial trend."
        ))
        
        trend_res_grid <- remove_polynomial_trend_trace(
          time = df_original$Time,
          signal = df_original$signal,
          degree = input$trend_degree,
          center_scale_t = isTRUE(input$trend_center_scale_t)
        )
        df_raw <- trend_res_grid$data
      } else {
        df_raw <- df_original
      }
      
      time <- df_raw$Time
      signal <- df_raw$signal
      dt <- if (length(time) > 1) median(diff(time), na.rm = TRUE) else 1
      lim_inf_num <- fft_grid_fixed_params$lim_inf
      lim_sup_num <- fft_grid_fixed_params$lim_sup
      
      candidate_results <- lapply(fft_candidate_fractions, function(f_val) {
        fft_res <- fft_lowpass_smooth_trace(
          time = time,
          signal = signal,
          f = f_val,
          keep_mean = isTRUE(input$keep_mean)
        )
        
        df_smoothed <- fft_res$data
        df_out <- data.frame(
          Time = time,
          t = seq_len(length(time)),
          x_t = signal,
          s_hat = df_smoothed$signal,
          resid = signal - df_smoothed$signal
        )
        
        df_freq <- build_fft_frequency_table(
          X = fft_res$fft,
          n = nrow(df_out),
          k_keep = fft_res$k_keep
        )
        
        # ------------------------------------------------------------
        # Dynamic Peak Height for FFT selection criteria only
        # For each candidate f, use the minimum intensity after FFT
        # as the minpeakheight threshold.
        # ------------------------------------------------------------
        dynamic_minpeakheight <- suppressWarnings(
          min(df_smoothed$signal, na.rm = TRUE)
        )
        
        if (!is.finite(dynamic_minpeakheight)) {
          dynamic_minpeakheight <- 0
        }
        
        # ------------------------------------------------------------
        # Peak detection parameters for FFT selection criteria only
        # ------------------------------------------------------------
        # Important:
        # This panel is intentionally independent from the user controls
        # in "Find Peaks Function Arguments". The goal is to evaluate each
        # FFT candidate f using the same exploratory peak-detection settings.
        #
        # Therefore, the following UI inputs are NOT used here:
        # input$minpeakheight2, input$nups2, input$ndowns2,
        # input$minpeakdistance2, input$min_FWHP, input$min_prominence.
        # ------------------------------------------------------------
        peak_metrics <- extract_peak_metrics_from_smoothed(
          data_smoothed = df_smoothed,
          
          # Fixed internal exploratory settings
          nups = 1,
          ndowns = 1,
          minpeakheight = dynamic_minpeakheight,
          minpeakdistance = 1,
          
          # Fixed baseline settings for FFT selection criteria
          baseline_mode = "2",
          lim_inf = 0,
          lim_sup = 20,
          own_baseline = 0,
          
          # No post-detection filtering in FFT selection criteria
          min_FWHP = 0,
          min_prominence = 0
        )
        
        corr_raw_smooth <- if (sd(signal, na.rm = TRUE) > 0 && sd(df_smoothed$signal, na.rm = TRUE) > 0) {
          cor(signal, df_smoothed$signal, use = "complete.obs")
        } else {
          NA_real_
        }
        
        residual_sd <- sd(df_out$resid, na.rm = TRUE)
        total_power <- sum(df_freq$Power, na.rm = TRUE)
        retained_power_frac <- if (is.finite(total_power) && total_power > 0) {
          sum(df_freq$Power[df_freq$kept], na.rm = TRUE) / total_power
        } else {
          NA_real_
        }
        
        list(
          f = f_val,
          f_label = sprintf("f = %.2f", f_val),
          fft_res = fft_res,
          df_smoothed = df_smoothed,
          df_out = df_out,
          df_freq = df_freq,
          peak_res = peak_metrics,
          minpeakheight_used = dynamic_minpeakheight,
          nups_used = 1,
          ndowns_used = 1,
          minpeakdistance_used = 1,
          min_FWHP_used = 0,
          min_prominence_used = 0,
          corr_raw_smooth = corr_raw_smooth,
          residual_sd = residual_sd,
          retained_power_frac = retained_power_frac,
          omega_cut = 2 * pi * (fft_res$k_keep - 1L) / nrow(df_out)
        )
      })
      
      ref_idx <- which.max(sapply(candidate_results, function(x) x$f))
      ref_res <- candidate_results[[ref_idx]]
      ref_residual_sd <- ref_res$residual_sd
      
      summary_df <- do.call(rbind, lapply(candidate_results, function(res) {
        summ <- res$peak_res$summary
        
        data.frame(
          f = res$f,
          f_label = res$f_label,
          k_keep = res$fft_res$k_keep,
          corr_raw_smooth = res$corr_raw_smooth,
          residual_sd = res$residual_sd,
          retained_power_frac = res$retained_power_frac,
          minpeakheight_used = res$minpeakheight_used,
          nups_used = res$nups_used,
          ndowns_used = res$ndowns_used,
          minpeakdistance_used = res$minpeakdistance_used,
          min_FWHP_used = res$min_FWHP_used,
          min_prominence_used = res$min_prominence_used,
          n_peaks = summ$n_peaks,
          mean_peak = summ$mean_peak,
          stringsAsFactors = FALSE
        )
      }))
      
      summary_df <- summary_df[order(summary_df$f), , drop = FALSE]
      
      pass_signal_vec <- is.finite(summary_df$corr_raw_smooth) & summary_df$corr_raw_smooth >= 0.90
      pass_residual_vec <- if (is.finite(ref_residual_sd) && ref_residual_sd > 0) {
        summary_df$residual_sd <= 1.5 * ref_residual_sd
      } else {
        rep(TRUE, nrow(summary_df))
      }
      pass_spectrum_vec <- is.finite(summary_df$retained_power_frac) & summary_df$retained_power_frac >= 0.90
      
      valid_mean_peak <- is.finite(summary_df$mean_peak) & is.finite(summary_df$n_peaks) & summary_df$n_peaks > 0
      if (any(valid_mean_peak)) {
        weighted_mean_peak <- stats::weighted.mean(
          x = summary_df$mean_peak[valid_mean_peak],
          w = summary_df$n_peaks[valid_mean_peak],
          na.rm = TRUE
        )
        mean_peak_abs_diff <- abs(summary_df$mean_peak - weighted_mean_peak)
        min_diff <- min(mean_peak_abs_diff[valid_mean_peak], na.rm = TRUE)
        pass_mean_peak_vec <- rep(FALSE, nrow(summary_df))
        pass_mean_peak_vec[valid_mean_peak] <- abs(mean_peak_abs_diff[valid_mean_peak] - min_diff) < 1e-12
      } else {
        weighted_mean_peak <- NA_real_
        mean_peak_abs_diff <- rep(NA_real_, nrow(summary_df))
        pass_mean_peak_vec <- rep(FALSE, nrow(summary_df))
      }
      
      summary_df$pass_signal <- pass_signal_vec
      summary_df$pass_residual <- pass_residual_vec
      summary_df$pass_spectrum <- pass_spectrum_vec
      summary_df$weighted_mean_peak <- weighted_mean_peak
      summary_df$mean_peak_abs_diff <- mean_peak_abs_diff
      summary_df$pass_mean_peak <- pass_mean_peak_vec
      summary_df$overall_pass <- pass_signal_vec & pass_residual_vec & pass_spectrum_vec & pass_mean_peak_vec
      summary_df$score_4criteria <- rowSums(cbind(
        pass_signal_vec,
        pass_residual_vec,
        pass_spectrum_vec,
        pass_mean_peak_vec
      ), na.rm = TRUE)
      
      if (any(summary_df$overall_pass, na.rm = TRUE)) {
        recommended_idx <- which(summary_df$overall_pass)[1]
        recommendation_mode <- "strict"
        recommendation_text <- paste0(
          "Recommended value: ", sprintf("%.2f", summary_df$f[recommended_idx]),
          ". This candidate has a mean_peak value closest to the weighted reference mean computed across all candidate fractions. This recommendation is based on Criterion 4 and should be interpreted as a quantitative guide; the final choice of f remains at the user's discretion."
        )
      } else {
        max_score <- max(summary_df$score_4criteria, na.rm = TRUE)
        recommended_idx <- which(summary_df$score_4criteria == max_score)[1]
        recommendation_mode <- "compromise"
        recommendation_text <- paste0(
          "Best compromise: ", sprintf("%.2f", summary_df$f[recommended_idx]),
          ". No candidate satisfied the four criteria simultaneously, so the app returns the smallest value with the highest criterion score."
        )
      }
      
      summary_df$recommended <- FALSE
      summary_df$recommended[recommended_idx] <- TRUE
      summary_df$recommendation_mode <- recommendation_mode
      
      overlay_df <- do.call(rbind, lapply(candidate_results, function(res) {
        rbind(
          data.frame(Time = res$df_out$Time, value = res$df_out$x_t, series = "Raw", f = res$f, f_label = res$f_label),
          data.frame(Time = res$df_out$Time, value = res$df_out$s_hat, series = "Smoothed", f = res$f, f_label = res$f_label)
        )
      }))
      
      residual_df <- do.call(rbind, lapply(candidate_results, function(res) {
        data.frame(Time = res$df_out$Time, resid = res$df_out$resid, f = res$f, f_label = res$f_label)
      }))
      
      periodogram_df <- do.call(rbind, lapply(candidate_results, function(res) {
        data.frame(
          omega = res$df_freq$omega,
          I = res$df_freq$I,
          kept = res$df_freq$kept,
          f = res$f,
          f_label = res$f_label,
          omega_cut = res$omega_cut
        )
      }))
      
      cutoff_df <- unique(periodogram_df[, c("f", "f_label", "omega_cut")])
      
      mean_peak_df <- data.frame(
        f = summary_df$f,
        f_label = summary_df$f_label,
        mean_peak = summary_df$mean_peak,
        weighted_mean_peak = summary_df$weighted_mean_peak,
        recommended = summary_df$recommended,
        pass_mean_peak = summary_df$pass_mean_peak
      )
      
      criterion4_note <- if (!is.finite(weighted_mean_peak)) {
        paste0(
          "Criterion 4 could not be computed because no candidate produced valid peak detections under the current thresholds. ",
          "As a result, the weighted mean of mean_peak is unavailable."
        )
      } else {
        paste0(
          "Criterion 4 uses the weighted mean of mean_peak, with the number of detected peaks as weights. ",
          "Weighted mean_peak = ", sprintf("%.4f", weighted_mean_peak), "."
        )
      }
      
      list(
        summary_df = summary_df,
        overlay_df = overlay_df,
        residual_df = residual_df,
        periodogram_df = periodogram_df,
        cutoff_df = cutoff_df,
        mean_peak_df = mean_peak_df,
        recommendation_text = recommendation_text,
        criterion4_note = criterion4_note,
        weighted_mean_peak = weighted_mean_peak,
        recommended_f = summary_df$f[recommended_idx],
        recommended_label = summary_df$f_label[recommended_idx]
      )
    })
    
    peaks_plot <- reactive({
      req(nrow(peaks_df()$table_peak) > 0)
      
      table_peak <- peaks_df()$table_peak
      table_positions_peaks <- peaks_df()$table_positions_peaks
      data_raw <- peaks_df()$data_raw
      data_smoothed <- peaks_df()$df_smoothed
      table_positions_peaks <- normalize_fft_peak_positions(table_positions_peaks)
      if (nrow(table_positions_peaks) == 0) return(NULL)
      peaks_idx <- table_positions_peaks$peak_index
      
      MSCPFP <- Time_of_the_first_peak(data1 = data_smoothed, peak = table_positions_peaks)$cambios_menor_que_pfp
      prom <- prominens2(data = data_smoothed, peak = table_positions_peaks, MSCPFP = MSCPFP)
      data_min <- prom$data_min
      df_peaks_parcia <- prom$df_peaks_parcia
      time_start_increasin_peak <- prom$time_start_increasin_peak
      
      Puntos_medios <- FWHP2(peaks = data_smoothed[, 1][peaks_idx], df_peaks_parcia = df_peaks_parcia)$Puntos_medios
      table_peak$prominence <- prom$prominens_amplitud
      table_peak$Prominence_Midpoint <- Puntos_medios$p_eak_mediun
      
      right_left <- right_left_FWHP(data1 = data_smoothed, peak = table_positions_peaks, P_M = Puntos_medios)
      left_FWHP <- right_left$df
      right_FWHP <- right_left$df2
      
      table_peak$Time_left_FWHP <- left_FWHP$Time_left_FWHP
      table_peak$Time_right_FWHP <- right_FWHP$Time_right_FWHP
      table_peak$FWHP <- right_FWHP$Time_right_FWHP - left_FWHP$Time_left_FWHP
      table_peak$Time_to_peak <- table_peak$posision_peaks - time_start_increasin_peak$Time
      table_peak$puntominimo_y <- prom$df_peaks_parcia$p_fin1
      
      baseline_info <- get_baseline_info(
        data_smoothed = data_smoothed,
        baseline_mode = input$Baseline,
        lim_inf = input$Lim_inf,
        lim_sup = input$Lim_sup,
        own_baseline = input$own_baseline,
        initial_fraction = input$initial_trace_fraction,
        minimal_fraction = input$minimal_lower_fraction,
        time_start_increasin_peak = time_start_increasin_peak,
        rolling_window = input$rolling_window,
        rolling_percentile = input$rolling_percentile
      )
      
      # Event-specific baseline: left-right FFT minima.
      # Peak detection/AUC use the main FFT-smoothed trace. The left/right
      # minima for this baseline are found on a separate FFT-smoothed trace
      # controlled by 6.4.1.
      if (identical(as.character(input$Baseline), "8")) {
        baseline_info <- compute_event_specific_minima_baseline(
          data_smoothed = data_smoothed,
          table_positions_peaks = table_positions_peaks,
          time_start_increasin_peak = time_start_increasin_peak,
          baseline_smoothed = peaks_df()$df_event_baseline_smoothed,
          baseline_fft_fraction = peaks_df()$event_baseline_fft_fraction
        )
      }
      
      if (identical(as.character(input$Baseline), "11")) {
        baseline_info <- compute_trace_start_to_onset_mean_baseline(
          data_smoothed = data_smoothed,
          time_start_increasin_peak = time_start_increasin_peak
        )
      }
      
      baseline1 <- baseline_info$scalar
      baseline_at_peak <- baseline_at_times(table_peak$posision_peaks, baseline_info)
      table_peak$Baseline_at_peak <- baseline_at_peak
      
      if (identical(as.character(input$auc2), "2")) {
        if (identical(as.character(input$Baseline), "8")) {
          AUC <- AUC_event_specific_minima(data_smoothed = data_smoothed, baseline_info = baseline_info)
        } else if (!is.null(baseline_info$trace) && nrow(baseline_info$trace) > 0) {
          AUC <- AUC2_dynamic(datos = data_smoothed, baseline_trace = baseline_info$trace)
        } else {
          AUC <- AUC2(datos = data_smoothed, Integration_Reference = baseline1)
        }
        tabla_AUC <- data.frame(AUC = AUC$area, P_min = AUC$P_min, P_max = AUC$P_max)
      } else {
        tabla_AUC <- data.frame()
      }
      
      df_raw_smoothed <- data.frame(data_smoothed = data_smoothed[, 2], data_raw = data_raw[, 2])
      modelo <- lm(data_smoothed ~ data_raw, data = df_raw_smoothed)
      r_cuadrado <- summary(modelo)$r.squared
      intercept <- coef(modelo)[1]
      slope_fit <- coef(modelo)[2]
      equation_text <- sprintf("y = %.2fx + %.2f", slope_fit, intercept)
      
      gg2 <- ggplot(df_raw_smoothed, aes(x = data_raw, y = data_smoothed)) +
        geom_point() +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
        labs(title = "Linear Regression",
             x = "Raw [Delta F/Baseline]",
             y = "Smoothed [Delta F/Baseline]") +
        theme_classic() +
        theme(
          plot.title = element_text(size = 28, face = "bold"),
          axis.title = element_text(size = 28, face = "bold"),
          axis.text = element_text(size = 16, face = "bold")
        ) +
        geom_text(
          x = mean(df_raw_smoothed$data_raw) - sd(df_raw_smoothed$data_raw),
          y = mean(df_raw_smoothed$data_smoothed) + sd(df_raw_smoothed$data_raw),
          label = equation_text, hjust = 0, vjust = 0, size = 5
        ) +
        geom_text(
          x = mean(df_raw_smoothed$data_raw) - sd(df_raw_smoothed$data_raw),
          y = mean(df_raw_smoothed$data_smoothed) + 2 * sd(df_raw_smoothed$data_raw),
          label = paste("R^2 =", round(r_cuadrado, 4)), hjust = 0, vjust = 0, size = 5
        )
      
      table_peak$Transient_Ocurrence_Time <- time_start_increasin_peak$Time
      
      deri1_vals <- prospectr::savitzkyGolay(
        X = data_smoothed$signal,
        m = 1,
        p = 2,
        w = 5
      )
      
      n_der <- length(deri1_vals)
      n_time <- nrow(data_smoothed)
      start_idx <- floor((n_time - n_der) / 2) + 1
      end_idx <- start_idx + n_der - 1
      
      primera_derivada1 <- data.frame(
        Time = data_smoothed$Time[start_idx:end_idx],
        deri1 = deri1_vals
      )
      
      data_minimos_crecientes <- data.frame(
        x1 = time_start_increasin_peak$Time,
        y1 = data_min$y,
        x2 = table_peak$posision_peaks,
        y2 = table_peak$absolute_amplitude
      )
      
      slope <- numeric(nrow(data_minimos_crecientes))
      for (i in seq_len(nrow(data_minimos_crecientes))) {
        resultados_filtrados <- primera_derivada1[
          primera_derivada1$Time >= data_minimos_crecientes$x1[i] &
            primera_derivada1$Time <= data_minimos_crecientes$x2[i],
          ,
          drop = FALSE
        ]
        slope[i] <- if (nrow(resultados_filtrados) > 0 && any(is.finite(resultados_filtrados$deri1))) {
          max(resultados_filtrados$deri1, na.rm = TRUE)
        } else {
          NA_real_
        }
      }
      table_peak$slope <- slope
      
      list(
        gg2 = gg2,
        table_peak = table_peak,
        tabla_AUC = tabla_AUC,
        baseline1 = baseline1,
        baseline_info = baseline_info,
        baseline_trace = baseline_info$trace,
        primera_derivada1 = primera_derivada1
      )
    })
    
    output$derivative <- renderPlot({
      data_derivative <- peaks_plot()$primera_derivada1
      ggplot(data_derivative, aes(x = Time, y = deri1)) +
        geom_line(linetype = "solid", linewidth = 1.5, color = "black") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "purple") +
        labs(title = "First Derivative",
             x = "Time [s]",
             y = "Delta F/Baseline * s^-1") +
        theme_classic() +
        theme(
          plot.title = element_text(size = 28, face = "bold"),
          axis.title.y = element_text(size = 28, face = "bold"),
          axis.title.x = element_text(size = 28, face = "bold"),
          axis.text.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold")
        )
    })
    
    peaks_FWHM <- reactive({
      table_positions_peaks <- peaks_df()$table_positions_peaks
      data_smoothed <- peaks_df()$df_smoothed
      baseline_info <- peaks_plot()$baseline_info
      
      if (is.null(table_positions_peaks) || NROW(table_positions_peaks) == 0) {
        return(list(df_FWHM = data.frame(), FWHM = numeric(0)))
      }
      
      table_positions_peaks <- normalize_fft_peak_positions(table_positions_peaks)
      if (nrow(table_positions_peaks) == 0) {
        return(list(df_FWHM = data.frame(), FWHM = numeric(0)))
      }
      peaks_idx <- table_positions_peaks$peak_index
      peak_times <- data_smoothed[, 1][peaks_idx]
      baseline_peak <- baseline_at_times(peak_times, baseline_info)
      p_eak_mediun <- c((table_positions_peaks[, 1] + baseline_peak) / 2)
      
      Puntos_medios <- data.frame(
        posiscion_medio = peak_times,
        p_eak_mediun = p_eak_mediun
      )
      
      right_left_FWHM <- right_left_FWHP(data1 = data_smoothed, peak = table_positions_peaks, P_M = Puntos_medios)
      left_FWHM <- right_left_FWHM$df
      right_FWHM <- right_left_FWHM$df2
      
      len_left <- length(left_FWHM$Time_left_FWHP)
      len_right <- length(right_FWHM$Time_right_FWHP)
      len_mid <- length(p_eak_mediun)
      min_len <- min(len_left, len_right, len_mid)
      
      if (!is.finite(min_len) || min_len == 0) {
        return(list(df_FWHM = data.frame(), FWHM = numeric(0)))
      }
      
      df_FWHM <- data.frame(
        Time_left_FWHM = left_FWHM$Time_left_FWHP[seq_len(min_len)],
        Time_right_FWHM = right_FWHM$Time_right_FWHP[seq_len(min_len)],
        Amplitude_Midpoint = p_eak_mediun[seq_len(min_len)]
      )
      FWHM <- df_FWHM$Time_right_FWHM - df_FWHM$Time_left_FWHM
      list(df_FWHM = df_FWHM, FWHM = FWHM)
    })
    
    Peaks_Data_Final <- reactive({
      df_p <- peaks_plot()$table_peak
      fwhm_res <- peaks_FWHM()
      df_FWHM1 <- fwhm_res$df_FWHM
      
      if (nrow(df_p) == 0 || nrow(df_FWHM1) == 0 || length(fwhm_res$FWHM) == 0) {
        return(list(df_p = data.frame()))
      }
      
      min_len <- min(nrow(df_p), nrow(df_FWHM1), length(fwhm_res$FWHM))
      df_p <- df_p[seq_len(min_len), , drop = FALSE]
      df_FWHM1 <- df_FWHM1[seq_len(min_len), , drop = FALSE]
      df_p$FWHM <- fwhm_res$FWHM[seq_len(min_len)]
      
      colnames(df_p) <- c(
        "Amplitude", "Peak_Occurence_Time", "L_inf", "L_sup",
        "Prominence", "Prominence_Midpoint", "Time_left_FWHP",
        "Time_right_FWHP", "FWHP", "Peak_Rise_Time",
        "puntominimo_y", "Baseline_at_peak", "Transient_Ocurrence_Time",
        "Rise_Rate", "FWHM"
      )
      
      df_FWHM2 <- cbind(df_p, df_FWHM1)
      
      # Stable ID for every original detected peak candidate.
      # compute_event_specific_minima_baseline() uses Event_ID = original
      # peak-candidate order. This ID is created before any filtering so
      # accepted peaks, sigmoid curves, event-specific baselines, and AUC
      # ribbons can be synchronized after filters remove invalid candidates.
      df_FWHM2$.Peak_Candidate_ID <- seq_len(nrow(df_FWHM2))
      
      # ------------------------------------------------------------
      # Freeze width metrics before optional sigmoid modeling.
      # These frozen columns are the source of truth for FWHP and FWHM.
      # Single sigmoid and double sigmoid can add their own metrics, but
      # they do not change these original width measurements.
      # ------------------------------------------------------------
      df_FWHM2$FWHP_NoSigmoid <- df_FWHM2$FWHP
      df_FWHM2$Time_left_FWHP_NoSigmoid <- df_FWHM2$Time_left_FWHP
      df_FWHM2$Time_right_FWHP_NoSigmoid <- df_FWHM2$Time_right_FWHP
      df_FWHM2$FWHP_Level_NoSigmoid <- df_FWHM2$Prominence_Midpoint
      
      df_FWHM2$FWHM_NoSigmoid <- df_FWHM2$FWHM
      df_FWHM2$Time_left_FWHM_NoSigmoid <- df_FWHM2$Time_left_FWHM
      df_FWHM2$Time_right_FWHM_NoSigmoid <- df_FWHM2$Time_right_FWHM
      df_FWHM2$FWHM_Level_NoSigmoid <- df_FWHM2$Amplitude_Midpoint
      
      # Add event-specific baseline diagnostic columns when available.
      baseline_info_for_table <- peaks_plot()$baseline_info
      if (!is.null(baseline_info_for_table$event_baselines) &&
          is.data.frame(baseline_info_for_table$event_baselines) &&
          nrow(baseline_info_for_table$event_baselines) > 0) {
        
        eb <- baseline_info_for_table$event_baselines
        
        # Match event-specific baseline objects back to final peak rows by
        # peak time, not by row number. Some candidate peaks can be removed
        # later by FWHP/prominence/baseline filters, so row-number matching can
        # leave green baseline segments on transients that do not have a red
        # accepted peak.
        dt_tol <- suppressWarnings(stats::median(diff(sort(unique(peaks_df()$df_smoothed$Time))), na.rm = TRUE))
        if (length(dt_tol) == 0 || !is.finite(dt_tol) || dt_tol <= 0) dt_tol <- 1e-8
        peak_match_tol <- max(dt_tol / 2, 1e-8)
        
        df_FWHM2$Baseline_Event_ID <- NA_integer_
        df_FWHM2$Baseline_Left_Min_Time <- NA_real_
        df_FWHM2$Baseline_Left_Min_Value <- NA_real_
        df_FWHM2$Baseline_Right_Min_Time <- NA_real_
        df_FWHM2$Baseline_Right_Min_Value <- NA_real_
        df_FWHM2$Baseline_AUC_Left_Time <- NA_real_
        df_FWHM2$Baseline_AUC_Right_Time <- NA_real_
        
        for (ii in seq_len(nrow(df_FWHM2))) {
          dd <- abs(eb$Peak_Time - df_FWHM2$Peak_Occurence_Time[ii])
          jj <- which.min(dd)
          
          if (length(jj) > 0 && is.finite(jj) && is.finite(dd[jj]) && dd[jj] <= peak_match_tol) {
            df_FWHM2$Baseline_Event_ID[ii] <- eb$Event_ID[jj]
            df_FWHM2$Baseline_Left_Min_Time[ii] <- eb$Left_Min_Time[jj]
            df_FWHM2$Baseline_Left_Min_Value[ii] <- eb$Left_Min_Value[jj]
            df_FWHM2$Baseline_Right_Min_Time[ii] <- eb$Right_Min_Time[jj]
            df_FWHM2$Baseline_Right_Min_Value[ii] <- eb$Right_Min_Value[jj]
            df_FWHM2$Baseline_AUC_Left_Time[ii] <- eb$Plot_Left_Time[jj]
            df_FWHM2$Baseline_AUC_Right_Time[ii] <- eb$Plot_Right_Time[jj]
          }
        }
        
        # In event-specific baseline mode, keep only peaks that truly have an
        # event-specific baseline object. If a candidate peak did not produce
        # a valid event baseline, it should not have a red peak, green baseline,
        # AUC ribbon, single sigmoid, or double sigmoid.
        if (identical(as.character(input$Baseline), "8")) {
          df_FWHM2 <- df_FWHM2[is.finite(df_FWHM2$Baseline_Event_ID), , drop = FALSE]
        }
      }
      
      df_FWHM2 <- df_FWHM2[df_FWHM2$FWHP > input$min_FWHP, , drop = FALSE]
      
      # ------------------------------------------------------------
      # Critical pre-sigmoid baseline filter
      # ------------------------------------------------------------
      # A local maximum that is at or below the selected baseline is not a
      # valid positive calcium transient. Therefore it must be removed before
      # optional single/double sigmoid modeling. This prevents the app from
      # drawing fitted sigmoid curves on sub-baseline bumps.
      if ("Baseline_at_peak" %in% names(df_FWHM2)) {
        df_FWHM2 <- df_FWHM2[
          is.finite(df_FWHM2$Baseline_at_peak) &
            is.finite(df_FWHM2$Amplitude) &
            df_FWHM2$Amplitude > df_FWHM2$Baseline_at_peak,
          ,
          drop = FALSE
        ]
      } else if (is.finite(peaks_plot()$baseline1)) {
        df_FWHM2 <- df_FWHM2[
          is.finite(df_FWHM2$Amplitude) & df_FWHM2$Amplitude > peaks_plot()$baseline1,
          ,
          drop = FALSE
        ]
      }
      
      # If double sigmoid is selected, prominence may be recalculated after
      # the double-sigmoid fit. The event-specific baseline itself remains
      # defined by the original calcium-transient left/right FFT minima.
      # For non-double-sigmoid runs, keep the original prominence filter here.
      if (!isTRUE(input$use_double_sigmoid_fft)) {
        df_FWHM2 <- df_FWHM2[df_FWHM2$Prominence > input$min_prominence, , drop = FALSE]
      }
      
      # .Peak_Candidate_ID was created before filtering and is kept through
      # sigmoid fitting. The fitted curves use this same stable ID.
      
      # Optional sigmoid metrics.
      sigmoid_out <- append_fft_sigmoid_metrics(
        metrics = df_FWHM2,
        time = peaks_df()$df_smoothed$Time,
        signal = peaks_df()$df_smoothed$signal,
        use_single_sigmoid = isTRUE(input$use_single_sigmoid_fft),
        use_double_sigmoid = isTRUE(input$use_double_sigmoid_fft),
        pre_event_window = input$sigmoid_pre_event_window_fft,
        single_left_pad = input$single_sigmoid_left_pad_fft,
        single_right_pad = input$single_sigmoid_right_pad_fft,
        double_left_pad = input$double_sigmoid_left_pad_fft,
        double_right_pad = input$double_sigmoid_right_pad_fft
      )
      
      metrics_final <- sigmoid_out$metrics
      sigmoid_baseline_info <- NULL
      
      # ============================================================
      # User-requested behavior for double sigmoid
      # ============================================================
      # If "10.2. Apply Double Sigmoid Fit for Full-Transient Metrics" is selected:
      #   1) Prominence is recalculated using the minimum of the fitted
      #      double-sigmoid curve on the left side of the detected peak.
      #   2) The event-specific baseline remains unchanged. When the selected
      #      baseline is "Event-Specific Baseline: Left-Right FFT Minima", the
      #      baseline line, AUC ribbon, Baseline_at_peak, and left/right minima
      #      continue to use the original calcium-transient FFT minima stored
      #      in peaks_plot()$baseline_info. Double-sigmoid minima are kept only
      #      in DoubleSigmoid_* diagnostic columns.
      # ============================================================
      if (isTRUE(input$use_double_sigmoid_fft) &&
          !is.null(metrics_final) &&
          nrow(metrics_final) > 0 &&
          "DoubleSigmoid_Fit_Status" %in% names(metrics_final)) {
        
        valid_double <- metrics_final$DoubleSigmoid_Fit_Status == "OK" &
          is.finite(metrics_final$DoubleSigmoid_Left_Min_Time) &
          is.finite(metrics_final$DoubleSigmoid_Left_Min_Value) &
          is.finite(metrics_final$DoubleSigmoid_Right_Min_Time) &
          is.finite(metrics_final$DoubleSigmoid_Right_Min_Value)
        
        if (any(valid_double)) {
          # Preserve original values for comparison.
          metrics_final$Original_Prominence <- metrics_final$Prominence
          metrics_final$Original_Prominence_Midpoint <- metrics_final$Prominence_Midpoint
          metrics_final$Original_puntominimo_y <- metrics_final$puntominimo_y
          if ("Baseline_at_peak" %in% names(metrics_final)) {
            metrics_final$Original_Baseline_at_peak <- metrics_final$Baseline_at_peak
          }
          if ("Baseline_Left_Min_Time" %in% names(metrics_final)) {
            metrics_final$Original_Baseline_Left_Min_Time <- metrics_final$Baseline_Left_Min_Time
            metrics_final$Original_Baseline_Left_Min_Value <- metrics_final$Baseline_Left_Min_Value
          }
          if ("Baseline_Right_Min_Time" %in% names(metrics_final)) {
            metrics_final$Original_Baseline_Right_Min_Time <- metrics_final$Baseline_Right_Min_Time
            metrics_final$Original_Baseline_Right_Min_Value <- metrics_final$Baseline_Right_Min_Value
          }
          
          # Prominence based on the left-side minimum of the double-sigmoid fit.
          # This does not alter the event-specific baseline line or AUC baseline.
          metrics_final$Prominence[valid_double] <-
            metrics_final$Amplitude[valid_double] -
            metrics_final$DoubleSigmoid_Left_Min_Value[valid_double]
          
          metrics_final$Prominence_Midpoint[valid_double] <-
            metrics_final$DoubleSigmoid_Left_Min_Value[valid_double] +
            metrics_final$Prominence[valid_double] / 2
          
          metrics_final$puntominimo_y[valid_double] <-
            metrics_final$DoubleSigmoid_Left_Min_Value[valid_double]
          
          metrics_final$Prominence_Method <- ifelse(
            valid_double,
            "Double sigmoid left minimum",
            "Original prominence"
          )
          
          if (identical(as.character(input$Baseline), "8")) {
            metrics_final$Baseline_Method_For_Event <-
              "Original calcium-transient left-right FFT minima; not double-sigmoid minima"
          }
        }
        
        # Re-apply filters after the double-sigmoid prominence update.
        # Baseline_at_peak remains the original event-specific calcium-transient
        # baseline when baseline mode 8 is selected.
        metrics_final <- metrics_final[
          is.finite(metrics_final$Prominence) & metrics_final$Prominence > input$min_prominence,
          ,
          drop = FALSE
        ]
        if ("Baseline_at_peak" %in% names(metrics_final)) {
          metrics_final <- metrics_final[
            is.finite(metrics_final$Baseline_at_peak) & metrics_final$Amplitude > metrics_final$Baseline_at_peak,
            ,
            drop = FALSE
          ]
        }
      }
      
      # Final safety restoration: FWHP and FWHM stay locked to the
      # original non-sigmoid geometry even when single/double sigmoid options
      # are active and even after post-sigmoid filtering.
      metrics_final <- restore_sigmoid_independent_width_metrics(metrics_final)
      
      # ------------------------------------------------------------
      # Critical post-filter for plotted sigmoid curves
      # ------------------------------------------------------------
      # append_fft_sigmoid_metrics() creates fitted curves before the final
      # post-sigmoid filters are applied. Therefore, an event can disappear
      # from the final metrics table/red peak points but its blue/purple curve
      # can remain in sigmoid_out$fitted_curves. Keep only curves whose
      # Event_ID is still present in the final accepted metrics.
      sigmoid_curves_final <- sigmoid_out$fitted_curves
      if (!is.null(sigmoid_curves_final) &&
          is.data.frame(sigmoid_curves_final) &&
          nrow(sigmoid_curves_final) > 0 &&
          "Event_ID" %in% names(sigmoid_curves_final) &&
          ".Peak_Candidate_ID" %in% names(metrics_final)) {
        keep_event_ids <- suppressWarnings(as.integer(metrics_final$.Peak_Candidate_ID))
        keep_event_ids <- keep_event_ids[is.finite(keep_event_ids)]
        sigmoid_curves_final <- sigmoid_curves_final[
          sigmoid_curves_final$Event_ID %in% keep_event_ids,
          ,
          drop = FALSE
        ]
      } else if (is.null(metrics_final) || nrow(metrics_final) == 0) {
        sigmoid_curves_final <- data.frame()
      }
      
      list(
        df_p = metrics_final,
        sigmoid_curves = sigmoid_curves_final,
        sigmoid_baseline_info = sigmoid_baseline_info
      )
    })
    
    Trance_Graph <- reactive({
      data_smoothed <- peaks_df()$df_smoothed
      data_raw <- peaks_df()$data_raw
      colnames(data_smoothed) <- c("Time", "Sing")
      df_p <- Peaks_Data_Final()$df_p
      
      # ------------------------------------------------------------
      # Accepted events only
      # ------------------------------------------------------------
      # peaks_plot()$baseline_info is computed before the final filters used
      # in Peaks_Data_Final(). Therefore, event-specific baseline segments/AUC
      # ribbons may still contain candidate events that were later removed from
      # the final metrics table. Keep the green baseline/AUC synchronized with
      # the final red accepted peaks.
      #
      # Green baseline objects use baseline_info$event_baselines$Event_ID.
      # Sigmoid curves use .Peak_Candidate_ID. Do not mix these two IDs.
      accepted_baseline_event_ids <- integer(0)
      if (!is.null(df_p) && nrow(df_p) > 0 && "Baseline_Event_ID" %in% names(df_p)) {
        accepted_baseline_event_ids <- suppressWarnings(as.integer(df_p$Baseline_Event_ID))
        accepted_baseline_event_ids <- accepted_baseline_event_ids[is.finite(accepted_baseline_event_ids)]
      }
      
      # Fallback for older cached objects: recover accepted baseline Event_IDs
      # by matching final peak times to event_baselines$Peak_Time.
      if (length(accepted_baseline_event_ids) == 0 &&
          !is.null(df_p) && nrow(df_p) > 0 &&
          "Peak_Occurence_Time" %in% names(df_p) &&
          !is.null(peaks_plot()$baseline_info$event_baselines) &&
          is.data.frame(peaks_plot()$baseline_info$event_baselines) &&
          nrow(peaks_plot()$baseline_info$event_baselines) > 0) {
        eb_tmp <- peaks_plot()$baseline_info$event_baselines
        dt_tol <- suppressWarnings(stats::median(diff(sort(unique(data_smoothed$Time))), na.rm = TRUE))
        if (length(dt_tol) == 0 || !is.finite(dt_tol) || dt_tol <= 0) dt_tol <- 1e-8
        peak_match_tol <- max(dt_tol / 2, 1e-8)
        for (pt in suppressWarnings(as.numeric(df_p$Peak_Occurence_Time))) {
          dd <- abs(eb_tmp$Peak_Time - pt)
          jj <- which.min(dd)
          if (length(jj) > 0 && is.finite(jj) && is.finite(dd[jj]) && dd[jj] <= peak_match_tol) {
            accepted_baseline_event_ids <- c(accepted_baseline_event_ids, eb_tmp$Event_ID[jj])
          }
        }
        accepted_baseline_event_ids <- unique(accepted_baseline_event_ids[is.finite(accepted_baseline_event_ids)])
      }
      
      filter_to_accepted_events <- function(x) {
        if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) return(x)
        if (length(accepted_baseline_event_ids) == 0) return(x[0, , drop = FALSE])
        if ("Event_ID" %in% names(x)) {
          return(x[x$Event_ID %in% accepted_baseline_event_ids, , drop = FALSE])
        }
        x
      }
      
      # Use the frozen non-sigmoid levels for plotting FWHP/FWHM.
      # This prevents the orange FWHP line or maroon FWHM line from moving
      # when the optional single/double sigmoid fits are activated.
      if (!"FWHP_Level_NoSigmoid" %in% names(df_p)) {
        df_p$FWHP_Level_NoSigmoid <- if ("Prominence_Midpoint" %in% names(df_p)) df_p$Prominence_Midpoint else numeric(nrow(df_p))
      }
      if (!"FWHM_Level_NoSigmoid" %in% names(df_p)) {
        df_p$FWHM_Level_NoSigmoid <- if ("Amplitude_Midpoint" %in% names(df_p)) df_p$Amplitude_Midpoint else numeric(nrow(df_p))
      }
      
      gg3 <- ggplot(data_smoothed, aes(x = Time, y = Sing)) +
        geom_line(linetype = "solid", linewidth = 1.5, color = "black") +
        geom_hline(yintercept = input$minpeakheight2, linetype = "dashed", color = "purple") +
        geom_point(data = df_p, aes(x = Peak_Occurence_Time, y = Amplitude), color = "red", size = 4) +
        geom_segment(data = df_p, aes(x = Peak_Occurence_Time, xend = Peak_Occurence_Time,
                                      y = Baseline_at_peak, yend = Amplitude),
                     linetype = "dashed", linewidth = 1.2, color = "red") +
        geom_segment(data = df_p, aes(x = Peak_Occurence_Time, xend = Peak_Occurence_Time,
                                      y = puntominimo_y, yend = Amplitude),
                     linetype = "dashed", linewidth = 1.2, color = "blue") +
        geom_segment(data = df_p, aes(x = Time_left_FWHP, xend = Time_right_FWHP,
                                      y = FWHP_Level_NoSigmoid, yend = FWHP_Level_NoSigmoid),
                     linetype = "solid", linewidth = 1.2, color = "orange") +
        labs(title = "Calcium Trace",
             x = "Time [s]",
             y = "Delta F/Baseline") +
        theme_classic() +
        theme(
          plot.title = element_text(size = 28, face = "bold"),
          axis.title.y = element_text(size = 28, face = "bold"),
          axis.title.x = element_text(size = 28, face = "bold"),
          axis.text.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold")
        )
      
      if (identical(as.character(input$raw_data), "2")) {
        gg3 <- gg3 + geom_line(data = data_raw, aes(x = Time, y = signal), color = "red", alpha = 0.6)
      }
      
      # Optional sigmoid curves and onset/offset annotations.
      sigmoid_curves <- Peaks_Data_Final()$sigmoid_curves
      if (!is.null(sigmoid_curves) && nrow(sigmoid_curves) > 0) {
        single_curves <- sigmoid_curves[sigmoid_curves$Model == "Single sigmoid", , drop = FALSE]
        double_curves <- sigmoid_curves[sigmoid_curves$Model == "Double sigmoid", , drop = FALSE]
        
        if (nrow(single_curves) > 0) {
          gg3 <- gg3 +
            geom_line(
              data = single_curves,
              aes(x = Time, y = Fitted, group = Event_ID),
              inherit.aes = FALSE,
              linewidth = 1.0,
              color = "#0072B2"
            )
        }
        
        if (nrow(double_curves) > 0) {
          gg3 <- gg3 +
            geom_line(
              data = double_curves,
              aes(x = Time, y = Fitted, group = Event_ID),
              inherit.aes = FALSE,
              linewidth = 1.0,
              color = "#7B3294"
            )
        }
      }
      
      if (nrow(df_p) > 0 && "SingleSigmoid_Onset_Time" %in% names(df_p)) {
        single_onsets <- df_p[is.finite(df_p$SingleSigmoid_Onset_Time), , drop = FALSE]
        if (nrow(single_onsets) > 0) {
          gg3 <- gg3 +
            geom_vline(
              data = single_onsets,
              aes(xintercept = SingleSigmoid_Onset_Time),
              inherit.aes = FALSE,
              linetype = "dotted",
              linewidth = 1.0,
              color = "#D55E00"
            )
        }
      }
      
      if (nrow(df_p) > 0 && "DoubleSigmoid_Offset_Time" %in% names(df_p)) {
        double_offsets <- df_p[is.finite(df_p$DoubleSigmoid_Offset_Time), , drop = FALSE]
        if (nrow(double_offsets) > 0) {
          gg3 <- gg3 +
            geom_vline(
              data = double_offsets,
              aes(xintercept = DoubleSigmoid_Offset_Time),
              inherit.aes = FALSE,
              linetype = "dotted",
              linewidth = 1.0,
              color = "#CC79A7"
            )
        }
      }
      
      if (identical(as.character(input$auc2), "2")) {
        Integration_Reference <- peaks_plot()$baseline1
        auc_df <- data_smoothed
        baseline_trace <- peaks_plot()$baseline_trace
        
        # Important: even when double sigmoid is active, the event-specific
        # AUC baseline remains the original calcium-transient left/right FFT-minima
        # baseline from peaks_plot()$baseline_trace. Double-sigmoid minima are not
        # used to draw the AUC ribbon or the event-specific baseline line.
        
        if (identical(as.character(input$Baseline), "8") &&
            !is.null(baseline_trace) &&
            nrow(baseline_trace) > 0) {
          
          # Event-specific baseline:
          # Draw the green AUC only on the restricted event-specific interval.
          # Do not interpolate the event-specific baseline over the full trace.
          event_auc_df <- baseline_trace
          event_auc_df <- filter_to_accepted_events(event_auc_df)
          event_auc_df$ymax_auc <- ifelse(
            event_auc_df$Signal > event_auc_df$Baseline,
            event_auc_df$Signal,
            NA_real_
          )
          
          gg3 <- gg3 +
            geom_ribbon(
              data = event_auc_df,
              aes(x = Time, ymin = Baseline, ymax = ymax_auc, group = Event_ID),
              fill = "green",
              alpha = 0.16,
              inherit.aes = FALSE
            )
          
        } else if (!is.null(baseline_trace) && nrow(baseline_trace) > 0) {
          auc_df$Baseline <- approx(
            x = baseline_trace$Time,
            y = baseline_trace$Baseline,
            xout = auc_df$Time,
            rule = 2
          )$y
          auc_df$ymax_auc <- ifelse(auc_df$Sing > auc_df$Baseline, auc_df$Sing, NA_real_)
          
          gg3 <- gg3 +
            geom_line(
              data = baseline_trace,
              aes(x = Time, y = Baseline),
              linetype = "dashed",
              linewidth = 1.2,
              color = "darkgreen",
              inherit.aes = FALSE
            ) +
            geom_ribbon(
              data = auc_df,
              aes(x = Time, ymin = Baseline, ymax = ymax_auc),
              fill = "green",
              alpha = 0.1,
              inherit.aes = FALSE
            )
        } else {
          auc_df$ymax_auc <- ifelse(auc_df$Sing > Integration_Reference, auc_df$Sing, NA_real_)
          gg3 <- gg3 +
            geom_hline(yintercept = Integration_Reference, linetype = "dashed", color = "green") +
            geom_ribbon(
              data = auc_df,
              aes(x = Time, ymin = Integration_Reference, ymax = ymax_auc),
              fill = "green",
              alpha = 0.1,
              inherit.aes = FALSE
            )
        }
      }
      
      if (as.character(input$Baseline) %in% c("6", "7", "8", "11") && !is.null(peaks_plot()$baseline_trace)) {
        baseline_trace_to_plot <- peaks_plot()$baseline_trace
        baseline_info_to_plot <- peaks_plot()$baseline_info
        
        # The event-specific baseline line must always use the original
        # calcium-transient left/right FFT minima. Do not replace these minima
        # with the fitted double-sigmoid minima when double sigmoid is active.
        if (identical(as.character(input$Baseline), "8")) {
          event_baseline_segments <- baseline_info_to_plot$event_baselines
          event_baseline_segments <- filter_to_accepted_events(event_baseline_segments)
          
          if (!is.null(event_baseline_segments) &&
              is.data.frame(event_baseline_segments) &&
              nrow(event_baseline_segments) > 0) {
            
            # Draw the event-specific baseline only on the clipped interval
            # that belongs to the accepted peak. Do not draw the full left-minimum
            # to right-minimum line, because that long line can pass through a
            # neighboring bump/transient that was not accepted as a peak.
            if (all(c("Plot_Left_Time", "Plot_Left_Value", "Plot_Right_Time", "Plot_Right_Value") %in% names(event_baseline_segments))) {
              event_baseline_segments$Segment_Left_Time <- event_baseline_segments$Plot_Left_Time
              event_baseline_segments$Segment_Left_Value <- event_baseline_segments$Plot_Left_Value
              event_baseline_segments$Segment_Right_Time <- event_baseline_segments$Plot_Right_Time
              event_baseline_segments$Segment_Right_Value <- event_baseline_segments$Plot_Right_Value
            } else if (all(c("Left_Min_Time", "Left_Min_Value", "Right_Min_Time", "Right_Min_Value") %in% names(event_baseline_segments))) {
              event_baseline_segments$Segment_Left_Time <- event_baseline_segments$Left_Min_Time
              event_baseline_segments$Segment_Left_Value <- event_baseline_segments$Left_Min_Value
              event_baseline_segments$Segment_Right_Time <- event_baseline_segments$Right_Min_Time
              event_baseline_segments$Segment_Right_Value <- event_baseline_segments$Right_Min_Value
            }
            
            if (all(c("Segment_Left_Time", "Segment_Left_Value", "Segment_Right_Time", "Segment_Right_Value") %in% names(event_baseline_segments))) {
              gg3 <- gg3 +
                geom_segment(
                  data = event_baseline_segments,
                  aes(
                    x = Segment_Left_Time,
                    xend = Segment_Right_Time,
                    y = Segment_Left_Value,
                    yend = Segment_Right_Value
                  ),
                  linetype = "dashed",
                  linewidth = 1.2,
                  color = "darkgreen",
                  inherit.aes = FALSE
                )
            }
          }
        } else {
          gg3 <- gg3 +
            geom_line(
              data = baseline_trace_to_plot,
              aes(x = Time, y = Baseline),
              linetype = "dashed",
              linewidth = 1.2,
              color = "darkgreen",
              inherit.aes = FALSE
            )
        }
      }
      
      # ------------------------------------------------------------
      # Show fixed consecutive baseline windows on calcium trace
      # First windows keep the same width; only the last one may be shorter.
      # Example: if Time = 0 to 100 and rolling_window = 20,
      # windows are 0-20, 20-40, 40-60, 60-80, 80-100.
      # If the total time is not divisible by rolling_window, the last window
      # ends at max(Time) and is allowed to be shorter.
      # ------------------------------------------------------------
      if (identical(as.character(input$Baseline), "7") && isTRUE(input$show_fixed_windows)) {
        window_size <- suppressWarnings(as.numeric(input$rolling_window))
        
        if (length(window_size) > 0 && is.finite(window_size) && window_size > 0) {
          time_range <- range(data_smoothed$Time, na.rm = TRUE)
          
          if (all(is.finite(time_range)) && diff(time_range) > 0) {
            window_starts <- seq(
              from = time_range[1],
              to = time_range[2],
              by = window_size
            )
            
            # Avoid creating a zero-width final rectangle if the last start
            # is exactly equal to max(Time).
            window_starts <- window_starts[window_starts < time_range[2]]
            
            if (length(window_starts) > 0) {
              window_df <- data.frame(
                xmin = window_starts,
                xmax = pmin(window_starts + window_size, time_range[2]),
                ymin = -Inf,
                ymax = Inf
              )
              
              gg3 <- gg3 +
                geom_rect(
                  data = window_df,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  inherit.aes = FALSE,
                  alpha = 0.04,
                  fill = "gray70",
                  color = "gray40",
                  linetype = "dashed",
                  linewidth = 0.7
                )
            }
          }
        }
      }
      
      if (identical(as.character(input$FWHM), "2") && nrow(df_p) > 0) {
        gg3 <- gg3 +
          geom_segment(
            data = df_p,
            aes(x = Time_left_FWHM, xend = Time_right_FWHM,
                y = FWHM_Level_NoSigmoid, yend = FWHM_Level_NoSigmoid),
            linetype = "solid",
            linewidth = 1.2,
            color = "maroon"
          )
      }
      
      list(gg3 = gg3)
    })
    
    output$plot_peak3 <- renderPlot({
      req(nrow(Peaks_Data_Final()$df_p) > 0)
      Trance_Graph()$gg3
    })
    
    output$plot_raw_smoothed <- renderPlot({
      peaks_plot()$gg2
    })
    
    output$plot_trend_removal <- renderPlot({
      
      shiny::validate(
        shiny::need(
          isTRUE(peaks_df()$trend_removed),
          "Trend removal is not selected."
        )
      )
      
      trend_df <- peaks_df()$trend_res$trend_data
      
      ggplot(trend_df, aes(x = Time)) +
        geom_line(
          aes(y = signal_original, color = "Original calcium trace"),
          linewidth = 1.2,
          alpha = 1
        ) +
        geom_line(
          aes(y = trend, color = "Fitted polynomial trend"),
          linewidth = 1.2,
          linetype = "dashed",
          alpha = 1
        ) +
        geom_line(
          aes(y = signal_detrended, color = "Detrended calcium trace"),
          linewidth = 1.2,
          alpha = 1
        ) +
        scale_color_manual(
          name = "Curves",
          values = c(
            "Original calcium trace" = "black",
            "Fitted polynomial trend" = "blue",
            "Detrended calcium trace" = "red"
          )
        ) +
        labs(
          title = "Polynomial Trend Removal",
          subtitle = "Original calcium trace, fitted polynomial trend, and detrended calcium trace",
          x = "Time [s]",
          y = "Delta F/Baseline"
        ) +
        theme_classic(base_size = 18) +
        theme(
          plot.title = element_text(
            size = 28,
            face = "bold",
            color = "black"
          ),
          plot.subtitle = element_text(
            size = 14,
            color = "black"
          ),
          axis.title.x = element_text(
            size = 28,
            face = "bold",
            color = "black"
          ),
          axis.title.y = element_text(
            size = 28,
            face = "bold",
            color = "black"
          ),
          axis.text.x = element_text(
            size = 18,
            face = "bold",
            color = "black"
          ),
          axis.text.y = element_text(
            size = 18,
            face = "bold",
            color = "black"
          ),
          axis.line = element_line(
            linewidth = 0.8,
            color = "black"
          ),
          axis.ticks = element_line(
            linewidth = 0.8,
            color = "black"
          ),
          legend.position = "bottom",
          legend.title = element_text(
            size = 14,
            face = "bold",
            color = "black"
          ),
          legend.text = element_text(
            size = 12,
            color = "black"
          ),
          panel.grid = element_blank()
        )
    })
    
    # ------------------------------------------------------------
    # Tables used by the interface and CSV downloads
    # ------------------------------------------------------------
    # These helpers are the single source of truth for the visible
    # Transient Metrics and Trace Metrics tables. The download handlers
    # below write exactly these same tables, so the CSV dimensions change
    # in the same way as the interface when single and/or double sigmoid
    # options are activated.
    # ------------------------------------------------------------
    round_numeric_columns <- function(df, digits = 3) {
      if (is.null(df) || nrow(df) == 0) return(df)
      numeric_cols <- vapply(df, is.numeric, logical(1))
      df[numeric_cols] <- lapply(df[numeric_cols], function(z) round(z, digits))
      df
    }
    
    build_transient_metrics_display_table <- function() {
      df_p <- Peaks_Data_Final()$df_p
      shiny::validate(shiny::need(
        !is.null(df_p) && nrow(df_p) > 0,
        "No peaks satisfy the selected thresholds."
      ))
      
      base_cols <- c(
        "Amplitude",
        "Baseline_at_peak",
        "Peak_Occurence_Time",
        "Prominence",
        "FWHP",
        "FWHM",
        "Transient_Ocurrence_Time",
        "Peak_Rise_Time",
        "Rise_Rate"
      )
      
      single_cols <- c(
        "SingleSigmoid_Slope_at_t0",
        "SingleSigmoid_R2"
      )
      
      double_cols <- c(
        "DoubleSigmoid_Decay_Time",
        "DoubleSigmoid_Decay_Slope"
      )
      
      selected_cols <- base_cols
      if (isTRUE(input$use_single_sigmoid_fft)) {
        selected_cols <- c(selected_cols, single_cols)
      }
      if (isTRUE(input$use_double_sigmoid_fft)) {
        selected_cols <- c(selected_cols, double_cols)
      }
      
      selected_cols <- selected_cols[selected_cols %in% names(df_p)]
      df_display <- df_p[, selected_cols, drop = FALSE]
      
      # The interface reports amplitude above the selected baseline.
      if (all(c("Amplitude", "Baseline_at_peak") %in% names(df_display))) {
        df_display$Amplitude <- df_display$Amplitude - df_display$Baseline_at_peak
      }
      
      round_numeric_columns(df_display, digits = 3)
    }
    
    build_trace_metrics_display_table <- function() {
      df_p <- Peaks_Data_Final()$df_p
      shiny::validate(shiny::need(
        !is.null(df_p) && nrow(df_p) > 0,
        "No trace metrics are available because no valid peaks were detected."
      ))
      
      time1 <- min(peaks_df()$data_raw$Time, na.rm = TRUE)
      time2 <- max(peaks_df()$data_raw$Time, na.rm = TRUE)
      duration <- time2 - time1
      
      Time_OnSet <- df_p$Transient_Ocurrence_Time[1]
      Frequency <- ifelse(is.finite(duration) && duration > 0, nrow(df_p) / duration, NA_real_)
      Baseline <- peaks_plot()$baseline1
      
      # Keep the trace-level Baseline equal to the selected event-specific
      # calcium-transient FFT-minima baseline. Double-sigmoid minima are not used
      # to redefine this baseline.
      
      trace_metrics <- data.frame(
        Time_Onset = Time_OnSet,
        Frequency = Frequency,
        Baseline = Baseline,
        Number_of_Peaks = nrow(df_p)
      )
      
      if (isTRUE(input$use_single_sigmoid_fft) &&
          "SingleSigmoid_Onset_Time" %in% names(df_p) &&
          any(is.finite(df_p$SingleSigmoid_Onset_Time))) {
        trace_metrics$Mean_SingleSigmoid_Onset_Time <- mean(df_p$SingleSigmoid_Onset_Time, na.rm = TRUE)
        trace_metrics$Mean_SingleSigmoid_Rise_Time <- mean(df_p$SingleSigmoid_Rise_Time, na.rm = TRUE)
        trace_metrics$Mean_SingleSigmoid_Rise_Rate <- mean(df_p$SingleSigmoid_Rise_Rate, na.rm = TRUE)
        trace_metrics$Mean_SingleSigmoid_Slope_at_t0 <- mean(df_p$SingleSigmoid_Slope_at_t0, na.rm = TRUE)
      }
      
      if (isTRUE(input$use_double_sigmoid_fft) &&
          "DoubleSigmoid_Decay_Time" %in% names(df_p) &&
          any(is.finite(df_p$DoubleSigmoid_Decay_Time))) {
        trace_metrics$Mean_DoubleSigmoid_Decay_Time <- mean(
          df_p$DoubleSigmoid_Decay_Time,
          na.rm = TRUE
        )
      }
      
      if (isTRUE(input$use_double_sigmoid_fft) &&
          "DoubleSigmoid_Decay_Slope" %in% names(df_p) &&
          any(is.finite(df_p$DoubleSigmoid_Decay_Slope))) {
        trace_metrics$Mean_DoubleSigmoid_Decay_Slope <- mean(
          df_p$DoubleSigmoid_Decay_Slope,
          na.rm = TRUE
        )
      }
      
      if (identical(as.character(input$auc2), "2")) {
        tabla_auc_use <- peaks_plot()$tabla_AUC
        
        # Keep AUC based on the selected event-specific calcium-transient
        # FFT-minima baseline. Double-sigmoid minima are not used for AUC.
        
        trace_metrics <- cbind(trace_metrics, tabla_auc_use)
      }
      
      round_numeric_columns(trace_metrics, digits = 3)
    }
    
    output$table_peaks2 <- renderDT({
      df_display <- build_transient_metrics_display_table()
      datatable(
        df_display,
        rownames = FALSE,
        options = list(scrollX = TRUE),
        caption = tags$caption(tags$strong("Transient Metrics"))
      )
    })
    
    output$table_peaks22 <- renderDT({
      trace_metrics <- build_trace_metrics_display_table()
      datatable(
        trace_metrics,
        rownames = FALSE,
        options = list(pagingType = "simple", dom = "t", autoWidth = TRUE, scrollX = TRUE),
        caption = tags$caption(tags$strong("Trace Metrics"))
      )
    })
    
    output$descargarP <- downloadHandler(
      filename = function() "Trace_Metrics.csv",
      content = function(file) {
        trace_metrics <- build_trace_metrics_display_table()
        write.csv(trace_metrics, file, row.names = FALSE)
      }
    )
    
    output$descargar <- downloadHandler(
      filename = function() "Transient_Metrics.csv",
      content = function(file) {
        transient_metrics <- build_transient_metrics_display_table()
        write.csv(transient_metrics, file, row.names = FALSE)
      }
    )
    
    
    output$Calcium_Trance_Graph <- downloadHandler(
      filename = function() paste0("calcium_trace_", Sys.Date(), ".png"),
      content = function(file) {
        ggsave(file, plot = Trance_Graph()$gg3, dpi = 300, width = 10, height = 6)
      }
    )
    
    output$fft_infoText <- renderPrint({
      x <- fft_panel_data()
      omega_cut <- 2 * pi * (x$k_keep - 1) / x$n
      
      cat("ROI selected:", input$Cell2, "\n")
      cat("n =", x$n, "\n")
      cat("Positive unique frequencies (including DC):", x$n_pos, "\n")
      cat("k_keep =", x$k_keep, "\n")
      cat("Approximate cutoff ω_c =", round(omega_cut, 4), "rad\n")
      cat("Low-frequency fraction =", input$fft_fraction, "\n")
      cat("Keep mean (DC) =", isTRUE(input$keep_mean), "\n")
    })
    
    output$fft_fitPlot <- renderPlot({
      df2 <- fft_panel_data()$df_out
      
      ggplot(df2, aes(x = t, y = x_t)) +
        geom_line(linewidth = 0.35, alpha = 0.7) +
        geom_point(alpha = 0.6) +
        geom_line(aes(y = s_hat), color = "blue", linewidth = 1.1) +
        labs(
          title = "ROI signal and FFT-smoothed signal",
          subtitle = "Low-pass reconstruction using retained low frequencies",
          x = "t (index)",
          y = "Signal"
        ) +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor = element_blank())
    })
    
    output$fft_residPlot <- renderPlot({
      df2 <- fft_panel_data()$df_out
      
      ggplot(df2, aes(x = t, y = resid)) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_line(linewidth = 0.9) +
        geom_point(alpha = 0.8, size = 1.8) +
        labs(
          title = "Residuals",
          subtitle = "Observed signal minus FFT-smoothed signal",
          x = "t (index)",
          y = "Residual"
        ) +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor = element_blank())
    })
    
    output$fft_specPlot <- renderPlot({
      x <- fft_panel_data()
      dfp <- x$df_freq
      omega_cut <- 2 * pi * (x$k_keep - 1) / x$n
      
      ggplot(dfp, aes(x = omega, y = I)) +
        geom_line(linewidth = 0.9) +
        geom_vline(xintercept = omega_cut, linetype = 2) +
        labs(
          title = "Periodogram and cutoff frequency",
          subtitle = paste0("Frequencies retained up to k = ", x$k_keep - 1),
          x = expression(omega),
          y = "I(omega)"
        ) +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor = element_blank())
    })
    
    output$fft_freqTable <- renderDT({
      df_freq <- fft_panel_data()$df_freq
      num_cols <- sapply(df_freq, is.numeric)
      df_freq[num_cols] <- lapply(df_freq[num_cols], round, 6)
      
      datatable(
        df_freq,
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 8)
      )
    })
    
    output$fft_formula_text <- renderText({
      fft_panel_data()$formula_txt
    })
    
    output$fft_formula_box <- renderUI({
      tagList(
        tags$div(
          class = "helper-card",
          tags$h4("1. What is shown here?"),
          tags$p("This panel explains how the selected ROI signal is decomposed into frequencies using the Fast Fourier Transform (FFT)."),
          tags$p("The smoothing is obtained by keeping only the low-frequency components and removing the high-frequency ones.")
        ),
        tags$div(
          class = "helper-card",
          tags$h4("2. Frequency decomposition"),
          tags$p(HTML("In R, the transform is computed with <code>X &lt;- fft(x)</code>.")),
          tags$p(HTML("Each coefficient <b>X_k</b> represents how much of frequency <b>k</b> is present in the signal.")),
          tags$p(HTML("<b>k</b>: frequency index.")),
          tags$p(HTML("<b>omega</b>: angular frequency in radians.")),
          tags$p(HTML("<b>freq</b>: cycles per observation.")),
          tags$p(HTML("<b>period</b>: approximate number of observations per cycle."))
        ),
        tags$div(
          class = "helper-card",
          tags$h4("3. Why smoothing works"),
          tags$p("Rapid oscillations and noise are usually represented by high frequencies."),
          tags$p("By keeping only the lowest frequencies, the reconstructed signal preserves the main shape of the calcium trace while reducing noise."),
          tags$p("The parameter 'FFT Low-Frequency Fraction' controls how many low frequencies are retained.")
        ),
        tags$div(
          class = "helper-card",
          tags$h4("4. Explicit formula"),
          tags$p("The reconstructed signal can be written explicitly as a finite sum of a constant term plus cosine and sine terms associated with the kept frequencies.")
        )
      )
    })
    
    
    output$trend_poly_info <- renderText({
      shiny::validate(
        shiny::need(
          isTRUE(peaks_df()$trend_removed),
          "Trend removal is not selected."
        )
      )
      
      trend_info <- build_polynomial_trend_info(peaks_df()$trend_res, digits = 6)
      trend_info$info_text
    })
    
    output$trend_coef_table <- renderDT({
      shiny::validate(
        shiny::need(
          isTRUE(peaks_df()$trend_removed),
          "Trend removal is not selected."
        )
      )
      
      trend_info <- build_polynomial_trend_info(peaks_df()$trend_res, digits = 6)
      
      datatable(
        trend_info$coef_table,
        options = list(
          dom = "t",
          paging = FALSE,
          searching = FALSE,
          info = FALSE,
          lengthChange = FALSE,
          scrollX = TRUE
        ),
        rownames = FALSE,
        caption = tags$caption(
          tags$strong("Polynomial trend coefficients")
        )
      )
    })
    
    output$fft_grid_current_params <- renderUI({
      baseline_label <- switch(
        as.character(fft_grid_fixed_params$baseline_mode),
        "2" = "Standard definition",
        "3" = "Low-fluorescence region baseline",
        "6" = "Rolling percentile baseline",
        "8" = "Event-specific baseline: left-right FFT minima",
        "11" = "Constant baseline: mean from trace start to Time_Onset",
        "9" = "Median signal (legacy)",
        "10" = "Mean signal (legacy)",
        "5" = "Min",
        "Unknown"
      )
      
      auc_label <- if (as.character(fft_grid_fixed_params$auc) == "2") "Yes" else "No"
      fwhm_label <- if (as.character(fft_grid_fixed_params$fwhm) == "2") "Yes" else "No"
      
      tagList(
        tags$p(HTML(paste0("<b>1. Peak Height (min):</b> ", fft_grid_fixed_params$minpeakheight))),
        tags$p(HTML(paste0("<b>2. Peak Ascent:</b> ", fft_grid_fixed_params$nups))),
        tags$p(HTML(paste0("<b>3. Peak Descent:</b> ", fft_grid_fixed_params$ndowns))),
        tags$p(HTML(paste0("<b>4. Min Peak Distance:</b> ", fft_grid_fixed_params$minpeakdistance))),
        tags$p(HTML(paste0("<b>5. FWHP (min):</b> ", fft_grid_fixed_params$min_FWHP))),
        tags$p(HTML(paste0("<b>6. Prominence (min):</b> ", fft_grid_fixed_params$min_prominence))),
        tags$p(HTML(paste0("<b>Baseline:</b> ", baseline_label))),
        tags$p(HTML(paste0("<b>Area Under the Curve (AUC):</b> ", auc_label))),
        tags$p(HTML(paste0("<b>Full Width at Half Maximum (FWHM):</b> ", fwhm_label)))
      )
    })
    
    output$fft_grid_recommendation <- renderUI({
      x <- fft_fraction_grid()
      tags$div(
        class = "helper-card",
        h4("Recommended FFT Low-Frequency Fraction"),
        tags$p(tags$strong(x$recommendation_text)),
        #tags$p("Thresholds used by the heuristic: correlation(raw, smoothed) >= 0.90, residual SD <= 1.5 x the residual SD at f = 0.40, retained spectral power >= 90%, and criterion 4 based on closeness of mean_peak to the weighted mean of the candidate mean_peak values.")
      )
    })
    
    output$fft_grid_summary <- renderDT({
      grid <- fft_fraction_grid()$summary_df
      
      grid_display <- grid[, c(
        "f",
        "k_keep",
        "corr_raw_smooth",
        "retained_power_frac",
        "n_peaks",
        "mean_peak"
      ), drop = FALSE]
      
      datatable(
        grid_display,
        options = list(
          dom = "t",
          paging = FALSE,
          searching = FALSE,
          info = FALSE,
          lengthChange = FALSE,
          scrollX = TRUE
        ),
        rownames = FALSE,
        caption = tags$caption(
          tags$strong("Grid assessment for candidate FFT fractions")
        )
      )
    })
    
    output$fft_grid_metrics_note <- renderUI({
      x <- fft_fraction_grid()
      tags$div(
        class = "helper-card",
        h4("Criterion 4 status"),
        tags$p(tags$strong(x$criterion4_note))
      )
    })
    
    output$fft_grid_signal_plot <- renderPlot({
      df <- fft_fraction_grid()$overlay_df
      ggplot(df, aes(x = Time, y = value, color = series)) +
        geom_line(linewidth = 0.7) +
        facet_wrap(~ f_label, ncol = 2, scales = "free_y") +
        labs(
          title = "Criterion 1: Raw vs. FFT-smoothed signal",
          subtitle = "Each facet shows the same ROI smoothed with a different FFT Low-Frequency Fraction",
          x = "Time [s]",
          y = "Delta F/Baseline"
        ) +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor = element_blank())
    })
    
    output$fft_grid_resid_plot <- renderPlot({
      df <- fft_fraction_grid()$residual_df
      ggplot(df, aes(x = Time, y = resid)) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_line(linewidth = 0.7) +
        facet_wrap(~ f_label, ncol = 2, scales = "free_y") +
        labs(
          title = "Criterion 2: Residuals",
          subtitle = "Residual = raw signal - smoothed signal",
          x = "Time [s]",
          y = "Residual"
        ) +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor = element_blank())
    })
    
    output$fft_grid_spec_plot <- renderPlot({
      x <- fft_fraction_grid()
      
      ann_df <- merge(
        x$cutoff_df,
        x$summary_df[, c("f", "k_keep")],
        by = "f",
        all.x = TRUE
      )
      
      ann_df$label <- paste0(
        "k_keep = ", ann_df$k_keep,
        "\n\u03C9_c = ", sprintf("%.3f", ann_df$omega_cut)
      )
      
      ggplot(x$periodogram_df, aes(x = omega, y = I)) +
        geom_line(linewidth = 0.7) +
        geom_vline(
          data = ann_df,
          aes(xintercept = omega_cut),
          linetype = 2,
          inherit.aes = FALSE
        ) +
        geom_text(
          data = ann_df,
          aes(x = Inf, y = Inf, label = label),
          inherit.aes = FALSE,
          hjust = 1.05,
          vjust = 1.15,
          size = 4
        ) +
        facet_wrap(~ f_label, ncol = 2, scales = "free_y") +
        labs(
          title = "Criterion 3: Periodogram and cutoff",
          subtitle = "Dashed line = approximate cutoff frequency for each candidate",
          x = expression(omega),
          y = "I(omega)"
        ) +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor = element_blank())
    })
    
    output$fft_grid_metrics_plot <- renderPlot({
      x <- fft_fraction_grid()
      df <- x$mean_peak_df
      
      if (!is.finite(x$weighted_mean_peak) || nrow(df) == 0 || all(!is.finite(df$mean_peak))) {
        plot.new()
        title(main = "Criterion 4: Mean peak closeness to weighted mean")
        text(0.5, 0.5, "No valid mean_peak values are available under the current thresholds.")
      } else {
        ggplot(df, aes(x = f, y = mean_peak)) +
          geom_line(linewidth = 0.8) +
          geom_point(aes(shape = recommended), size = 3) +
          geom_hline(yintercept = x$weighted_mean_peak, linetype = 3, color = "blue") +
          geom_vline(xintercept = x$recommended_f, linetype = 2, color = "red") +
          scale_x_continuous(breaks = sort(unique(df$f))) +
          scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), guide = "none") +
          labs(
            title = "Criterion 4: Mean peak closeness to weighted mean",
            subtitle = paste0(
              "Blue dotted line = weighted mean_peak (", sprintf("%.4f", x$weighted_mean_peak),
              "); red dashed line = recommended fraction (", sprintf("%.2f", x$recommended_f), ")"
            ),
            x = "FFT Low-Frequency Fraction",
            y = "mean_peak"
          ) +
          theme_minimal(base_size = 13) +
          theme(panel.grid.minor = element_blank())
      }
    })  })
}


