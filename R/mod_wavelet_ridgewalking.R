# Dependencies are declared in DESCRIPTION/NAMESPACE; helpers are loaded through DESCRIPTION Collate.
# R/mod_wavelet_ridgewalking.R
# ============================================================
# Wavelet ridgewalking-inspired event detection module
# CalciumInsights
#
# This version adds:
# 1) Event_Status classification: Accepted / Low confidence / Rejected
# 2) Biological post-filtering after wavelet ridge detection
# 3) Event_AUC shaded in orange in the figure
# 4) Total_AUC_Global_Above_Baseline shaded in blue in the figure
# 5) AUC summary table with total area between the curve and baseline
# 6) Wavelet-derived metrics comparable with FFT + Baseline Analysis.
# 7) Amplitude and prominence are drawn on the detected-events figure.
# 8) Prominence is calculated using FFT smoothing with f = event_specific_fft_f only for prominence.
# 9) Prominence remains independent of single-sigmoid and double-sigmoid fits.
# 10) Trace-start-to-onset baseline is synchronized after sigmoid onset refinement.
# 11) Duplicate wavelet events with the same Peak_Time are removed before metrics are reported.
# ============================================================

# ------------------------------------------------------------
# Example data
# ------------------------------------------------------------
example_calcium_wavelet_data <- function() {
  set.seed(123)
  t <- seq(0, 100, by = 0.5)
  
  mk_signal <- function(shift = 0, noise = 0.12) {
    baseline <- 0.12 * sin(2 * pi * t / 40) + 0.04 * t / max(t)
    peaks <- exp(-0.5 * ((t - (18 + shift)) / 1.8)^2) +
      1.1 * exp(-0.5 * ((t - (46 + shift)) / 2.5)^2) +
      0.85 * exp(-0.5 * ((t - (78 + shift)) / 3.0)^2)
    0.45 + baseline + peaks + rnorm(length(t), 0, noise)
  }
  
  data.frame(
    Time = t,
    ROI_1 = mk_signal(0),
    ROI_2 = mk_signal(3),
    ROI_3 = mk_signal(-2),
    ROI_4 = mk_signal(5)
  )
}

# ------------------------------------------------------------
# Basic utilities
# ------------------------------------------------------------
standardize_signal <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  sx <- stats::sd(x, na.rm = TRUE)
  mx <- mean(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) return(x - mx)
  as.numeric((x - mx) / sx)
}

safe_num <- function(x, fallback = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || !is.finite(x[1])) return(fallback)
  x[1]
}

safe_int <- function(x, fallback = 1L) {
  x <- suppressWarnings(as.integer(x))
  if (length(x) == 0 || !is.finite(x[1])) return(fallback)
  x[1]
}


# ============================================================
# Optional sigmoid modeling for the Wavelet Ridgewalking module
# ============================================================
# The wavelet module keeps wavelet ridgewalking as the event detector.
# When the user activates the optional sigmoid models, the selected ROI trace
# is first smoothed with a Fourier low-pass reconstruction and the sigmoid
# models are fitted on that FFT-smoothed trace. This mirrors the logic used in
# the FFT module while preserving the wavelet-based event detection step.
# ============================================================

wavelet_fft_lowpass_smooth_trace <- function(time,
                                             signal,
                                             f = 0.20,
                                             keep_mean = TRUE) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  n <- length(signal)
  
  empty <- list(
    data = data.frame(Time = time, signal = signal),
    fft = complex(length = 0),
    fft_filt = complex(length = 0),
    k_keep = NA_integer_,
    n_pos = NA_integer_,
    keep_idx = integer(0)
  )
  
  if (length(time) != n || n < 4 || all(!is.finite(signal))) {
    return(empty)
  }
  
  if (any(!is.finite(signal))) {
    idx <- seq_along(signal)
    good <- is.finite(signal)
    if (sum(good) >= 2) {
      signal[!good] <- stats::approx(
        x = idx[good],
        y = signal[good],
        xout = idx[!good],
        rule = 2
      )$y
    } else {
      return(empty)
    }
  }
  
  f <- safe_num(f, 0.20)
  f <- max(0.01, min(1, f))
  
  X <- stats::fft(signal)
  n_pos <- floor(n / 2) + 1L
  k_keep <- max(1L, floor(f * n_pos))
  neg_count <- max(0L, k_keep - 1L)
  neg_idx <- if (neg_count > 0L) {
    seq.int(from = n - neg_count + 1L, to = n)
  } else {
    integer(0)
  }
  keep_idx <- unique(c(seq.int(1L, k_keep), neg_idx))
  
  Xf <- complex(length = n)
  Xf[keep_idx] <- X[keep_idx]
  if (!isTRUE(keep_mean)) Xf[1] <- 0
  
  smoothed <- Re(stats::fft(Xf, inverse = TRUE) / n)
  
  list(
    data = data.frame(Time = time, signal = smoothed),
    fft = X,
    fft_filt = Xf,
    k_keep = k_keep,
    n_pos = n_pos,
    keep_idx = keep_idx
  )
}

wavelet_sigmoid_simple_fun <- function(t, Baseline, A, k, t0) {
  Baseline + A * stats::plogis(k * (t - t0))
}

wavelet_double_sigmoid_fun <- function(t, Baseline, A, k1, t1, k2, t2) {
  s1 <- stats::plogis(k1 * (t - t1))
  s2 <- stats::plogis(k2 * (t - t2))
  Baseline + A * (s1 - s2)
}

wavelet_safe_r2 <- function(y, yhat) {
  y <- suppressWarnings(as.numeric(y))
  yhat <- suppressWarnings(as.numeric(yhat))
  ok <- is.finite(y) & is.finite(yhat)
  if (sum(ok) < 3) return(NA_real_)
  ss_res <- sum((y[ok] - yhat[ok])^2, na.rm = TRUE)
  ss_tot <- sum((y[ok] - mean(y[ok], na.rm = TRUE))^2, na.rm = TRUE)
  if (!is.finite(ss_tot) || ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}

wavelet_compute_event_baseline_rms <- function(time,
                                               signal,
                                               reference_time,
                                               pre_event_window = 5,
                                               fallback_baseline = NA_real_) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  reference_time <- safe_num(reference_time, NA_real_)
  pre_event_window <- safe_num(pre_event_window, 5)
  fallback_baseline <- safe_num(fallback_baseline, NA_real_)
  
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

wavelet_find_adjacent_left_minimum <- function(time,
                                               signal,
                                               peak_time,
                                               prev_peak_time = NA_real_) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  peak_time <- safe_num(peak_time, NA_real_)
  prev_peak_time <- safe_num(prev_peak_time, NA_real_)
  
  empty <- list(time = NA_real_, value = NA_real_, index = NA_integer_)
  if (length(time) != length(signal) || length(time) < 3 || !is.finite(peak_time)) return(empty)
  
  left_bound_time <- min(time, na.rm = TRUE)
  if (is.finite(prev_peak_time) && prev_peak_time < peak_time) {
    left_bound_time <- max(left_bound_time, (prev_peak_time + peak_time) / 2)
  }
  
  idx <- which(time >= left_bound_time & time < peak_time)
  idx <- idx[is.finite(signal[idx]) & is.finite(time[idx])]
  if (length(idx) == 0) return(empty)
  
  local_idx <- idx[idx > 1 & idx < length(signal)]
  local_idx <- local_idx[
    signal[local_idx] <= signal[local_idx - 1L] &
      signal[local_idx] <= signal[local_idx + 1L]
  ]
  
  candidate <- if (length(local_idx) > 0) {
    max(local_idx, na.rm = TRUE)
  } else {
    idx[which.min(signal[idx])]
  }
  
  list(time = time[candidate], value = signal[candidate], index = candidate)
}

fit_wavelet_single_sigmoid_event <- function(time,
                                             signal,
                                             left_reference_time,
                                             peak_time,
                                             peak_value,
                                             fallback_baseline = NA_real_,
                                             pre_event_window = 5,
                                             fit_left_pad = 0,
                                             fit_right_pad = 0,
                                             dense_n = 300,
                                             left_reference_value = NA_real_) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  left_reference_time <- safe_num(left_reference_time, NA_real_)
  peak_time <- safe_num(peak_time, NA_real_)
  peak_value <- safe_num(peak_value, NA_real_)
  fallback_baseline <- safe_num(fallback_baseline, NA_real_)
  left_reference_value <- safe_num(left_reference_value, NA_real_)
  
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
  
  baseline_info <- wavelet_compute_event_baseline_rms(
    time = time,
    signal = signal,
    reference_time = left_reference_time,
    pre_event_window = pre_event_window,
    fallback_baseline = fallback_baseline
  )
  
  Baseline_fixed <- baseline_info$Baseline
  if (is.finite(left_reference_value)) Baseline_fixed <- min(Baseline_fixed, left_reference_value)
  threshold <- Baseline_fixed + baseline_info$RMS
  fit_left_pad <- safe_num(fit_left_pad, 0)
  fit_right_pad <- safe_num(fit_right_pad, 0)
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
    peak_value <- stats::approx(time, signal, xout = peak_time, rule = 2)$y
  }
  
  A_start <- max(peak_value - Baseline_fixed, diff(range(fit_data$y, na.rm = TRUE)), 1e-6)
  rise_span <- peak_time - left_reference_time
  if (!is.finite(rise_span) || rise_span <= 0) rise_span <- stats::median(diff(time), na.rm = TRUE)
  if (!is.finite(rise_span) || rise_span <= 0) rise_span <- 1
  
  k_start <- 4 / rise_span
  t0_start <- left_reference_time + 0.5 * (peak_time - left_reference_time)
  
  fit_obj <- tryCatch(
    stats::nls(
      y ~ Baseline_fixed + A * stats::plogis(k * (t - t0)),
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
  yhat_train <- wavelet_sigmoid_simple_fun(fit_data$t, Baseline_fixed, pars[["A"]], pars[["k"]], pars[["t0"]])
  r2 <- wavelet_safe_r2(fit_data$y, yhat_train)
  
  dense_t <- seq(fit_start, fit_end, length.out = dense_n)
  dense_y <- wavelet_sigmoid_simple_fun(dense_t, Baseline_fixed, pars[["A"]], pars[["k"]], pars[["t0"]])
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

wavelet_find_double_sigmoid_decay_end <- function(time,
                                                  signal,
                                                  peak_time,
                                                  left_reference_time = NA_real_,
                                                  fallback_baseline = NA_real_,
                                                  pre_event_window = 5,
                                                  next_peak_time = NA_real_,
                                                  min_right_time = NA_real_) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  peak_time <- safe_num(peak_time, NA_real_)
  left_reference_time <- safe_num(left_reference_time, NA_real_)
  fallback_baseline <- safe_num(fallback_baseline, NA_real_)
  pre_event_window <- safe_num(pre_event_window, 5)
  next_peak_time <- safe_num(next_peak_time, NA_real_)
  min_right_time <- safe_num(min_right_time, NA_real_)
  
  if (length(time) != length(signal) || length(time) < 4 || !is.finite(peak_time)) {
    return(NA_real_)
  }
  if (all(!is.finite(time)) || all(!is.finite(signal))) return(NA_real_)
  
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
  if (length(pidx) == 0 || !is.finite(pidx) || pidx >= length(time)) return(NA_real_)
  
  time_min <- min(time, na.rm = TRUE)
  time_max <- max(time, na.rm = TRUE)
  if (!is.finite(time_min) || !is.finite(time_max) || time_max <= time_min) return(NA_real_)
  
  right_bound_time <- if (is.finite(next_peak_time) && next_peak_time > peak_time) {
    min(time_max, (peak_time + next_peak_time) / 2)
  } else {
    time_max
  }
  if (!is.finite(right_bound_time) || right_bound_time <= peak_time) right_bound_time <- time_max
  
  right_bound_idx <- max(which(time <= right_bound_time))
  if (length(right_bound_idx) == 0 || !is.finite(right_bound_idx) || right_bound_idx <= pidx) {
    right_bound_idx <- length(time)
  }
  
  if (!is.finite(min_right_time) || min_right_time <= peak_time) min_right_time <- peak_time
  min_right_time <- min(min_right_time, time[right_bound_idx])
  
  baseline_reference_time <- if (is.finite(left_reference_time)) left_reference_time else peak_time
  local_baseline <- wavelet_compute_event_baseline_rms(
    time = time,
    signal = signal,
    reference_time = baseline_reference_time,
    pre_event_window = pre_event_window,
    fallback_baseline = fallback_baseline
  )
  
  baseline_level <- local_baseline$Baseline
  rms_level <- local_baseline$RMS
  peak_value <- signal[pidx]
  
  if (!is.finite(baseline_level)) baseline_level <- stats::median(signal, na.rm = TRUE)
  if (!is.finite(rms_level) || rms_level < 0) rms_level <- 0
  
  amplitude_above_baseline <- peak_value - baseline_level
  if (!is.finite(amplitude_above_baseline)) amplitude_above_baseline <- diff(range(signal, na.rm = TRUE))
  if (!is.finite(amplitude_above_baseline) || amplitude_above_baseline < 0) amplitude_above_baseline <- 0
  
  end_threshold <- baseline_level + max(rms_level, 0.05 * amplitude_above_baseline, 1e-8)
  
  search_idx <- seq.int(pidx + 1L, right_bound_idx)
  search_idx <- search_idx[time[search_idx] >= min_right_time]
  if (length(search_idx) == 0) return(time[right_bound_idx])
  
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
    after_return <- search_idx[search_idx >= min(below_threshold)]
    local_after_return <- local_minima_in(after_return)
    if (length(local_after_return) > 0) {
      candidate_idx <- local_after_return[1]
    } else {
      candidate_idx <- below_threshold[1]
    }
  } else {
    local_candidates <- local_minima_in(search_idx)
    if (length(local_candidates) > 0) {
      candidate_idx <- local_candidates[which.min(signal[local_candidates])]
    } else {
      candidate_idx <- search_idx[which.min(signal[search_idx])]
    }
  }
  
  if (length(candidate_idx) == 0 || !is.finite(candidate_idx)) return(time[right_bound_idx])
  
  candidate_time <- time[candidate_idx]
  if (!is.finite(candidate_time) || candidate_time <= peak_time) candidate_time <- time[right_bound_idx]
  
  candidate_time
}

fit_wavelet_double_sigmoid_event <- function(time,
                                             signal,
                                             left_reference_time,
                                             peak_time,
                                             right_reference_time,
                                             peak_value,
                                             fallback_baseline = NA_real_,
                                             pre_event_window = 5,
                                             fit_left_pad = 0,
                                             fit_right_pad = 0,
                                             dense_n = 400,
                                             left_reference_value = NA_real_,
                                             right_reference_value = NA_real_) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  left_reference_time <- safe_num(left_reference_time, NA_real_)
  peak_time <- safe_num(peak_time, NA_real_)
  right_reference_time <- safe_num(right_reference_time, NA_real_)
  peak_value <- safe_num(peak_value, NA_real_)
  fallback_baseline <- safe_num(fallback_baseline, NA_real_)
  left_reference_value <- safe_num(left_reference_value, NA_real_)
  right_reference_value <- safe_num(right_reference_value, NA_real_)
  
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
  
  baseline_info <- wavelet_compute_event_baseline_rms(
    time = time,
    signal = signal,
    reference_time = left_reference_time,
    pre_event_window = pre_event_window,
    fallback_baseline = fallback_baseline
  )
  
  Baseline_fixed <- baseline_info$Baseline
  finite_refs <- c(left_reference_value, right_reference_value)
  finite_refs <- finite_refs[is.finite(finite_refs)]
  if (length(finite_refs) > 0) {
    Baseline_fixed <- min(c(Baseline_fixed, finite_refs), na.rm = TRUE)
  }
  threshold <- Baseline_fixed + baseline_info$RMS
  fit_left_pad <- safe_num(fit_left_pad, 0)
  fit_right_pad <- safe_num(fit_right_pad, 0)
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
  
  if (!is.finite(peak_value)) peak_value <- stats::approx(time, signal, xout = peak_time, rule = 2)$y
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
        stats::plogis(k1 * (t - t1)) - stats::plogis(k2 * (t - t2))
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
  yhat_train <- wavelet_double_sigmoid_fun(
    fit_data$t,
    Baseline_fixed,
    pars[["A"]],
    pars[["k1"]],
    pars[["t1"]],
    pars[["k2"]],
    pars[["t2"]]
  )
  r2 <- wavelet_safe_r2(fit_data$y, yhat_train)
  
  dense_t <- seq(fit_start, fit_end, length.out = dense_n)
  dense_y <- wavelet_double_sigmoid_fun(
    dense_t,
    Baseline_fixed,
    pars[["A"]],
    pars[["k1"]],
    pars[["t1"]],
    pars[["k2"]],
    pars[["t2"]]
  )
  
  model_peak_idx <- which.max(dense_y)
  model_peak_time <- dense_t[model_peak_idx]
  model_peak_value <- dense_y[model_peak_idx]
  
  left_min_time <- left_reference_time
  left_min_value <- if (is.finite(left_reference_value)) left_reference_value else stats::approx(time, signal, xout = left_reference_time, rule = 2)$y
  right_min_time <- right_reference_time
  right_min_value <- if (is.finite(right_reference_value)) right_reference_value else stats::approx(time, signal, xout = right_reference_time, rule = 2)$y
  ds_baseline_at_peak <- NA_real_
  ds_prominence <- NA_real_
  
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
  
  ds_onset <- if (length(onset_idx) > 0) dense_t[min(onset_idx)] else left_reference_time
  ds_offset <- if (length(offset_idx) > 0) dense_t[min(offset_idx)] else right_reference_time
  
  ds_duration <- ifelse(is.finite(ds_onset) && is.finite(ds_offset), ds_offset - ds_onset, NA_real_)
  ds_rise_time <- ifelse(is.finite(ds_onset) && is.finite(peak_time), peak_time - ds_onset, NA_real_)
  ds_decay_time <- ifelse(is.finite(ds_offset) && is.finite(peak_time), ds_offset - peak_time, NA_real_)
  
  if (!is.finite(ds_decay_time) && is.finite(right_min_time) && is.finite(peak_time) && right_min_time > peak_time) {
    ds_decay_time <- right_min_time - peak_time
  }
  if (!is.finite(ds_duration) && is.finite(ds_onset) && is.finite(right_min_time) && right_min_time > ds_onset) {
    ds_duration <- right_min_time - ds_onset
  }
  
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
  
  if (!is.finite(decay_slope)) {
    dy <- diff(dense_y) / diff(dense_t)
    t_mid <- (head(dense_t, -1) + tail(dense_t, -1)) / 2
    right_deriv_idx <- which(t_mid >= model_peak_time & is.finite(dy))
    if (length(right_deriv_idx) > 0) decay_slope <- min(dy[right_deriv_idx], na.rm = TRUE)
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

append_wavelet_sigmoid_metrics <- function(metrics,
                                           time,
                                           signal,
                                           use_single_sigmoid = FALSE,
                                           use_double_sigmoid = FALSE,
                                           pre_event_window = 5,
                                           single_left_pad = 0,
                                           single_right_pad = 0,
                                           double_left_pad = 0,
                                           double_right_pad = 0,
                                           fft_fraction = NA_real_) {
  if (is.null(metrics) || nrow(metrics) == 0) {
    return(list(metrics = metrics, fitted_curves = data.frame()))
  }
  
  if (!isTRUE(use_single_sigmoid) && !isTRUE(use_double_sigmoid)) {
    return(list(metrics = metrics, fitted_curves = data.frame()))
  }
  
  curve_list <- list()
  metrics_out <- metrics
  
  # ------------------------------------------------------------
  # Freeze FWHP and FWHM before optional sigmoid modeling.
  # ------------------------------------------------------------
  # The single sigmoid can refine onset/rise-time metrics, and the double
  # sigmoid can add decay/full-transient metrics and update prominence or
  # event-specific baseline when requested. However, the width metrics FWHP
  # and FWHM must remain the original values computed before any sigmoid fit.
  # These columns are restored at the end of this function so they never
  # depend on either the single sigmoid or the double sigmoid.
  width_cols_no_sigmoid <- intersect(
    c(
      "FWHP", "FWHP_Left_Time", "FWHP_Right_Time", "FWHP_Level",
      "FWHM", "FWHM_Left_Time", "FWHM_Right_Time", "FWHM_Level"
    ),
    names(metrics_out)
  )
  width_metrics_no_sigmoid <- if (length(width_cols_no_sigmoid) > 0) {
    metrics_out[, width_cols_no_sigmoid, drop = FALSE]
  } else {
    data.frame()
  }
  
  # ------------------------------------------------------------
  # Freeze prominence before optional sigmoid modeling.
  # ------------------------------------------------------------
  # User-requested behavior:
  # Prominence must ALWAYS be calculated from the FFT-smoothed left valley
  # of the original calcium transient, even when the single sigmoid and/or
  # double sigmoid options are activated.
  #
  # Therefore, the main columns below are saved before any sigmoid fit and
  # restored at the end of this function:
  #   Prominence
  #   Prominence_Reference
  #   Prominence_Reference_Time
  #   Left_Valley_FFT_Value
  #   Left_Valley_FFT_Time
  #
  # Double-sigmoid-specific prominence is still kept only as an optional
  # diagnostic column named DoubleSigmoid_Prominence, but it never replaces
  # the main Prominence metric.
  prominence_cols_no_sigmoid <- intersect(
    c(
      "Prominence",
      "Prominence_Reference",
      "Prominence_Reference_Time",
      "Left_Valley_FFT_Value",
      "Left_Valley_FFT_Time"
    ),
    names(metrics_out)
  )
  prominence_metrics_no_sigmoid <- if (length(prominence_cols_no_sigmoid) > 0) {
    metrics_out[, prominence_cols_no_sigmoid, drop = FALSE]
  } else {
    data.frame()
  }
  
  metrics_out$Sigmoid_Fit_Input <- "FFT-smoothed trace"
  metrics_out$Sigmoid_FFT_Fraction <- safe_num(fft_fraction, NA_real_)
  
  if (isTRUE(use_single_sigmoid)) {
    single_rows <- vector("list", nrow(metrics_out))
  }
  if (isTRUE(use_double_sigmoid)) {
    double_rows <- vector("list", nrow(metrics_out))
  }
  
  peak_times_all <- suppressWarnings(as.numeric(metrics_out$Peak_Time))
  
  for (i in seq_len(nrow(metrics_out))) {
    peak_time <- safe_num(metrics_out$Peak_Time[i], NA_real_)
    peak_value_fit <- if (is.finite(peak_time)) {
      stats::approx(time, signal, xout = peak_time, rule = 2)$y
    } else {
      NA_real_
    }
    
    fallback_baseline <- if ("Baseline" %in% names(metrics_out)) metrics_out$Baseline[i] else NA_real_
    
    prev_peaks <- peak_times_all[is.finite(peak_times_all) & peak_times_all < peak_time]
    prev_peak_time <- if (length(prev_peaks) > 0) max(prev_peaks, na.rm = TRUE) else NA_real_
    
    left_hint <- if ("Event_Baseline_Left_Time" %in% names(metrics_out) &&
                     is.finite(metrics_out$Event_Baseline_Left_Time[i])) {
      metrics_out$Event_Baseline_Left_Time[i]
    } else if ("Onset_Time" %in% names(metrics_out) &&
               is.finite(metrics_out$Onset_Time[i])) {
      metrics_out$Onset_Time[i]
    } else if ("Transient_Occurrence_Time" %in% names(metrics_out) &&
               is.finite(metrics_out$Transient_Occurrence_Time[i])) {
      metrics_out$Transient_Occurrence_Time[i]
    } else {
      NA_real_
    }
    
    left_min_info <- wavelet_find_adjacent_left_minimum(
      time = time,
      signal = signal,
      peak_time = peak_time,
      prev_peak_time = prev_peak_time
    )
    
    left_ref <- if (is.finite(left_min_info$time)) left_min_info$time else left_hint
    left_ref_value <- if (is.finite(left_min_info$value)) left_min_info$value else if (is.finite(left_ref)) stats::approx(time, signal, xout = left_ref, rule = 2)$y else NA_real_
    
    if (isTRUE(use_single_sigmoid)) {
      single_fit <- fit_wavelet_single_sigmoid_event(
        time = time,
        signal = signal,
        left_reference_time = left_ref,
        peak_time = peak_time,
        peak_value = peak_value_fit,
        fallback_baseline = fallback_baseline,
        pre_event_window = pre_event_window,
        fit_left_pad = single_left_pad,
        fit_right_pad = single_right_pad,
        left_reference_value = left_ref_value
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
        single_fit$fit_df$Event_ID <- metrics_out$Event_ID[i]
        curve_list[[length(curve_list) + 1]] <- single_fit$fit_df
      }
    }
    
    if (isTRUE(use_double_sigmoid)) {
      future_peaks <- peak_times_all[is.finite(peak_times_all) & peak_times_all > peak_time]
      next_peak_time <- if (length(future_peaks) > 0) min(future_peaks, na.rm = TRUE) else NA_real_
      
      min_right_candidates <- c(peak_time)
      if ("Offset_Time" %in% names(metrics_out) && is.finite(metrics_out$Offset_Time[i])) {
        min_right_candidates <- c(min_right_candidates, metrics_out$Offset_Time[i])
      }
      if ("FWHM_Right_Time" %in% names(metrics_out) && is.finite(metrics_out$FWHM_Right_Time[i])) {
        min_right_candidates <- c(min_right_candidates, metrics_out$FWHM_Right_Time[i])
      }
      if ("FWHP_Right_Time" %in% names(metrics_out) && is.finite(metrics_out$FWHP_Right_Time[i])) {
        min_right_candidates <- c(min_right_candidates, metrics_out$FWHP_Right_Time[i])
      }
      min_right_time <- max(min_right_candidates, na.rm = TRUE)
      
      right_ref_auto <- wavelet_find_double_sigmoid_decay_end(
        time = time,
        signal = signal,
        peak_time = peak_time,
        left_reference_time = left_ref,
        fallback_baseline = fallback_baseline,
        pre_event_window = pre_event_window,
        next_peak_time = next_peak_time,
        min_right_time = min_right_time
      )
      
      right_ref_fallback <- if ("Event_Baseline_Right_Time" %in% names(metrics_out) &&
                                is.finite(metrics_out$Event_Baseline_Right_Time[i])) {
        metrics_out$Event_Baseline_Right_Time[i]
      } else if ("Offset_Time" %in% names(metrics_out) &&
                 is.finite(metrics_out$Offset_Time[i])) {
        metrics_out$Offset_Time[i]
      } else if ("FWHM_Right_Time" %in% names(metrics_out) &&
                 is.finite(metrics_out$FWHM_Right_Time[i])) {
        metrics_out$FWHM_Right_Time[i]
      } else if ("FWHP_Right_Time" %in% names(metrics_out) &&
                 is.finite(metrics_out$FWHP_Right_Time[i])) {
        metrics_out$FWHP_Right_Time[i]
      } else {
        NA_real_
      }
      
      right_ref <- if (is.finite(right_ref_auto) && right_ref_auto > peak_time) {
        right_ref_auto
      } else {
        right_ref_fallback
      }
      right_ref_value <- if (is.finite(right_ref)) stats::approx(time, signal, xout = right_ref, rule = 2)$y else NA_real_
      
      double_fit <- fit_wavelet_double_sigmoid_event(
        time = time,
        signal = signal,
        left_reference_time = left_ref,
        peak_time = peak_time,
        right_reference_time = right_ref,
        peak_value = peak_value_fit,
        fallback_baseline = fallback_baseline,
        pre_event_window = pre_event_window,
        fit_left_pad = double_left_pad,
        fit_right_pad = double_right_pad,
        left_reference_value = left_ref_value,
        right_reference_value = right_ref_value
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
        double_fit$fit_df$Event_ID <- metrics_out$Event_ID[i]
        curve_list[[length(curve_list) + 1]] <- double_fit$fit_df
      }
    }
  }
  
  if (isTRUE(use_single_sigmoid)) {
    single_df <- do.call(rbind, single_rows)
    metrics_out <- cbind(metrics_out, single_df)
    
    valid_single <- is.finite(metrics_out$SingleSigmoid_Onset_Time) &
      metrics_out$SingleSigmoid_Fit_Status == "OK"
    
    metrics_out$Original_Onset_Time <- metrics_out$Onset_Time
    metrics_out$Original_Transient_Occurrence_Time <- metrics_out$Transient_Occurrence_Time
    metrics_out$Original_Peak_Rise_Time <- metrics_out$Peak_Rise_Time
    metrics_out$Original_Rise_Rate <- metrics_out$Rise_Rate
    
    metrics_out$Onset_Method <- ifelse(
      valid_single,
      "Single sigmoid + RMS after FFT smoothing",
      "Original wavelet local-minimum onset"
    )
    
    metrics_out$Onset_Time[valid_single] <- metrics_out$SingleSigmoid_Onset_Time[valid_single]
    metrics_out$Transient_Occurrence_Time[valid_single] <- metrics_out$SingleSigmoid_Onset_Time[valid_single]
    metrics_out$Onset[valid_single] <- metrics_out$SingleSigmoid_Onset_Time[valid_single]
    
    metrics_out$Peak_Rise_Time[valid_single] <-
      metrics_out$Peak_Time[valid_single] - metrics_out$Onset_Time[valid_single]
    
    if ("Amplitude_Above_Baseline" %in% names(metrics_out)) {
      metrics_out$Rise_Rate[valid_single] <- ifelse(
        is.finite(metrics_out$Peak_Rise_Time[valid_single]) &
          metrics_out$Peak_Rise_Time[valid_single] > 0,
        metrics_out$Amplitude_Above_Baseline[valid_single] /
          metrics_out$Peak_Rise_Time[valid_single],
        NA_real_
      )
    }
    
    metrics_out$SingleSigmoid_Rise_Time[valid_single] <- metrics_out$Peak_Rise_Time[valid_single]
    metrics_out$SingleSigmoid_Rise_Rate[valid_single] <- metrics_out$Rise_Rate[valid_single]
  }
  
  if (isTRUE(use_double_sigmoid)) {
    double_df <- do.call(rbind, double_rows)
    metrics_out <- cbind(metrics_out, double_df)
    
    valid_double <- metrics_out$DoubleSigmoid_Fit_Status == "OK" &
      is.finite(metrics_out$DoubleSigmoid_Left_Min_Time) &
      is.finite(metrics_out$DoubleSigmoid_Left_Min_Value) &
      is.finite(metrics_out$DoubleSigmoid_Right_Min_Time) &
      is.finite(metrics_out$DoubleSigmoid_Right_Min_Value) &
      is.finite(metrics_out$DoubleSigmoid_Prominence)
    
    # Keep copies of the original FFT-based prominence for transparency.
    metrics_out$Original_Prominence <- prominence_metrics_no_sigmoid$Prominence
    metrics_out$Original_Prominence_Reference <- prominence_metrics_no_sigmoid$Prominence_Reference
    metrics_out$Original_Prominence_Reference_Time <- prominence_metrics_no_sigmoid$Prominence_Reference_Time
    
    # IMPORTANT:
    # Do NOT replace the main Prominence, Prominence_Reference,
    # Prominence_Reference_Time, Left_Valley_FFT_Value or Left_Valley_FFT_Time
    # with double-sigmoid minima.
    #
    # The fitted double sigmoid can still provide its own optional diagnostic
    # columns, including DoubleSigmoid_Prominence and DoubleSigmoid_*_Min_*,
    # but the main Prominence metric must remain tied to the FFT-smoothed
    # minimum of the calcium transient.
    if (any(valid_double)) {
      metrics_out$DoubleSigmoid_Prominence_Method <-
        "Diagnostic only; main Prominence remains FFT-smoothed left-valley prominence"
      
      # Keep event-specific baseline geometry unchanged here as well. The
      # event-specific baseline line and AUC should remain based on the
      # left/right FFT minima of the calcium transient, not on sigmoid minima.
      event_specific_rows <- rep(FALSE, nrow(metrics_out))
      if ("Baseline_Type" %in% names(metrics_out)) {
        event_specific_rows <- valid_double & metrics_out$Baseline_Type == "event_specific"
      }
      if (any(event_specific_rows)) {
        metrics_out$Baseline_Method_For_Event[event_specific_rows] <-
          "Event-specific baseline: left-right FFT minima from calcium transient"
      }
    }
    
    metrics_out$Prominence_Method <-
      "FFT-smoothed left-valley prominence; independent of single and double sigmoid"
  }
  
  # Final safety restoration: prominence stays locked to the original
  # FFT-smoothed left-valley calculation even when single/double sigmoid
  # options are active.
  if (length(prominence_cols_no_sigmoid) > 0 &&
      nrow(metrics_out) == nrow(prominence_metrics_no_sigmoid)) {
    metrics_out[, prominence_cols_no_sigmoid] <-
      prominence_metrics_no_sigmoid[, prominence_cols_no_sigmoid, drop = FALSE]
  }
  
  if (!"Prominence_Method" %in% names(metrics_out)) {
    metrics_out$Prominence_Method <-
      "FFT-smoothed left-valley prominence; independent of single and double sigmoid"
  }
  metrics_out$Prominence_Method <-
    "FFT-smoothed left-valley prominence; independent of single and double sigmoid"
  
  # Final safety restoration: FWHP and FWHM stay locked to the original
  # non-sigmoid geometry even when single/double sigmoid options are active.
  if (length(width_cols_no_sigmoid) > 0 &&
      nrow(metrics_out) == nrow(width_metrics_no_sigmoid)) {
    metrics_out[, width_cols_no_sigmoid] <- width_metrics_no_sigmoid[, width_cols_no_sigmoid, drop = FALSE]
  }
  
  fitted_curves <- if (length(curve_list) > 0) do.call(rbind, curve_list) else data.frame()
  
  list(
    metrics = metrics_out,
    fitted_curves = fitted_curves
  )
}

update_wavelet_event_auc_from_metrics <- function(metrics, time, signal) {
  if (is.null(metrics) || nrow(metrics) == 0) return(metrics)
  
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  
  if (!all(c("Onset_Time", "Offset_Time", "Baseline") %in% names(metrics))) {
    return(metrics)
  }
  
  event_auc <- vapply(seq_len(nrow(metrics)), function(i) {
    onset_time <- safe_num(metrics$Onset_Time[i], NA_real_)
    offset_time <- safe_num(metrics$Offset_Time[i], NA_real_)
    if (!is.finite(onset_time) || !is.finite(offset_time) || offset_time <= onset_time) return(NA_real_)
    
    idx <- which(time >= onset_time & time <= offset_time)
    if (length(idx) < 2) return(NA_real_)
    
    if (all(c(
      "Event_Baseline_Left_Time",
      "Event_Baseline_Left_Value",
      "Event_Baseline_Right_Time",
      "Event_Baseline_Right_Value"
    ) %in% names(metrics)) &&
    is.finite(metrics$Event_Baseline_Left_Time[i]) &&
    is.finite(metrics$Event_Baseline_Right_Time[i]) &&
    metrics$Event_Baseline_Right_Time[i] > metrics$Event_Baseline_Left_Time[i]) {
      baseline_vec <- stats::approx(
        x = c(metrics$Event_Baseline_Left_Time[i], metrics$Event_Baseline_Right_Time[i]),
        y = c(metrics$Event_Baseline_Left_Value[i], metrics$Event_Baseline_Right_Value[i]),
        xout = time[idx],
        rule = 2
      )$y
    } else {
      baseline_vec <- rep(metrics$Baseline[i], length(idx))
    }
    
    trapezoid_auc_above_baseline(
      time = time[idx],
      signal = signal[idx],
      baseline = baseline_vec
    )
  }, numeric(1))
  
  metrics$Event_AUC <- event_auc
  metrics
}


# ------------------------------------------------------------
# Keep the trace-start-to-onset baseline synchronized after optional sigmoids
# ------------------------------------------------------------
# For the baseline "Constant baseline: mean from trace start to Time_Onset",
# the selected baseline shown in the plot must be the same value stored in
# every event row. When the single-sigmoid option is active, Time_Onset can be
# refined after the first baseline calculation. This helper recalculates the
# scalar baseline using the final onset column and updates all baseline-
# dependent event columns so the blue amplitude segment starts exactly at the
# black dashed selected-baseline line.
sync_wavelet_onset_mean_baseline_with_metrics <- function(metrics, time, signal) {
  if (is.null(metrics) || nrow(metrics) == 0) return(metrics)
  
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  
  onset_values <- numeric(0)
  if ("Transient_Occurrence_Time" %in% names(metrics)) {
    onset_values <- suppressWarnings(as.numeric(metrics$Transient_Occurrence_Time))
  } else if ("Onset_Time" %in% names(metrics)) {
    onset_values <- suppressWarnings(as.numeric(metrics$Onset_Time))
  }
  onset_values <- onset_values[is.finite(onset_values)]
  
  time_onset <- if (length(onset_values) > 0) min(onset_values, na.rm = TRUE) else NA_real_
  baseline_value <- compute_wavelet_onset_mean_baseline(
    time = time,
    signal = signal,
    time_onset = time_onset
  )
  
  if (!is.finite(baseline_value)) return(metrics)
  
  # Store the same selected-baseline value in every transient row.
  metrics$Baseline <- baseline_value
  metrics$Baseline_Method_For_Event <- "Constant baseline: mean from trace start to Time_Onset"
  metrics$Baseline_Type <- "static"
  metrics$Time_Onset_For_Baseline <- time_onset
  metrics$Baseline_Sync_Method <- "Synchronized after optional sigmoid onset refinement"
  
  # For a static baseline, the local baseline segment used by Event_AUC is
  # horizontal and should match the selected baseline line.
  if ("Onset_Time" %in% names(metrics)) {
    metrics$Event_Baseline_Left_Time <- suppressWarnings(as.numeric(metrics$Onset_Time))
  }
  if ("Offset_Time" %in% names(metrics)) {
    metrics$Event_Baseline_Right_Time <- suppressWarnings(as.numeric(metrics$Offset_Time))
  }
  metrics$Event_Baseline_Left_Value <- baseline_value
  metrics$Event_Baseline_Right_Value <- baseline_value
  
  # Update amplitude-dependent quantities so the blue vertical line and the
  # exported tables use the same baseline value as the dashed black line.
  if ("Peak_Amplitude" %in% names(metrics)) {
    peak_amp <- suppressWarnings(as.numeric(metrics$Peak_Amplitude))
    metrics$Amplitude_Above_Baseline <- peak_amp - baseline_value
    metrics$DeltaF_over_F0_Peak <- ifelse(
      is.finite(baseline_value) && baseline_value != 0,
      metrics$Amplitude_Above_Baseline / baseline_value,
      NA_real_
    )
  }
  
  # Rise rate depends on amplitude and rise time. Width metrics are intentionally
  # not recomputed here, preserving the previous requirement that FWHM/FWHP stay
  # independent of optional single/double sigmoid modeling.
  if (all(c("Amplitude_Above_Baseline", "Peak_Rise_Time") %in% names(metrics))) {
    rise_time <- suppressWarnings(as.numeric(metrics$Peak_Rise_Time))
    amp <- suppressWarnings(as.numeric(metrics$Amplitude_Above_Baseline))
    metrics$Rise_Rate <- ifelse(
      is.finite(rise_time) & rise_time > 0,
      amp / rise_time,
      NA_real_
    )
  }
  
  if ("SingleSigmoid_Rise_Rate" %in% names(metrics) && "Rise_Rate" %in% names(metrics)) {
    valid_single <- "SingleSigmoid_Fit_Status" %in% names(metrics) &
      metrics$SingleSigmoid_Fit_Status == "OK" &
      is.finite(suppressWarnings(as.numeric(metrics$Rise_Rate)))
    metrics$SingleSigmoid_Rise_Rate[valid_single] <- metrics$Rise_Rate[valid_single]
  }
  
  metrics
}

# Real Mexican-hat/Ricker wavelet approximation.
# This avoids complex-valued rendering issues in Shiny.
ricker_wavelet <- function(u) {
  (1 - u^2) * exp(-(u^2) / 2)
}

# ------------------------------------------------------------
# FFT smoothing used ONLY for prominence calculation
# ------------------------------------------------------------
fft_smooth_for_prominence <- function(signal, f = event_specific_fft_f, keep_mean = TRUE) {
  signal <- suppressWarnings(as.numeric(signal))
  n <- length(signal)
  
  if (n < 4 || all(!is.finite(signal))) return(signal)
  
  # Replace missing values by linear interpolation if needed.
  if (any(!is.finite(signal))) {
    idx <- seq_along(signal)
    good <- is.finite(signal)
    signal[!good] <- stats::approx(
      x = idx[good],
      y = signal[good],
      xout = idx[!good],
      rule = 2
    )$y
  }
  
  f <- safe_num(f, 0.20)
  f <- max(0.01, min(1, f))
  
  x_fft <- stats::fft(signal)
  
  # Keep a fraction f of the lowest-frequency coefficients.
  # The DC component is coefficient 1.
  k_keep <- max(1, floor((n * f) / 2))
  
  keep <- rep(FALSE, n)
  keep[1] <- keep_mean
  
  # Positive low frequencies
  keep[2:min(n, k_keep + 1)] <- TRUE
  
  # Negative low frequencies at the end of the FFT vector
  if (k_keep > 0) {
    keep[(n - k_keep + 1):n] <- TRUE
  }
  
  x_fft_filtered <- x_fft
  x_fft_filtered[!keep] <- 0
  
  smoothed <- Re(stats::fft(x_fft_filtered, inverse = TRUE) / n)
  as.numeric(smoothed)
}

find_left_valley_fft_before_peak <- function(time, signal, peak_idx, f = event_specific_fft_f) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  
  n <- length(signal)
  peak_idx <- suppressWarnings(as.integer(peak_idx))
  
  if (!is.finite(peak_idx) || peak_idx <= 1 || peak_idx > n) {
    return(list(
      valley_idx = NA_integer_,
      valley_time = NA_real_,
      valley_value = NA_real_,
      smoothed_signal = fft_smooth_for_prominence(signal, f = f)
    ))
  }
  
  smoothed <- fft_smooth_for_prominence(signal, f = f)
  
  # Search the complete region before the peak, not only the event window.
  y <- smoothed[1:peak_idx]
  
  if (length(y) < 3) {
    valley_idx <- which.min(y)
  } else {
    local_min_idx <- which(
      y[2:(length(y) - 1)] <= y[1:(length(y) - 2)] &
        y[2:(length(y) - 1)] <= y[3:length(y)]
    ) + 1
    
    if (length(local_min_idx) > 0) {
      # Use the closest local minimum before the peak.
      valley_idx <- tail(local_min_idx, 1)
    } else {
      # Fallback: global minimum before the peak.
      valley_idx <- which.min(y)
    }
  }
  
  list(
    valley_idx = valley_idx,
    valley_time = time[valley_idx],
    valley_value = smoothed[valley_idx],
    smoothed_signal = smoothed
  )
}


find_right_valley_fft_after_peak <- function(time, signal, peak_idx, f = event_specific_fft_f) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  
  n <- length(signal)
  peak_idx <- suppressWarnings(as.integer(peak_idx))
  
  if (!is.finite(peak_idx) || peak_idx < 1 || peak_idx >= n) {
    return(list(
      valley_idx = NA_integer_,
      valley_time = NA_real_,
      valley_value = NA_real_,
      smoothed_signal = fft_smooth_for_prominence(signal, f = f)
    ))
  }
  
  smoothed <- fft_smooth_for_prominence(signal, f = f)
  
  # Search the complete region after the peak, not only the event window.
  y <- smoothed[peak_idx:n]
  
  if (length(y) < 3) {
    local_rel_idx <- which.min(y)
  } else {
    local_min_idx <- which(
      y[2:(length(y) - 1)] <= y[1:(length(y) - 2)] &
        y[2:(length(y) - 1)] <= y[3:length(y)]
    ) + 1
    
    if (length(local_min_idx) > 0) {
      # Use the closest local minimum after the peak.
      local_rel_idx <- local_min_idx[1]
    } else {
      # Fallback: global minimum after the peak.
      local_rel_idx <- which.min(y)
    }
  }
  
  valley_idx <- peak_idx + local_rel_idx - 1
  
  list(
    valley_idx = valley_idx,
    valley_time = time[valley_idx],
    valley_value = smoothed[valley_idx],
    smoothed_signal = smoothed
  )
}

# ------------------------------------------------------------
# Continuous wavelet transform approximation
# ------------------------------------------------------------
compute_cwt_ricker <- function(time, signal, scales) {
  time <- suppressWarnings(as.numeric(time))
  signal_std <- standardize_signal(signal)
  n <- length(signal_std)
  
  dt <- if (n > 1) stats::median(diff(time), na.rm = TRUE) else 1
  if (!is.finite(dt) || dt <= 0) dt <- 1
  
  scales <- suppressWarnings(as.numeric(scales))
  scales <- scales[is.finite(scales) & scales > 0]
  if (length(scales) == 0) scales <- seq(2, 40, length.out = 32)
  
  coef <- matrix(0, nrow = length(scales), ncol = n)
  
  for (si in seq_along(scales)) {
    s <- scales[si]
    for (ti in seq_len(n)) {
      u <- (time - time[ti]) / s
      psi <- ricker_wavelet(u) / sqrt(s)
      coef[si, ti] <- sum(signal_std * psi, na.rm = TRUE) * dt
    }
  }
  
  amp <- abs(coef)
  
  list(
    coef = coef,
    amplitude = amp,
    scales = scales,
    time = time
  )
}

# ------------------------------------------------------------
# Local maxima in the CWT map
# ------------------------------------------------------------
find_local_maxima_by_scale <- function(amp, time, scales, min_coef_quantile = 0.90) {
  if (is.null(amp) || length(amp) == 0 || ncol(amp) < 3) {
    return(data.frame(
      scale_index = integer(0),
      scale = numeric(0),
      time_index = integer(0),
      Time = numeric(0),
      coef_amp = numeric(0)
    ))
  }
  
  min_coef_quantile <- safe_num(min_coef_quantile, 0.90)
  min_coef_quantile <- max(0, min(1, min_coef_quantile))
  
  threshold <- as.numeric(stats::quantile(as.vector(amp), probs = min_coef_quantile, na.rm = TRUE))
  maxima <- list()
  
  for (si in seq_len(nrow(amp))) {
    y <- amp[si, ]
    if (length(y) < 3 || all(!is.finite(y))) next
    
    center <- y[2:(length(y) - 1)]
    left <- y[1:(length(y) - 2)]
    right <- y[3:length(y)]
    
    idx <- which(center > left & center >= right & center >= threshold) + 1
    
    if (length(idx) > 0) {
      maxima[[length(maxima) + 1]] <- data.frame(
        scale_index = si,
        scale = scales[si],
        time_index = idx,
        Time = time[idx],
        coef_amp = y[idx]
      )
    }
  }
  
  if (length(maxima) == 0) {
    return(data.frame(
      scale_index = integer(0),
      scale = numeric(0),
      time_index = integer(0),
      Time = numeric(0),
      coef_amp = numeric(0)
    ))
  }
  
  do.call(rbind, maxima)
}

# ------------------------------------------------------------
# Ridge construction
# ------------------------------------------------------------
build_wavelet_ridges <- function(maxima_df, time_tolerance = 2, min_ridge_length = 5) {
  if (is.null(maxima_df) || nrow(maxima_df) == 0) {
    return(data.frame(
      ridge_id = integer(0),
      ridge_length = integer(0),
      event_time = numeric(0),
      event_time_index = integer(0),
      max_coef_amp = numeric(0),
      mean_scale = numeric(0),
      min_scale = numeric(0),
      max_scale = numeric(0)
    ))
  }
  
  time_tolerance <- safe_num(time_tolerance, 2)
  min_ridge_length <- safe_num(min_ridge_length, 5)
  if (!is.finite(time_tolerance) || time_tolerance < 0) time_tolerance <- 2
  if (!is.finite(min_ridge_length) || min_ridge_length < 1) min_ridge_length <- 5
  
  maxima_df <- maxima_df[order(-maxima_df$scale_index, maxima_df$Time), , drop = FALSE]
  maxima_df$assigned <- FALSE
  ridges <- list()
  
  for (i in seq_len(nrow(maxima_df))) {
    if (isTRUE(maxima_df$assigned[i])) next
    
    current <- maxima_df[i, , drop = FALSE]
    maxima_df$assigned[i] <- TRUE
    last_time <- current$Time[1]
    last_scale_index <- current$scale_index[1]
    
    repeat {
      candidate_rows <- which(
        !maxima_df$assigned &
          maxima_df$scale_index == (last_scale_index - 1) &
          abs(maxima_df$Time - last_time) <= time_tolerance
      )
      
      if (length(candidate_rows) == 0) break
      
      best_row <- candidate_rows[which.max(maxima_df$coef_amp[candidate_rows])]
      current <- rbind(current, maxima_df[best_row, , drop = FALSE])
      maxima_df$assigned[best_row] <- TRUE
      last_time <- maxima_df$Time[best_row]
      last_scale_index <- maxima_df$scale_index[best_row]
    }
    
    current$ridge_id <- length(ridges) + 1L
    ridges[[length(ridges) + 1L]] <- current
  }
  
  if (length(ridges) == 0) {
    return(data.frame(
      ridge_id = integer(0),
      ridge_length = integer(0),
      event_time = numeric(0),
      event_time_index = integer(0),
      max_coef_amp = numeric(0),
      mean_scale = numeric(0),
      min_scale = numeric(0),
      max_scale = numeric(0)
    ))
  }
  
  ridges_df <- do.call(rbind, ridges)
  split_ridges <- split(ridges_df, ridges_df$ridge_id)
  
  ridge_summary <- do.call(rbind, lapply(split_ridges, function(z) {
    best_idx <- which.max(z$coef_amp)
    
    data.frame(
      ridge_id = z$ridge_id[1],
      ridge_length = nrow(z),
      event_time = z$Time[best_idx],
      event_time_index = z$time_index[best_idx],
      max_coef_amp = max(z$coef_amp, na.rm = TRUE),
      mean_scale = mean(z$scale, na.rm = TRUE),
      min_scale = min(z$scale, na.rm = TRUE),
      max_scale = max(z$scale, na.rm = TRUE)
    )
  }))
  
  ridge_summary <- ridge_summary[ridge_summary$ridge_length >= min_ridge_length, , drop = FALSE]
  if (nrow(ridge_summary) == 0) return(ridge_summary)
  ridge_summary[order(ridge_summary$event_time), , drop = FALSE]
}

merge_close_wavelet_events <- function(events, min_distance_time = 1) {
  if (is.null(events) || nrow(events) <= 1) return(events)
  
  min_distance_time <- safe_num(min_distance_time, 1)
  if (!is.finite(min_distance_time) || min_distance_time < 0) min_distance_time <- 1
  
  events <- events[order(events$event_time), , drop = FALSE]
  keep <- rep(TRUE, nrow(events))
  
  for (i in seq_len(nrow(events))) {
    if (!keep[i]) next
    
    close_idx <- which(
      keep &
        seq_len(nrow(events)) != i &
        abs(events$event_time - events$event_time[i]) <= min_distance_time
    )
    
    if (length(close_idx) > 0) {
      group_idx <- c(i, close_idx)
      best <- group_idx[which.max(events$max_coef_amp[group_idx])]
      keep[group_idx[group_idx != best]] <- FALSE
    }
  }
  
  events[keep, , drop = FALSE]
}

# ------------------------------------------------------------
# Signal geometry utilities for event metrics
# ------------------------------------------------------------
interp_crossing_time <- function(x1, y1, x2, y2, level) {
  if (!is.finite(y1) || !is.finite(y2) || !is.finite(level)) return(NA_real_)
  if (y1 == y2) return(x1)
  x1 + (level - y1) * (x2 - x1) / (y2 - y1)
}

find_left_crossing <- function(time, signal, peak_idx, level, left_limit_idx = 1) {
  if (!is.finite(level) || is.na(peak_idx) || peak_idx <= left_limit_idx) return(NA_real_)
  
  for (j in seq(from = peak_idx - 1, to = left_limit_idx, by = -1)) {
    y1 <- signal[j]
    y2 <- signal[j + 1]
    if ((y1 <= level && y2 >= level) || (y1 >= level && y2 <= level)) {
      return(interp_crossing_time(time[j], y1, time[j + 1], y2, level))
    }
  }
  
  NA_real_
}

find_right_crossing <- function(time, signal, peak_idx, level, right_limit_idx = length(time)) {
  if (!is.finite(level) || is.na(peak_idx) || peak_idx >= right_limit_idx) return(NA_real_)
  
  for (j in seq(from = peak_idx, to = right_limit_idx - 1, by = 1)) {
    y1 <- signal[j]
    y2 <- signal[j + 1]
    if ((y1 >= level && y2 <= level) || (y1 <= level && y2 >= level)) {
      return(interp_crossing_time(time[j], y1, time[j + 1], y2, level))
    }
  }
  
  NA_real_
}

trapezoid_auc_above_baseline <- function(time, signal, baseline) {
  if (length(time) < 2) return(NA_real_)
  
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  
  if (length(baseline) == 1) {
    baseline <- rep(baseline, length(signal))
  } else {
    baseline <- suppressWarnings(as.numeric(baseline))
  }
  
  if (length(baseline) != length(signal)) {
    baseline <- rep(stats::median(baseline, na.rm = TRUE), length(signal))
  }
  
  y <- pmax(signal - baseline, 0)
  sum(diff(time) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
}

rolling_quantile_baseline <- function(time, signal, window_size = 20, percentile = 2) {
  # Time-window rolling percentile baseline.
  # This matches the FFT + Baseline Analysis definition: for each time point,
  # use a centered moving window in the same time units as the Time column.
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  n <- length(signal)
  
  if (length(time) != n || n == 0) {
    return(rep(NA_real_, max(length(time), n)))
  }
  
  window_size <- safe_num(window_size, 20)
  if (!is.finite(window_size) || window_size <= 0) {
    window_size <- diff(range(time, na.rm = TRUE))
  }
  if (!is.finite(window_size) || window_size <= 0) {
    window_size <- 1
  }
  
  percentile <- safe_num(percentile, 2)
  percentile <- max(0, min(100, percentile))
  probs <- percentile / 100
  
  time_min <- min(time, na.rm = TRUE)
  time_max <- max(time, na.rm = TRUE)
  half_window <- window_size / 2
  
  out <- vapply(seq_len(n), function(i) {
    t0 <- time[i]
    xmin <- max(t0 - half_window, time_min)
    xmax <- min(t0 + half_window, time_max)
    idx <- which(is.finite(time) & is.finite(signal) & time >= xmin & time <= xmax)
    
    if (length(idx) == 0) {
      return(NA_real_)
    }
    
    as.numeric(stats::quantile(signal[idx], probs = probs, na.rm = TRUE, names = FALSE, type = 7))
  }, numeric(1))
  
  missing_idx <- which(!is.finite(out))
  if (length(missing_idx) > 0) {
    fallback <- as.numeric(stats::quantile(signal, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
    out[missing_idx] <- fallback
  }
  
  as.numeric(out)
}


estimate_wavelet_time_onset_from_events <- function(time,
                                                    signal,
                                                    events,
                                                    event_window = 2) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  
  if (is.null(events) || nrow(events) == 0 || length(time) < 3 || length(signal) != length(time)) {
    return(NA_real_)
  }
  
  event_window <- safe_num(event_window, 2)
  if (!is.finite(event_window) || event_window <= 0) event_window <- 2
  
  onset_times <- vapply(seq_len(nrow(events)), function(i) {
    center <- safe_num(events$event_time[i], NA_real_)
    if (!is.finite(center)) return(NA_real_)
    
    window_idx <- which(time >= center - event_window & time <= center + event_window)
    if (length(window_idx) < 3) {
      closest_idx <- which.min(abs(time - center))
      if (length(closest_idx) == 0 || !is.finite(closest_idx)) return(NA_real_)
      window_idx <- seq(
        max(1, closest_idx - 2),
        min(length(time), closest_idx + 2)
      )
    }
    
    if (length(window_idx) == 0 || all(!is.finite(signal[window_idx]))) return(NA_real_)
    
    peak_local <- which.max(signal[window_idx])
    peak_global_idx <- window_idx[peak_local]
    left_window_idx <- window_idx[window_idx <= peak_global_idx]
    
    onset_global_idx <- if (length(left_window_idx) > 0 && any(is.finite(signal[left_window_idx]))) {
      left_window_idx[which.min(signal[left_window_idx])]
    } else {
      peak_global_idx
    }
    
    time[onset_global_idx]
  }, numeric(1))
  
  onset_times <- onset_times[is.finite(onset_times)]
  if (length(onset_times) == 0) return(NA_real_)
  
  min(onset_times, na.rm = TRUE)
}

compute_wavelet_onset_mean_baseline <- function(time,
                                                signal,
                                                time_onset = NA_real_) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  time_onset <- safe_num(time_onset, NA_real_)
  
  safe_trace_mean <- function() {
    if (length(signal) == 0 || all(!is.finite(signal))) return(0)
    mean(signal, na.rm = TRUE)
  }
  
  if (length(time) != length(signal) || length(signal) == 0 || all(!is.finite(signal))) {
    return(0)
  }
  
  if (!is.finite(time_onset)) {
    return(safe_trace_mean())
  }
  
  idx <- which(is.finite(time) & is.finite(signal) & time <= time_onset)
  if (length(idx) == 0) {
    idx <- which(is.finite(signal))
    if (length(idx) > 10) idx <- head(idx, 10)
  }
  
  if (length(idx) == 0) return(safe_trace_mean())
  
  baseline <- mean(signal[idx], na.rm = TRUE)
  if (!is.finite(baseline)) baseline <- safe_trace_mean()
  baseline
}

compute_wavelet_baseline_info <- function(time,
                                          signal,
                                          baseline_method = "f0_initial",
                                          initial_fraction = input$initial_f0_fraction,
                                          minimal_fraction = 0.10,
                                          rolling_window = 20,
                                          rolling_percentile = 2,
                                          time_onset = NA_real_) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  
  baseline_method <- as.character(baseline_method)
  
  initial_fraction <- safe_num(initial_fraction, 0.10)
  initial_fraction <- max(0.01, min(1, initial_fraction))
  
  minimal_fraction <- safe_num(minimal_fraction, 0.10)
  minimal_fraction <- max(0.01, min(1, minimal_fraction))
  
  # Backward-compatible options from previous versions.
  if (baseline_method == "median") {
    scalar <- stats::median(signal, na.rm = TRUE)
    return(list(method = "Median signal", type = "static", scalar = scalar, trace = NULL))
  }
  
  if (baseline_method == "mean") {
    scalar <- mean(signal, na.rm = TRUE)
    return(list(method = "Mean signal", type = "static", scalar = scalar, trace = NULL))
  }
  
  if (baseline_method == "minimum") {
    scalar <- min(signal, na.rm = TRUE)
    return(list(method = "Min", type = "static", scalar = scalar, trace = NULL))
  }
  
  if (baseline_method == "f0_initial") {
    # Original standard baseline definition:
    # use an initial percentage/fraction of the trace and calculate its mean.
    n_initial <- max(3, floor(length(signal) * initial_fraction))
    n_initial <- min(n_initial, length(signal))
    
    scalar <- mean(signal[seq_len(n_initial)], na.rm = TRUE)
    
    return(list(
      method = "Standard definition",
      type = "static",
      scalar = scalar,
      trace = NULL,
      initial_fraction = initial_fraction
    ))
  }
  
  if (baseline_method == "f0_onset_mean") {
    # Constant baseline requested by the user:
    # calculate one scalar mean from the beginning of the calcium trace
    # up to the first detected Time_Onset, then use that same value for all events.
    time_onset <- safe_num(time_onset, NA_real_)
    scalar <- compute_wavelet_onset_mean_baseline(
      time = time,
      signal = signal,
      time_onset = time_onset
    )
    
    return(list(
      method = "Constant baseline: mean from trace start to Time_Onset",
      type = "static",
      scalar = scalar,
      trace = NULL,
      time_onset = time_onset
    ))
  }
  
  if (baseline_method == "f0_minimal") {
    # Match FFT + Baseline Analysis:
    # sort the finite fluorescence values and calculate the mean of the lowest fraction.
    y_sorted <- sort(signal[is.finite(signal)])
    if (length(y_sorted) == 0) {
      scalar <- mean(signal, na.rm = TRUE)
    } else {
      n_low <- max(1L, floor(length(y_sorted) * minimal_fraction))
      scalar <- mean(y_sorted[seq_len(n_low)], na.rm = TRUE)
    }
    return(list(
      method = "Low-fluorescence region baseline",
      type = "static",
      scalar = scalar,
      trace = NULL,
      minimal_fraction = minimal_fraction
    ))
  }
  
  if (baseline_method == "f0_smooth") {
    baseline_trace <- rolling_quantile_baseline(
      time = time,
      signal = signal,
      window_size = rolling_window,
      percentile = rolling_percentile
    )
    
    return(list(
      method = "Rolling percentile baseline",
      type = "dynamic",
      scalar = stats::median(baseline_trace, na.rm = TRUE),
      trace = data.frame(Time = time, Baseline = baseline_trace),
      rolling_window = rolling_window,
      rolling_percentile = rolling_percentile
    ))
  }
  
  if (baseline_method == "f0_event_specific") {
    return(list(
      method = "Event-specific baseline: left-right FFT minima",
      type = "event_specific",
      scalar = stats::median(signal, na.rm = TRUE),
      trace = NULL
    ))
  }
  
  scalar <- stats::median(signal, na.rm = TRUE)
  list(method = "Median signal", type = "static", scalar = scalar, trace = NULL)
}

baseline_at_times <- function(time_points, baseline_info) {
  time_points <- suppressWarnings(as.numeric(time_points))
  
  if (is.null(baseline_info$trace)) {
    return(rep(baseline_info$scalar, length(time_points)))
  }
  
  stats::approx(
    x = baseline_info$trace$Time,
    y = baseline_info$trace$Baseline,
    xout = time_points,
    rule = 2
  )$y
}

compute_local_event_f0 <- function(time,
                                   signal,
                                   onset_time,
                                   local_pre_event_window = 5,
                                   statistic = "median") {
  local_pre_event_window <- safe_num(local_pre_event_window, 5)
  if (!is.finite(local_pre_event_window) || local_pre_event_window <= 0) {
    local_pre_event_window <- 5
  }
  
  idx <- which(time >= (onset_time - local_pre_event_window) & time <= onset_time)
  
  if (length(idx) < 3) {
    idx <- which(time <= onset_time)
    if (length(idx) > 10) idx <- tail(idx, 10)
  }
  
  if (length(idx) == 0) return(stats::median(signal, na.rm = TRUE))
  
  statistic <- as.character(statistic)
  
  if (statistic == "p10") {
    return(as.numeric(stats::quantile(signal[idx], probs = 0.10, na.rm = TRUE, names = FALSE)))
  }
  
  stats::median(signal[idx], na.rm = TRUE)
}

compute_global_auc_above_baseline <- function(time, signal, baseline_info) {
  if (is.null(baseline_info)) return(NA_real_)
  
  # For Event-specific baseline, there is no single biological baseline curve for the
  # entire trace. For visualization and the green total-AUC background,
  # use the scalar fallback stored in baseline_info$scalar.
  if (!is.null(baseline_info$type) && baseline_info$type == "event_specific") {
    baseline_vec <- rep(baseline_info$scalar, length(signal))
  } else {
    baseline_vec <- baseline_at_times(time, baseline_info)
  }
  
  trapezoid_auc_above_baseline(time, signal, baseline_vec)
}

build_global_auc_plot_data <- function(time, signal, baseline_info) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  
  empty_df <- data.frame(
    Time = numeric(0),
    Signal = numeric(0),
    Baseline = numeric(0),
    AUC_Upper = numeric(0),
    Above_Baseline = logical(0)
  )
  
  if (is.null(baseline_info)) {
    return(empty_df)
  }
  
  # For Event-specific baseline, draw the total AUC background using a scalar
  # fallback baseline. The event-specific yellow areas still use the
  # left-right FFT-minima baseline for each transient.
  if (!is.null(baseline_info$type) && baseline_info$type == "event_specific") {
    baseline_vec <- rep(baseline_info$scalar, length(signal))
  } else {
    baseline_vec <- baseline_at_times(time, baseline_info)
  }
  
  data.frame(
    Time = time,
    Signal = signal,
    Baseline = baseline_vec,
    AUC_Upper = pmax(signal, baseline_vec),
    Above_Baseline = signal > baseline_vec
  )
}

build_event_auc_plot_data <- function(time, signal, metrics) {
  empty_df <- data.frame(
    Event_ID = integer(0),
    Time = numeric(0),
    Signal = numeric(0),
    Baseline = numeric(0),
    AUC_Upper = numeric(0),
    Event_AUC = numeric(0),
    Event_Status = character(0),
    stringsAsFactors = FALSE
  )
  
  if (is.null(metrics) || nrow(metrics) == 0) return(empty_df)
  
  auc_list <- lapply(seq_len(nrow(metrics)), function(i) {
    
    # The visible Event_AUC region is restricted to the detected transient:
    # onset -> offset. This avoids a very wide event area when the FFT-smoothed
    # left/right minima are far from the peak.
    plot_left_time <- metrics$Onset_Time[i]
    plot_right_time <- metrics$Offset_Time[i]
    
    # The baseline used inside that narrow event window can still be the
    # event-specific line estimated from FFT-smoothed left/right minima.
    if (all(c(
      "Event_Baseline_Left_Time",
      "Event_Baseline_Left_Value",
      "Event_Baseline_Right_Time",
      "Event_Baseline_Right_Value"
    ) %in% names(metrics)) &&
    is.finite(metrics$Event_Baseline_Left_Time[i]) &&
    is.finite(metrics$Event_Baseline_Right_Time[i]) &&
    metrics$Event_Baseline_Right_Time[i] > metrics$Event_Baseline_Left_Time[i]) {
      
      baseline_left_time <- metrics$Event_Baseline_Left_Time[i]
      baseline_right_time <- metrics$Event_Baseline_Right_Time[i]
      baseline_left_value <- metrics$Event_Baseline_Left_Value[i]
      baseline_right_value <- metrics$Event_Baseline_Right_Value[i]
      
    } else {
      baseline_left_time <- plot_left_time
      baseline_right_time <- plot_right_time
      baseline_left_value <- metrics$Baseline[i]
      baseline_right_value <- metrics$Baseline[i]
    }
    
    idx <- which(time >= plot_left_time & time <= plot_right_time)
    if (length(idx) < 2) return(NULL)
    
    baseline_vec <- stats::approx(
      x = c(baseline_left_time, baseline_right_time),
      y = c(baseline_left_value, baseline_right_value),
      xout = time[idx],
      rule = 2
    )$y
    
    data.frame(
      Event_ID = metrics$Event_ID[i],
      Time = time[idx],
      Signal = signal[idx],
      Baseline = baseline_vec,
      AUC_Upper = pmax(signal[idx], baseline_vec),
      Event_AUC = metrics$Event_AUC[i],
      Event_Status = metrics$Event_Status[i],
      stringsAsFactors = FALSE
    )
  })
  
  auc_list <- auc_list[!sapply(auc_list, is.null)]
  if (length(auc_list) == 0) return(empty_df)
  do.call(rbind, auc_list)
}


build_event_specific_f0_interval_plot_data <- function(time, signal, metrics) {
  empty_df <- data.frame(
    Event_ID = integer(0),
    Time = numeric(0),
    Signal = numeric(0),
    Baseline = numeric(0),
    AUC_Upper = numeric(0),
    Event_Status = character(0),
    stringsAsFactors = FALSE
  )
  
  if (is.null(metrics) || nrow(metrics) == 0) return(empty_df)
  
  required_cols <- c(
    "Event_Baseline_Left_Time",
    "Event_Baseline_Left_Value",
    "Event_Baseline_Right_Time",
    "Event_Baseline_Right_Value",
    "Peak_Time"
  )
  
  if (!all(required_cols %in% names(metrics))) return(empty_df)
  
  # Sort by peak time so each event gets its own local territory.
  metrics <- metrics[order(metrics$Peak_Time), , drop = FALSE]
  
  interval_list <- lapply(seq_len(nrow(metrics)), function(i) {
    left_time <- metrics$Event_Baseline_Left_Time[i]
    right_time <- metrics$Event_Baseline_Right_Time[i]
    left_value <- metrics$Event_Baseline_Left_Value[i]
    right_value <- metrics$Event_Baseline_Right_Value[i]
    
    if (!is.finite(left_time) ||
        !is.finite(right_time) ||
        !is.finite(left_value) ||
        !is.finite(right_value) ||
        right_time <= left_time) {
      return(NULL)
    }
    
    # ------------------------------------------------------------
    # Prevent the green event-specific baseline area from crossing into
    # neighboring transients.
    #
    # The FFT-smoothed right/left minima can be far away when the
    # signal between events is noisy or slowly varying. Therefore,
    # for visualization we clip each green interval using the midpoint
    # between neighboring peak times.
    # ------------------------------------------------------------
    current_peak <- metrics$Peak_Time[i]
    
    left_limit <- min(time, na.rm = TRUE)
    right_limit <- max(time, na.rm = TRUE)
    
    if (i > 1 && is.finite(metrics$Peak_Time[i - 1])) {
      left_limit <- (metrics$Peak_Time[i - 1] + current_peak) / 2
    }
    
    if (i < nrow(metrics) && is.finite(metrics$Peak_Time[i + 1])) {
      right_limit <- (current_peak + metrics$Peak_Time[i + 1]) / 2
    }
    
    plot_left_time <- max(left_time, left_limit, na.rm = TRUE)
    plot_right_time <- min(right_time, right_limit, na.rm = TRUE)
    
    # Fallback: if clipping makes the region invalid, use the detected
    # onset-offset interval if available.
    if (!is.finite(plot_left_time) ||
        !is.finite(plot_right_time) ||
        plot_right_time <= plot_left_time) {
      if (all(c("Onset_Time", "Offset_Time") %in% names(metrics)) &&
          is.finite(metrics$Onset_Time[i]) &&
          is.finite(metrics$Offset_Time[i]) &&
          metrics$Offset_Time[i] > metrics$Onset_Time[i]) {
        plot_left_time <- metrics$Onset_Time[i]
        plot_right_time <- metrics$Offset_Time[i]
      } else {
        return(NULL)
      }
    }
    
    idx <- which(time >= plot_left_time & time <= plot_right_time)
    if (length(idx) < 2) return(NULL)
    
    # Use the original left-right FFT-minima baseline to define the
    # baseline value inside the clipped interval.
    baseline_vec <- stats::approx(
      x = c(left_time, right_time),
      y = c(left_value, right_value),
      xout = time[idx],
      rule = 2
    )$y
    
    data.frame(
      Event_ID = metrics$Event_ID[i],
      Time = time[idx],
      Signal = signal[idx],
      Baseline = baseline_vec,
      AUC_Upper = pmax(signal[idx], baseline_vec),
      Event_Status = metrics$Event_Status[i],
      stringsAsFactors = FALSE
    )
  })
  
  interval_list <- interval_list[!sapply(interval_list, is.null)]
  if (length(interval_list) == 0) return(empty_df)
  do.call(rbind, interval_list)
}


build_baseline_line_plot_data <- function(time, baseline_info, metrics = NULL) {
  empty_df <- data.frame(
    Time = numeric(0),
    Baseline = numeric(0),
    Baseline_Group = character(0),
    stringsAsFactors = FALSE
  )
  
  if (is.null(baseline_info)) return(empty_df)
  
  # Event-specific baseline: draw one dashed baseline segment for each event.
  if (!is.null(baseline_info$type) && baseline_info$type == "event_specific") {
    required_cols <- c(
      "Event_ID",
      "Event_Baseline_Left_Time",
      "Event_Baseline_Left_Value",
      "Event_Baseline_Right_Time",
      "Event_Baseline_Right_Value"
    )
    
    if (is.null(metrics) || nrow(metrics) == 0 || !all(required_cols %in% names(metrics))) {
      return(empty_df)
    }
    
    seg_list <- lapply(seq_len(nrow(metrics)), function(i) {
      left_time <- metrics$Event_Baseline_Left_Time[i]
      right_time <- metrics$Event_Baseline_Right_Time[i]
      left_value <- metrics$Event_Baseline_Left_Value[i]
      right_value <- metrics$Event_Baseline_Right_Value[i]
      
      if (!is.finite(left_time) ||
          !is.finite(right_time) ||
          !is.finite(left_value) ||
          !is.finite(right_value) ||
          right_time <= left_time) {
        return(NULL)
      }
      
      data.frame(
        Time = c(left_time, right_time),
        Baseline = c(left_value, right_value),
        Baseline_Group = paste0("Event_", metrics$Event_ID[i]),
        stringsAsFactors = FALSE
      )
    })
    
    seg_list <- seg_list[!sapply(seg_list, is.null)]
    if (length(seg_list) == 0) return(empty_df)
    return(do.call(rbind, seg_list))
  }
  
  # Dynamic Rolling percentile baseline: draw the full time-varying baseline curve.
  if (!is.null(baseline_info$type) && baseline_info$type == "dynamic") {
    baseline_vec <- baseline_at_times(time, baseline_info)
    return(data.frame(
      Time = time,
      Baseline = baseline_vec,
      Baseline_Group = "Dynamic_Baseline",
      stringsAsFactors = FALSE
    ))
  }
  
  # Static baseline: draw a horizontal baseline across the full trace.
  baseline_vec <- rep(baseline_info$scalar, length(time))
  
  data.frame(
    Time = time,
    Baseline = baseline_vec,
    Baseline_Group = "Static_Baseline",
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# FFT + Baseline style tables for Wavelet module
# ------------------------------------------------------------
wavelet_fft_style_baseline_scalar <- function(metrics, baseline_info) {
  baseline_value <- NA_real_
  
  if (!is.null(baseline_info) && !is.null(baseline_info$type) &&
      baseline_info$type == "event_specific" &&
      !is.null(metrics) && nrow(metrics) > 0 &&
      "Baseline" %in% names(metrics)) {
    vals <- suppressWarnings(as.numeric(metrics$Baseline))
    vals <- vals[is.finite(vals)]
    if (length(vals) > 0) baseline_value <- mean(vals, na.rm = TRUE)
  }
  
  if (!is.finite(baseline_value) && !is.null(baseline_info) && !is.null(baseline_info$scalar)) {
    baseline_value <- safe_num(baseline_info$scalar, NA_real_)
  }
  
  baseline_value
}

wavelet_fft_style_auc_summary <- function(metrics, time, signal, baseline_info) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  
  empty <- data.frame(AUC = NA_real_, P_min = NA_real_, P_max = NA_real_)
  if (length(time) < 2 || length(signal) != length(time)) return(empty)
  
  # Event-specific baseline has no single full-trace baseline curve. For the
  # FFT-style Trace_Metrics table, report the sum of event AUC values and the
  # time span covered by the accepted event intervals.
  if (!is.null(baseline_info) && !is.null(baseline_info$type) && baseline_info$type == "event_specific") {
    if (!is.null(metrics) && nrow(metrics) > 0) {
      auc_val <- if ("Event_AUC" %in% names(metrics)) {
        sum(suppressWarnings(as.numeric(metrics$Event_AUC)), na.rm = TRUE)
      } else {
        NA_real_
      }
      onset <- if ("Onset_Time" %in% names(metrics)) {
        suppressWarnings(as.numeric(metrics$Onset_Time))
      } else if ("Transient_Occurrence_Time" %in% names(metrics)) {
        suppressWarnings(as.numeric(metrics$Transient_Occurrence_Time))
      } else {
        numeric(0)
      }
      offset <- if ("Offset_Time" %in% names(metrics)) {
        suppressWarnings(as.numeric(metrics$Offset_Time))
      } else {
        numeric(0)
      }
      onset <- onset[is.finite(onset)]
      offset <- offset[is.finite(offset)]
      return(data.frame(
        AUC = auc_val,
        P_min = if (length(onset) > 0) min(onset, na.rm = TRUE) else NA_real_,
        P_max = if (length(offset) > 0) max(offset, na.rm = TRUE) else NA_real_
      ))
    }
    return(data.frame(AUC = 0, P_min = NA_real_, P_max = NA_real_))
  }
  
  baseline_vec <- if (!is.null(baseline_info)) {
    baseline_at_times(time, baseline_info)
  } else {
    rep(NA_real_, length(time))
  }
  
  if (length(baseline_vec) != length(signal) || all(!is.finite(baseline_vec))) return(empty)
  
  positive_signal <- pmax(signal - baseline_vec, 0)
  auc_val <- if (length(positive_signal) >= 2) {
    sum(diff(time) * (head(positive_signal, -1) + tail(positive_signal, -1)) / 2, na.rm = TRUE)
  } else {
    NA_real_
  }
  above <- which(positive_signal > 0 & is.finite(positive_signal) & is.finite(time))
  
  data.frame(
    AUC = auc_val,
    P_min = if (length(above) > 0) min(time[above], na.rm = TRUE) else NA_real_,
    P_max = if (length(above) > 0) max(time[above], na.rm = TRUE) else NA_real_
  )
}

build_wavelet_transient_metrics_table <- function(metrics) {
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
  
  if (is.null(metrics) || nrow(metrics) == 0) {
    out <- as.data.frame(setNames(replicate(length(base_cols), numeric(0), simplify = FALSE), base_cols))
    return(out)
  }
  
  amplitude_above_baseline <- if ("Amplitude_Above_Baseline" %in% names(metrics)) {
    suppressWarnings(as.numeric(metrics$Amplitude_Above_Baseline))
  } else {
    suppressWarnings(as.numeric(metrics$Peak_Amplitude)) - suppressWarnings(as.numeric(metrics$Baseline))
  }
  
  transient_df <- data.frame(
    Amplitude = amplitude_above_baseline,
    Baseline_at_peak = metrics$Baseline,
    Peak_Occurence_Time = metrics$Peak_Time,
    Prominence = metrics$Prominence,
    FWHP = metrics$FWHP,
    FWHM = metrics$FWHM,
    Transient_Ocurrence_Time = metrics$Transient_Occurrence_Time,
    Peak_Rise_Time = metrics$Peak_Rise_Time,
    Rise_Rate = metrics$Rise_Rate,
    stringsAsFactors = FALSE
  )
  
  # Same compact optional sigmoid columns used by the FFT module downloads.
  # They appear only when the corresponding optional model was requested and
  # therefore exists in the metrics object.
  single_cols <- c(
    "SingleSigmoid_Slope_at_t0",
    "SingleSigmoid_R2"
  )
  
  double_cols <- c(
    "DoubleSigmoid_Decay_Time",
    "DoubleSigmoid_Decay_Slope"
  )
  
  optional_cols <- c(single_cols, double_cols)
  optional_cols <- optional_cols[optional_cols %in% names(metrics)]
  
  if (length(optional_cols) > 0) {
    transient_df <- cbind(transient_df, metrics[, optional_cols, drop = FALSE])
  }
  
  transient_df
}

build_wavelet_trace_metrics_table <- function(metrics, time, signal, baseline_info) {
  duration <- max(time, na.rm = TRUE) - min(time, na.rm = TRUE)
  baseline_value <- wavelet_fft_style_baseline_scalar(metrics, baseline_info)
  auc_df <- wavelet_fft_style_auc_summary(metrics, time, signal, baseline_info)
  
  if (is.null(metrics) || nrow(metrics) == 0) {
    trace_df <- data.frame(
      Time_Onset = NA_real_,
      Frequency = 0,
      Baseline = baseline_value,
      Number_of_Peaks = 0,
      stringsAsFactors = FALSE
    )
    trace_df <- cbind(trace_df, auc_df)
    return(trace_df)
  }
  
  onset_values <- if ("Transient_Occurrence_Time" %in% names(metrics)) {
    suppressWarnings(as.numeric(metrics$Transient_Occurrence_Time))
  } else if ("Onset_Time" %in% names(metrics)) {
    suppressWarnings(as.numeric(metrics$Onset_Time))
  } else {
    NA_real_
  }
  onset_values <- onset_values[is.finite(onset_values)]
  
  trace_df <- data.frame(
    Time_Onset = if (length(onset_values) > 0) min(onset_values, na.rm = TRUE) else NA_real_,
    Frequency = ifelse(is.finite(duration) && duration > 0, nrow(metrics) / duration, NA_real_),
    Baseline = baseline_value,
    Number_of_Peaks = nrow(metrics),
    stringsAsFactors = FALSE
  )
  
  # Same trace-level optional sigmoid summaries used by the FFT module.
  # Do not add Mean_SingleSigmoid_R2, Mean_DoubleSigmoid_Duration, or
  # Sigmoid_FFT_Fraction here because those are not part of the FFT downloads.
  if ("SingleSigmoid_Onset_Time" %in% names(metrics) && any(is.finite(metrics$SingleSigmoid_Onset_Time))) {
    trace_df$Mean_SingleSigmoid_Onset_Time <- mean(metrics$SingleSigmoid_Onset_Time, na.rm = TRUE)
    trace_df$Mean_SingleSigmoid_Rise_Time <- mean(metrics$SingleSigmoid_Rise_Time, na.rm = TRUE)
    trace_df$Mean_SingleSigmoid_Rise_Rate <- mean(metrics$SingleSigmoid_Rise_Rate, na.rm = TRUE)
    trace_df$Mean_SingleSigmoid_Slope_at_t0 <- mean(metrics$SingleSigmoid_Slope_at_t0, na.rm = TRUE)
  }
  
  if ("DoubleSigmoid_Decay_Time" %in% names(metrics) && any(is.finite(metrics$DoubleSigmoid_Decay_Time))) {
    trace_df$Mean_DoubleSigmoid_Decay_Time <- mean(metrics$DoubleSigmoid_Decay_Time, na.rm = TRUE)
  }
  
  if ("DoubleSigmoid_Decay_Slope" %in% names(metrics) && any(is.finite(metrics$DoubleSigmoid_Decay_Slope))) {
    trace_df$Mean_DoubleSigmoid_Decay_Slope <- mean(metrics$DoubleSigmoid_Decay_Slope, na.rm = TRUE)
  }
  
  # FFT-style AUC columns at the end of Trace_Metrics.
  trace_df <- cbind(trace_df, auc_df)
  
  trace_df
}

build_auc_summary_table <- function(metrics, time, signal, baseline_info) {
  duration <- max(time, na.rm = TRUE) - min(time, na.rm = TRUE)
  total_global_auc <- compute_global_auc_above_baseline(time, signal, baseline_info)
  
  baseline_display <- ifelse(
    !is.null(baseline_info$type) && baseline_info$type == "event_specific",
    "Event-specific baseline",
    round(baseline_info$scalar, 4)
  )
  
  if (is.null(metrics) || nrow(metrics) == 0) {
    return(data.frame(
      Events_Used = 0,
      Baseline_Method = baseline_info$method,
      Baseline = baseline_display,
      Total_AUC_Global_Above_Baseline = total_global_auc,
      Total_Event_AUC = 0,
      Mean_Event_AUC = NA_real_,
      Median_Event_AUC = NA_real_,
      Min_Event_AUC = NA_real_,
      Max_Event_AUC = NA_real_,
      Frequency = 0,
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    Events_Used = nrow(metrics),
    Baseline_Method = baseline_info$method,
    Baseline = baseline_display,
    Total_AUC_Global_Above_Baseline = total_global_auc,
    Total_Event_AUC = sum(metrics$Event_AUC, na.rm = TRUE),
    Mean_Event_AUC = mean(metrics$Event_AUC, na.rm = TRUE),
    Median_Event_AUC = stats::median(metrics$Event_AUC, na.rm = TRUE),
    Min_Event_AUC = min(metrics$Event_AUC, na.rm = TRUE),
    Max_Event_AUC = max(metrics$Event_AUC, na.rm = TRUE),
    Frequency = ifelse(is.finite(duration) && duration > 0, nrow(metrics) / duration, NA_real_),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# Remove duplicate wavelet events that resolve to the same calcium peak
# ------------------------------------------------------------
# Wavelet ridgewalking can occasionally detect two independent ridges that,
# after local peak refinement, point to the same Peak_Time in the calcium trace.
# In that case the graph shows one biological peak, but the metrics table can
# count two events. This function keeps only one row per Peak_Time.
#
# Priority for keeping a duplicate:
#   1) Accepted event over Low confidence / Rejected
#   2) Highest Max_CWT_Amplitude
#   3) Longest Ridge_Length
#   4) Highest Amplitude_Above_Baseline
#   5) Highest Event_AUC
#
# Event_ID is renumbered after duplicate removal so the UI and downloaded tables
# show consecutive biological events. Original_Event_ID is kept for diagnosis.
deduplicate_wavelet_metrics_by_peak_time <- function(metrics,
                                                     digits = 10,
                                                     renumber_event_id = TRUE) {
  if (is.null(metrics) || nrow(metrics) <= 1 || !"Peak_Time" %in% names(metrics)) {
    return(metrics)
  }
  
  peak_time <- suppressWarnings(as.numeric(metrics$Peak_Time))
  if (all(!is.finite(peak_time))) return(metrics)
  
  # Keep the original Event_ID for diagnostic purposes before renumbering.
  if ("Event_ID" %in% names(metrics) && !"Original_Event_ID" %in% names(metrics)) {
    metrics$Original_Event_ID <- metrics$Event_ID
  }
  
  # Use a rounded key to avoid tiny floating point differences, while still
  # treating duplicated discrete calcium time points as the same biological peak.
  peak_key <- ifelse(
    is.finite(peak_time),
    formatC(round(peak_time, digits = digits), digits = digits, format = "f"),
    paste0("NA_row_", seq_len(nrow(metrics)))
  )
  
  status_score <- rep(0, nrow(metrics))
  if ("Event_Status" %in% names(metrics)) {
    status_chr <- as.character(metrics$Event_Status)
    status_score <- dplyr::case_when(
      status_chr == "Accepted" ~ 3,
      grepl("^Low confidence", status_chr) ~ 2,
      grepl("^Rejected", status_chr) ~ 1,
      TRUE ~ 0
    )
  }
  
  get_numeric_col <- function(colname, default = NA_real_) {
    if (colname %in% names(metrics)) {
      z <- suppressWarnings(as.numeric(metrics[[colname]]))
    } else {
      z <- rep(default, nrow(metrics))
    }
    z[!is.finite(z)] <- -Inf
    z
  }
  
  max_cwt <- get_numeric_col("Max_CWT_Amplitude")
  ridge_length <- get_numeric_col("Ridge_Length")
  amplitude <- get_numeric_col("Amplitude_Above_Baseline")
  event_auc <- get_numeric_col("Event_AUC")
  
  keep_idx <- unlist(lapply(split(seq_len(nrow(metrics)), peak_key), function(idx) {
    if (length(idx) == 1) return(idx)
    score_df <- data.frame(
      row_index = idx,
      status_score = status_score[idx],
      max_cwt = max_cwt[idx],
      ridge_length = ridge_length[idx],
      amplitude = amplitude[idx],
      event_auc = event_auc[idx],
      peak_time = peak_time[idx]
    )
    score_df <- score_df[order(
      -score_df$status_score,
      -score_df$max_cwt,
      -score_df$ridge_length,
      -score_df$amplitude,
      -score_df$event_auc,
      score_df$row_index
    ), , drop = FALSE]
    score_df$row_index[1]
  }), use.names = FALSE)
  
  keep_idx <- sort(unique(keep_idx))
  out <- metrics[keep_idx, , drop = FALSE]
  
  # Sort by Peak_Time for readability in the plot and tables.
  if ("Peak_Time" %in% names(out)) {
    ord <- order(suppressWarnings(as.numeric(out$Peak_Time)), na.last = TRUE)
    out <- out[ord, , drop = FALSE]
  }
  
  if (isTRUE(renumber_event_id) && "Event_ID" %in% names(out)) {
    out$Event_ID <- seq_len(nrow(out))
  }
  
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------
# Wavelet-derived event metrics
# ------------------------------------------------------------
compute_wavelet_event_metrics <- function(time,
                                          signal,
                                          events,
                                          event_window = 2,
                                          baseline_method = "f0_initial",
                                          initial_fraction = input$initial_f0_fraction,
                                          minimal_fraction = 0.10,
                                          rolling_window = 20,
                                          rolling_percentile = 2,
                                          event_specific_fft_f = event_specific_fft_f,
                                          local_pre_event_window = 5,
                                          local_f0_statistic = "median",
                                          min_amp_above_baseline = 0,
                                          accepted_min_ridge_length = 7,
                                          min_cwt_amplitude = 0) {
  time <- suppressWarnings(as.numeric(time))
  signal <- suppressWarnings(as.numeric(signal))
  
  empty <- data.frame(
    Event_ID = integer(0),
    Wavelet_Event_Time = numeric(0),
    Onset_Time = numeric(0),
    Peak_Time = numeric(0),
    Offset_Time = numeric(0),
    Peak_Amplitude = numeric(0),
    Baseline = numeric(0),
    Baseline_Method_For_Event = character(0),
    Baseline_Type = character(0),
    Event_Baseline_Left_Time = numeric(0),
    Event_Baseline_Left_Value = numeric(0),
    Event_Baseline_Right_Time = numeric(0),
    Event_Baseline_Right_Value = numeric(0),
    Amplitude_Above_Baseline = numeric(0),
    DeltaF_over_F0_Peak = numeric(0),
    Prominence = numeric(0),
    Prominence_Reference = numeric(0),
    Prominence_Reference_Time = numeric(0),
    Left_Valley_FFT_Value = numeric(0),
    Left_Valley_FFT_Time = numeric(0),
    FWHP = numeric(0),
    FWHP_Left_Time = numeric(0),
    FWHP_Right_Time = numeric(0),
    FWHP_Level = numeric(0),
    FWHM = numeric(0),
    FWHM_Left_Time = numeric(0),
    FWHM_Right_Time = numeric(0),
    FWHM_Level = numeric(0),
    Peak_Rise_Time = numeric(0),
    Transient_Occurrence_Time = numeric(0),
    Onset = numeric(0),
    Rise_Rate = numeric(0),
    Event_AUC = numeric(0),
    Ridge_Length = integer(0),
    Max_CWT_Amplitude = numeric(0),
    Mean_Scale = numeric(0),
    Min_Scale = numeric(0),
    Max_Scale = numeric(0),
    Event_Status = character(0)
  )
  
  if (is.null(events) || nrow(events) == 0 || length(time) < 3) return(empty)
  
  event_window <- safe_num(event_window, 2)
  if (!is.finite(event_window) || event_window <= 0) event_window <- 2
  
  min_amp_above_baseline <- safe_num(min_amp_above_baseline, 0)
  accepted_min_ridge_length <- safe_num(accepted_min_ridge_length, 7)
  min_cwt_amplitude <- safe_num(min_cwt_amplitude, 0)
  event_specific_fft_f <- safe_num(event_specific_fft_f, 0.20)
  event_specific_fft_f <- max(0.01, min(1, event_specific_fft_f))
  
  time_onset_for_baseline <- if (identical(as.character(baseline_method), "f0_onset_mean")) {
    estimate_wavelet_time_onset_from_events(
      time = time,
      signal = signal,
      events = events,
      event_window = event_window
    )
  } else {
    NA_real_
  }
  
  baseline_info <- compute_wavelet_baseline_info(
    time = time,
    signal = signal,
    baseline_method = baseline_method,
    initial_fraction = initial_fraction,
    minimal_fraction = minimal_fraction,
    rolling_window = rolling_window,
    rolling_percentile = rolling_percentile,
    time_onset = time_onset_for_baseline
  )
  
  
  out <- lapply(seq_len(nrow(events)), function(i) {
    center <- events$event_time[i]
    
    window_idx <- which(time >= center - event_window & time <= center + event_window)
    if (length(window_idx) < 3) {
      closest_idx <- which.min(abs(time - center))
      window_idx <- seq(
        max(1, closest_idx - 2),
        min(length(time), closest_idx + 2)
      )
    }
    
    local_time <- time[window_idx]
    local_signal <- signal[window_idx]
    
    peak_local <- which.max(local_signal)
    peak_global_idx <- window_idx[peak_local]
    peak_time <- time[peak_global_idx]
    peak_amp <- signal[peak_global_idx]
    
    # Onset: local minimum before peak inside the event window.
    left_window_idx <- window_idx[window_idx <= peak_global_idx]
    if (length(left_window_idx) > 0) {
      onset_global_idx <- left_window_idx[which.min(signal[left_window_idx])]
    } else {
      onset_global_idx <- peak_global_idx
    }
    
    onset_time <- time[onset_global_idx]
    onset_amp <- signal[onset_global_idx]
    
    # Offset: local minimum after peak inside the event window.
    right_window_idx <- window_idx[window_idx >= peak_global_idx]
    if (length(right_window_idx) > 0) {
      offset_global_idx <- right_window_idx[which.min(signal[right_window_idx])]
    } else {
      offset_global_idx <- peak_global_idx
    }
    
    offset_time <- time[offset_global_idx]
    offset_amp <- signal[offset_global_idx]
    
    # ------------------------------------------------------------
    # Baseline after wavelet event detection
    # ------------------------------------------------------------
    # For event-specific baseline, the baseline is a local linear baseline
    # drawn from the FFT-smoothed minimum to the left of the peak
    # to the FFT-smoothed minimum to the right of the peak.
    left_event_valley_fft <- find_left_valley_fft_before_peak(
      time = time,
      signal = signal,
      peak_idx = peak_global_idx,
      f = event_specific_fft_f
    )
    
    right_event_valley_fft <- find_right_valley_fft_after_peak(
      time = time,
      signal = signal,
      peak_idx = peak_global_idx,
      f = event_specific_fft_f
    )
    
    if (!is.null(baseline_info$type) && baseline_info$type == "event_specific") {
      event_baseline_left_time <- left_event_valley_fft$valley_time
      event_baseline_left_value <- left_event_valley_fft$valley_value
      event_baseline_right_time <- right_event_valley_fft$valley_time
      event_baseline_right_value <- right_event_valley_fft$valley_value
      
      if (is.finite(event_baseline_left_time) &&
          is.finite(event_baseline_right_time) &&
          event_baseline_right_time > event_baseline_left_time) {
        
        baseline <- stats::approx(
          x = c(event_baseline_left_time, event_baseline_right_time),
          y = c(event_baseline_left_value, event_baseline_right_value),
          xout = peak_time,
          rule = 2
        )$y
        
      } else {
        baseline <- compute_local_event_f0(
          time = time,
          signal = signal,
          onset_time = onset_time,
          local_pre_event_window = local_pre_event_window,
          statistic = local_f0_statistic
        )
        
        event_baseline_left_time <- onset_time
        event_baseline_left_value <- baseline
        event_baseline_right_time <- offset_time
        event_baseline_right_value <- baseline
      }
      
    } else {
      baseline <- baseline_at_times(peak_time, baseline_info)[1]
      event_baseline_left_time <- onset_time
      event_baseline_left_value <- baseline
      event_baseline_right_time <- offset_time
      event_baseline_right_value <- baseline
    }
    
    amp_above_baseline <- peak_amp - baseline
    
    deltaf_over_f0_peak <- ifelse(
      is.finite(baseline) && baseline != 0,
      amp_above_baseline / baseline,
      NA_real_
    )
    
    # ------------------------------------------------------------
    # Prominence calculation requested by the user:
    # 1) Smooth the whole calcium trace using FFT with f = event_specific_fft_f.
    # 2) Search the left local minimum before the peak using the full trace.
    # 3) Prominence = Peak_Amplitude - Left_Valley_FFT_Value
    # ------------------------------------------------------------
    left_valley_fft <- left_event_valley_fft
    
    valley_reference <- left_valley_fft$valley_value
    valley_reference_time <- left_valley_fft$valley_time
    
    prominence <- peak_amp - valley_reference
    if (!is.finite(prominence)) prominence <- NA_real_
    
    # FWHM: full width at half maximum relative to event/global baseline.
    fwhm_level <- baseline + amp_above_baseline / 2
    
    # FWHP: full width at half prominence.
    fwhp_level <- peak_amp - prominence / 2
    
    # Search crossing points from the peak outward across the full trace.
    left_limit_idx <- 1
    right_limit_idx <- length(time)
    
    fwhm_left <- find_left_crossing(time, signal, peak_global_idx, fwhm_level, left_limit_idx)
    fwhm_right <- find_right_crossing(time, signal, peak_global_idx, fwhm_level, right_limit_idx)
    fwhm <- ifelse(is.finite(fwhm_left) && is.finite(fwhm_right), fwhm_right - fwhm_left, NA_real_)
    
    fwhp_left <- find_left_crossing(time, signal, peak_global_idx, fwhp_level, left_limit_idx)
    fwhp_right <- find_right_crossing(time, signal, peak_global_idx, fwhp_level, right_limit_idx)
    fwhp <- ifelse(is.finite(fwhp_left) && is.finite(fwhp_right), fwhp_right - fwhp_left, NA_real_)
    
    peak_rise_time <- peak_time - onset_time
    rise_rate <- ifelse(
      is.finite(peak_rise_time) && peak_rise_time > 0,
      amp_above_baseline / peak_rise_time,
      NA_real_
    )
    
    # Event_AUC is integrated only across the detected transient interval.
    # For event-specific baseline, the local baseline is still estimated from the
    # FFT-smoothed left/right minima, but the shaded/event AUC region is
    # restricted to onset -> offset so the event area stays narrow.
    auc_idx <- which(time >= onset_time & time <= offset_time)
    
    event_auc <- if (length(auc_idx) >= 2) {
      if (!is.null(baseline_info$type) && baseline_info$type == "dynamic") {
        baseline_vec <- baseline_at_times(time[auc_idx], baseline_info)
      } else if (!is.null(baseline_info$type) && baseline_info$type == "event_specific" &&
                 is.finite(event_baseline_left_time) &&
                 is.finite(event_baseline_right_time) &&
                 event_baseline_right_time > event_baseline_left_time) {
        baseline_vec <- stats::approx(
          x = c(event_baseline_left_time, event_baseline_right_time),
          y = c(event_baseline_left_value, event_baseline_right_value),
          xout = time[auc_idx],
          rule = 2
        )$y
      } else {
        baseline_vec <- rep(baseline, length(auc_idx))
      }
      
      trapezoid_auc_above_baseline(
        time = time[auc_idx],
        signal = signal[auc_idx],
        baseline = baseline_vec
      )
    } else {
      NA_real_
    }
    
    status <- dplyr::case_when(
      !is.finite(amp_above_baseline) ~ "Rejected: invalid amplitude",
      amp_above_baseline <= 0 ~ "Rejected: below baseline",
      amp_above_baseline < min_amp_above_baseline ~ "Low confidence: weak amplitude",
      events$ridge_length[i] < accepted_min_ridge_length ~ "Low confidence: short ridge",
      events$max_coef_amp[i] < min_cwt_amplitude ~ "Low confidence: weak CWT ridge",
      TRUE ~ "Accepted"
    )
    
    data.frame(
      Event_ID = i,
      Wavelet_Event_Time = center,
      Onset_Time = onset_time,
      Peak_Time = peak_time,
      Offset_Time = offset_time,
      Peak_Amplitude = peak_amp,
      Baseline = baseline,
      Baseline_Method_For_Event = baseline_info$method,
      Baseline_Type = baseline_info$type,
      Event_Baseline_Left_Time = event_baseline_left_time,
      Event_Baseline_Left_Value = event_baseline_left_value,
      Event_Baseline_Right_Time = event_baseline_right_time,
      Event_Baseline_Right_Value = event_baseline_right_value,
      Amplitude_Above_Baseline = amp_above_baseline,
      DeltaF_over_F0_Peak = deltaf_over_f0_peak,
      Prominence = prominence,
      Prominence_Reference = valley_reference,
      Prominence_Reference_Time = valley_reference_time,
      Left_Valley_FFT_Value = valley_reference,
      Left_Valley_FFT_Time = valley_reference_time,
      FWHP = fwhp,
      FWHP_Left_Time = fwhp_left,
      FWHP_Right_Time = fwhp_right,
      FWHP_Level = fwhp_level,
      FWHM = fwhm,
      FWHM_Left_Time = fwhm_left,
      FWHM_Right_Time = fwhm_right,
      FWHM_Level = fwhm_level,
      Peak_Rise_Time = peak_rise_time,
      Transient_Occurrence_Time = onset_time,
      Onset = onset_time,
      Rise_Rate = rise_rate,
      Event_AUC = event_auc,
      Ridge_Length = events$ridge_length[i],
      Max_CWT_Amplitude = events$max_coef_amp[i],
      Mean_Scale = events$mean_scale[i],
      Min_Scale = events$min_scale[i],
      Max_Scale = events$max_scale[i],
      Event_Status = status,
      stringsAsFactors = FALSE
    )
  })
  
  metrics_out <- do.call(rbind, out)
  metrics_out <- deduplicate_wavelet_metrics_by_peak_time(metrics_out)
  metrics_out
}

summarize_wavelet_metrics <- function(metrics, time, signal, baseline_info) {
  duration <- max(time, na.rm = TRUE) - min(time, na.rm = TRUE)
  total_global_auc <- compute_global_auc_above_baseline(time, signal, baseline_info)
  
  if (is.null(metrics) || nrow(metrics) == 0) {
    return(data.frame(
      Number_of_Wavelet_Events = 0,
      Frequency = 0,
      Baseline_Method = baseline_info$method,
      Mean_Amplitude_Above_Baseline = NA_real_,
      Mean_DeltaF_over_F0_Peak = NA_real_,
      Mean_Prominence = NA_real_,
      Mean_FWHP = NA_real_,
      Mean_FWHM = NA_real_,
      Mean_Peak_Rise_Time = NA_real_,
      Mean_Peak_Occurrence_Time = NA_real_,
      Mean_Transient_Occurrence_Time = NA_real_,
      Mean_Onset = NA_real_,
      Mean_Rise_Rate = NA_real_,
      Total_AUC_Global_Above_Baseline = total_global_auc,
      Total_Event_AUC = 0,
      Mean_Event_AUC = NA_real_,
      Mean_Ridge_Length = NA_real_,
      Mean_CWT_Amplitude = NA_real_,
      Mean_Scale = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  summary_df <- data.frame(
    Number_of_Wavelet_Events = nrow(metrics),
    Frequency = ifelse(is.finite(duration) && duration > 0, nrow(metrics) / duration, NA_real_),
    Baseline_Method = baseline_info$method,
    Mean_Amplitude_Above_Baseline = mean(metrics$Amplitude_Above_Baseline, na.rm = TRUE),
    Mean_DeltaF_over_F0_Peak = mean(metrics$DeltaF_over_F0_Peak, na.rm = TRUE),
    Mean_Prominence = mean(metrics$Prominence, na.rm = TRUE),
    Mean_FWHP = mean(metrics$FWHP, na.rm = TRUE),
    Mean_FWHM = mean(metrics$FWHM, na.rm = TRUE),
    Mean_Peak_Rise_Time = mean(metrics$Peak_Rise_Time, na.rm = TRUE),
    Mean_Peak_Occurrence_Time = mean(metrics$Peak_Time, na.rm = TRUE),
    Mean_Transient_Occurrence_Time = mean(metrics$Transient_Occurrence_Time, na.rm = TRUE),
    Mean_Onset = mean(metrics$Onset, na.rm = TRUE),
    Mean_Rise_Rate = mean(metrics$Rise_Rate, na.rm = TRUE),
    Total_AUC_Global_Above_Baseline = total_global_auc,
    Total_Event_AUC = sum(metrics$Event_AUC, na.rm = TRUE),
    Mean_Event_AUC = mean(metrics$Event_AUC, na.rm = TRUE),
    Mean_Ridge_Length = mean(metrics$Ridge_Length, na.rm = TRUE),
    Mean_CWT_Amplitude = mean(metrics$Max_CWT_Amplitude, na.rm = TRUE),
    Mean_Scale = mean(metrics$Mean_Scale, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  if ("SingleSigmoid_Onset_Time" %in% names(metrics) && any(is.finite(metrics$SingleSigmoid_Onset_Time))) {
    summary_df$Mean_SingleSigmoid_Onset_Time <- mean(metrics$SingleSigmoid_Onset_Time, na.rm = TRUE)
    summary_df$Mean_SingleSigmoid_Rise_Time <- mean(metrics$SingleSigmoid_Rise_Time, na.rm = TRUE)
    summary_df$Mean_SingleSigmoid_Rise_Rate <- mean(metrics$SingleSigmoid_Rise_Rate, na.rm = TRUE)
    summary_df$Mean_SingleSigmoid_Slope_at_t0 <- mean(metrics$SingleSigmoid_Slope_at_t0, na.rm = TRUE)
    summary_df$Mean_SingleSigmoid_R2 <- mean(metrics$SingleSigmoid_R2, na.rm = TRUE)
  }
  
  if ("DoubleSigmoid_Decay_Time" %in% names(metrics) && any(is.finite(metrics$DoubleSigmoid_Decay_Time))) {
    summary_df$Mean_DoubleSigmoid_Decay_Time <- mean(metrics$DoubleSigmoid_Decay_Time, na.rm = TRUE)
  }
  
  if ("DoubleSigmoid_Decay_Slope" %in% names(metrics) && any(is.finite(metrics$DoubleSigmoid_Decay_Slope))) {
    summary_df$Mean_DoubleSigmoid_Decay_Slope <- mean(metrics$DoubleSigmoid_Decay_Slope, na.rm = TRUE)
  }
  
  if ("DoubleSigmoid_Duration" %in% names(metrics) && any(is.finite(metrics$DoubleSigmoid_Duration))) {
    summary_df$Mean_DoubleSigmoid_Duration <- mean(metrics$DoubleSigmoid_Duration, na.rm = TRUE)
  }
  
  if ("Sigmoid_FFT_Fraction" %in% names(metrics) && any(is.finite(metrics$Sigmoid_FFT_Fraction))) {
    summary_df$Sigmoid_FFT_Fraction <- mean(metrics$Sigmoid_FFT_Fraction, na.rm = TRUE)
  }
  
  summary_df
}


# ------------------------------------------------------------
# Baseline sensitivity analysis for Wavelet module
# ------------------------------------------------------------
# This function keeps the same wavelet-detected candidate events and recomputes
# post-detection metrics under several baseline definitions. The goal is not to
# change ridge detection, but to evaluate how amplitude, AUC, FWHM, prominence,
# and related metrics change when the baseline definition changes.
build_wavelet_baseline_sensitivity_table <- function(time,
                                                     signal,
                                                     events,
                                                     event_window = 2,
                                                     initial_fraction = 0.10,
                                                     minimal_fraction = 0.10,
                                                     rolling_window = 20,
                                                     rolling_percentile = 2,
                                                     event_specific_fft_f = 0.20,
                                                     local_pre_event_window = 5,
                                                     local_f0_statistic = "median",
                                                     min_amp_above_baseline = 0,
                                                     accepted_min_ridge_length = 7,
                                                     min_cwt_amplitude = 0,
                                                     summary_event_set = "accepted") {
  baseline_grid <- data.frame(
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
    Baseline_Code = c(
      "f0_initial",
      "f0_minimal",
      "f0_smooth",
      "f0_event_specific",
      "f0_onset_mean",
      "median",
      "mean",
      "minimum"
    ),
    stringsAsFactors = FALSE
  )
  
  time_onset_for_baseline_grid <- estimate_wavelet_time_onset_from_events(
    time = time,
    signal = signal,
    events = events,
    event_window = event_window
  )
  
  out <- lapply(seq_len(nrow(baseline_grid)), function(i) {
    method_code <- baseline_grid$Baseline_Code[i]
    
    baseline_info <- compute_wavelet_baseline_info(
      time = time,
      signal = signal,
      baseline_method = method_code,
      initial_fraction = initial_fraction,
      minimal_fraction = minimal_fraction,
      rolling_window = rolling_window,
      rolling_percentile = rolling_percentile,
      time_onset = if (identical(method_code, "f0_onset_mean")) time_onset_for_baseline_grid else NA_real_
    )
    
    metrics_all <- compute_wavelet_event_metrics(
      time = time,
      signal = signal,
      events = events,
      event_window = event_window,
      baseline_method = method_code,
      initial_fraction = initial_fraction,
      minimal_fraction = minimal_fraction,
      rolling_window = rolling_window,
      rolling_percentile = rolling_percentile,
      event_specific_fft_f = event_specific_fft_f,
      local_pre_event_window = local_pre_event_window,
      local_f0_statistic = local_f0_statistic,
      min_amp_above_baseline = min_amp_above_baseline,
      accepted_min_ridge_length = accepted_min_ridge_length,
      min_cwt_amplitude = min_cwt_amplitude
    )
    
    metrics <- metrics_all
    if (!is.null(metrics) && nrow(metrics) > 0) {
      if (summary_event_set == "accepted") {
        metrics <- metrics[metrics$Event_Status == "Accepted", , drop = FALSE]
      }
      
      if (summary_event_set == "accepted_low") {
        metrics <- metrics[
          metrics$Event_Status == "Accepted" |
            grepl("^Low confidence", metrics$Event_Status),
          ,
          drop = FALSE
        ]
      }
    }
    
    baseline_values <- NULL
    
    if (!is.null(baseline_info$type) && baseline_info$type == "dynamic" &&
        !is.null(baseline_info$trace) && nrow(baseline_info$trace) > 0) {
      baseline_values <- suppressWarnings(as.numeric(baseline_info$trace$Baseline))
    } else if (!is.null(baseline_info$type) && baseline_info$type == "event_specific") {
      if (!is.null(metrics_all) && nrow(metrics_all) > 0 && "Baseline" %in% names(metrics_all)) {
        baseline_values <- suppressWarnings(as.numeric(metrics_all$Baseline))
      } else {
        baseline_values <- suppressWarnings(as.numeric(baseline_info$scalar))
      }
    } else {
      baseline_values <- suppressWarnings(as.numeric(baseline_info$scalar))
    }
    
    if (length(baseline_values) == 0 || all(!is.finite(baseline_values))) {
      baseline_values <- NA_real_
    }
    
    summary_df <- summarize_wavelet_metrics(
      metrics = metrics,
      time = time,
      signal = signal,
      baseline_info = baseline_info
    )
    
    auc_df <- build_auc_summary_table(
      metrics = metrics,
      time = time,
      signal = signal,
      baseline_info = baseline_info
    )
    
    data.frame(
      Baseline_Method = baseline_grid$Baseline_Method[i],
      Baseline_Type = ifelse(
        baseline_info$type == "dynamic",
        "Dynamic / time-varying",
        ifelse(baseline_info$type == "event_specific", "Event-specific", "Static / scalar")
      ),
      Baseline_Mean = mean(baseline_values, na.rm = TRUE),
      Baseline_Median = stats::median(baseline_values, na.rm = TRUE),
      Baseline_Min = min(baseline_values, na.rm = TRUE),
      Baseline_Max = max(baseline_values, na.rm = TRUE),
      Events_Used = ifelse(is.null(metrics), 0, nrow(metrics)),
      Number_of_Wavelet_Events = summary_df$Number_of_Wavelet_Events[1],
      Frequency = summary_df$Frequency[1],
      Mean_Amplitude_Above_Baseline = summary_df$Mean_Amplitude_Above_Baseline[1],
      Mean_DeltaF_over_Baseline_Peak = summary_df$Mean_DeltaF_over_F0_Peak[1],
      Mean_Prominence = summary_df$Mean_Prominence[1],
      Mean_FWHM = summary_df$Mean_FWHM[1],
      Mean_FWHP = summary_df$Mean_FWHP[1],
      Mean_Event_AUC = auc_df$Mean_Event_AUC[1],
      Total_Event_AUC = auc_df$Total_Event_AUC[1],
      Total_AUC_Global_Above_Baseline = auc_df$Total_AUC_Global_Above_Baseline[1],
      Mean_CWT_Amplitude = summary_df$Mean_CWT_Amplitude[1],
      Mean_Ridge_Length = summary_df$Mean_Ridge_Length[1],
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, out)
}

# ------------------------------------------------------------
# UI
# ------------------------------------------------------------
mod_wavelet_ridgewalking_ui <- function(id) {
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
        #   ns("wavelet_info_button"), "Help",
        #   class = "btn-sm",
        #   style = "position:absolute; top:0; right:15px; margin:5px;"
        # ),
        
        radioButtons(
          ns("data_simulate"), "1. Example Data",
          choices = c("Yes" = 1, "No" = 0),
          selected = 1
        ),
        
        fileInput(
          ns("file"),
          accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".tsv"),
          label = h5("2. Dataset")
        ),
        
        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),
        
        numericInput(
          ns("roi"),
          "3. Region of Interest (ROI):",
          value = 1,
          min = 1
        ),
        
        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),
        
        tags$h4("4. Wavelet Ridgewalking Function Arguments", style = "color: gray; margin-top: 10px;"),
        
        fluidRow(
          column(
            width = 6,
            sliderInput(
              ns("scale_min"),
              "4.1. Minimum Scale",
              min = 1,
              max = 30,
              value = 2,
              step = 1
            ),
            numericInput(
              ns("n_scales"),
              "4.2. Number of Scales",
              value = 32,
              min = 8,
              max = 100,
              step = 1
            ),
            numericInput(
              ns("time_tolerance"),
              "4.3. Ridge Connection Tolerance",
              value = 2,
              min = 0.1,
              step = 0.5
            ),
            numericInput(
              ns("min_event_distance"),
              "4.4. Min Event Distance",
              value = 1,
              min = 0,
              step = 0.5
            )
          ),
          column(
            width = 6,
            sliderInput(
              ns("scale_max"),
              "4.5. Maximum Scale",
              min = 5,
              max = 100,
              value = 40,
              step = 1
            ),
            sliderInput(
              ns("coef_quantile"),
              "4.6. CWT Amplitude Quantile",
              min = 0.50,
              max = 0.99,
              value = 0.90,
              step = 0.01
            ),
            numericInput(
              ns("min_ridge_length"),
              "4.7. Min Ridge Length",
              value = 5,
              min = 1,
              step = 1
            ),
            numericInput(
              ns("event_window"),
              "4.8. Metric Window Around Event",
              value = 2,
              min = 0.1,
              step = 0.5
            )
          )
        ),
        
        helpText(
          "These arguments control the multiscale ridge detection step. They do not change the core wavelet ridgewalking logic."
        ),
        
        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),
        
        selectInput(
          ns("baseline_method"),
          "5. Baseline:",
          choices = c(
            "5.1. Standard definition" = "f0_initial",
            "5.2. Low-fluorescence region baseline" = "f0_minimal",
            "5.3. Rolling percentile baseline" = "f0_smooth",
            "5.4. Event-specific baseline: left-right FFT minima" = "f0_event_specific",
            "5.5. Constant baseline: Mean from trace start to Time_Onset" = "f0_onset_mean",
            "5.6. Median signal (legacy)" = "median",
            "5.7. Mean signal (legacy)" = "mean",
            "5.8. Min" = "minimum"
          ),
          selected = "f0_initial"
        ),
        
        conditionalPanel(
          condition = "input.baseline_method == 'f0_initial'",
          ns = ns,
          numericInput(
            ns("initial_f0_fraction"),
            "5.1.1. Initial Trace Fraction:",
            value = 0.10,
            min = 0.01,
            max = 1,
            step = 0.01
          ),
          tags$div(
            class = "helper-card",
            tags$strong("Standard definition"),
            tags$p(
              "Uses an initial fraction of the trace to calculate a global baseline for wavelet-derived metrics."
            )
          )
        ),
        
        conditionalPanel(
          condition = "input.baseline_method == 'f0_minimal'",
          ns = ns,
          numericInput(
            ns("minimal_f0_fraction"),
            "5.2.1. Minimal Lower Fraction:",
            value = 0.10,
            min = 0.01,
            max = 1,
            step = 0.01
          )
        ),
        
        conditionalPanel(
          condition = "input.baseline_method == 'f0_smooth'",
          ns = ns,
          fluidRow(
            column(
              6,
              numericInput(
                ns("rolling_window"),
                "5.3.1. Moving Window Size:",
                value = 20,
                min = 0.0001,
                step = 1
              )
            ),
            column(
              6,
              numericInput(
                ns("rolling_percentile"),
                "5.3.2. Percentile:",
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
          condition = "input.baseline_method == 'f0_event_specific'",
          ns = ns,
          tags$div(
            class = "helper-card",
            tags$strong("Event-specific baseline"),
            tags$p(
              "For each transient, the baseline is drawn between the FFT-smoothed minimum to the left of the peak ",
              "and the FFT-smoothed minimum to the right of the peak."
            )
          ),
          sliderInput(
            ns("event_specific_fft_f"),
            "5.4.1. FFT Smoothing Fraction for Event-Specific Baseline:",
            min = 0.05,
            max = 0.50,
            value = 0.20,
            step = 0.01
          ),
          numericInput(
            ns("local_pre_event_window"),
            "5.4.2. Fallback Pre-Event Window:",
            value = 5,
            min = 0.1,
            step = 0.5
          ),
          selectInput(
            ns("local_f0_statistic"),
            "5.4.3. Fallback Baseline Statistic:",
            choices = c(
              "Median of pre-event window" = "median",
              "10th percentile of pre-event window" = "p10"
            ),
            selected = "median"
          )
        ),
        
        conditionalPanel(
          condition = "input.baseline_method == 'f0_onset_mean'",
          ns = ns,
          tags$div(
            class = "helper-card",
            tags$strong("Constant baseline: mean from trace start to Time_Onset"),
            tags$p(
              "This option finds the first Time_Onset from the wavelet-detected transients, takes the calcium trace from the beginning up to that Time_Onset, ",
              "calculates the mean of that portion, and uses that single constant value as the baseline for all transients."
            )
          )
        ),
        
        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),
        
        tags$h4("6. Biological Event Filters", style = "color: gray; margin-top: 10px;"),
        
        fluidRow(
          column(
            width = 6,
            numericInput(
              ns("min_amp_above_baseline"),
              "6.1. Amplitude Above Baseline (min)",
              value = 0,
              min = 0,
              step = 0.05
            ),
            numericInput(
              ns("min_cwt_amplitude"),
              "6.2. CWT Amplitude (min)",
              value = 0,
              min = 0,
              step = 0.1
            )
          ),
          column(
            width = 6,
            numericInput(
              ns("accepted_min_ridge_length"),
              "6.3. Ridge Length for Accepted Status",
              value = 7,
              min = 1,
              step = 1
            ),
            selectInput(
              ns("summary_event_set"),
              "6.4. Events Used in Summary Metrics:",
              choices = c(
                "Accepted events only" = "accepted",
                "Accepted + low-confidence events" = "accepted_low",
                "All detected candidates" = "all"
              ),
              selected = "accepted"
            )
          )
        ),
        
        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),
        
        tags$h4("7. Graph Display Options", style = "color: gray; margin-top: 10px;"),
        
        radioButtons(
          ns("show_fwhp"),
          "7.1. Full Width at Half Prominence (FWHP):",
          choices = c("No" = 0, "Yes" = 1),
          selected = 0
        ),
        
        radioButtons(
          ns("show_nonaccepted_event_lines"),
          "7.2. Show non-accepted event lines:",
          choices = c("No" = 0, "Yes" = 1),
          selected = 0
        ),
        
        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),
        
        tags$h4("8. Optional Sigmoid Modeling", style = "color: gray; margin-top: 10px;"),
        tags$div(
          class = "helper-card",
          tags$strong("Sigmoid fits after Fourier smoothing"),
          tags$p(
            "Wavelet ridgewalking continues to detect the candidate transients. ",
            "When these options are activated, the selected ROI trace is first smoothed with an FFT low-pass reconstruction, ",
            "and the single or double sigmoid is fitted on that FFT-smoothed trace. ",
            "FWHP and FWHM remain fixed from the original non-sigmoid wavelet metrics."
          )
        ),
        checkboxInput(
          ns("use_single_sigmoid_wavelet"),
          "8.1. Apply Single Sigmoid Fit for Onset and Rise-Time Metrics",
          value = FALSE
        ),
        checkboxInput(
          ns("use_double_sigmoid_wavelet"),
          "8.2. Apply Double Sigmoid Fit for Full-Transient Metrics",
          value = FALSE
        ),
        conditionalPanel(
          condition = "input.use_single_sigmoid_wavelet == true || input.use_double_sigmoid_wavelet == true",
          ns = ns,
          sliderInput(
            ns("sigmoid_fft_fraction_wavelet"),
            "8.3. FFT Low-Frequency Fraction for Sigmoid Fits:",
            min = 0.01,
            max = 1,
            value = 0.20,
            step = 0.01
          ),
          checkboxInput(
            ns("sigmoid_keep_mean_wavelet"),
            "8.4. Keep Mean (DC Component) for Sigmoid FFT Smoothing",
            value = TRUE
          ),
          numericInput(
            ns("sigmoid_pre_event_window_wavelet"),
            "8.5. Pre-Event Window for Baseline RMS:",
            value = 5,
            min = 0.1,
            step = 0.5
          )
        ),
        conditionalPanel(
          condition = "input.use_single_sigmoid_wavelet == true",
          ns = ns,
          fluidRow(
            column(
              6,
              numericInput(
                ns("single_sigmoid_left_pad_wavelet"),
                "8.6. Single Sigmoid Left Padding:",
                value = 0,
                min = 0,
                step = 0.25
              )
            ),
            column(
              6,
              numericInput(
                ns("single_sigmoid_right_pad_wavelet"),
                "8.7. Single Sigmoid Right Padding:",
                value = 0,
                min = 0,
                step = 0.25
              )
            )
          )
        ),
        conditionalPanel(
          condition = "input.use_double_sigmoid_wavelet == true",
          ns = ns,
          fluidRow(
            column(
              6,
              numericInput(
                ns("double_sigmoid_left_pad_wavelet"),
                "8.8. Double Sigmoid Left Padding:",
                value = 0,
                min = 0,
                step = 0.25
              )
            ),
            column(
              6,
              numericInput(
                ns("double_sigmoid_right_pad_wavelet"),
                "8.9. Double Sigmoid Right Padding:",
                value = 0,
                min = 0,
                step = 0.25
              )
            )
          )
        ),
        
        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),
        
        tags$h4("9. Download Options", style = "color: gray; margin-top: 10px;"),
        downloadButton(ns("descargarP"), "9.1. Trace Metrics"),
        downloadButton(ns("descargar"), "9.2. Transient Metrics"),
        downloadButton(ns("Calcium_Trance_Graph"), "9.3. Calcium Trace Graph")
      ),
      
      mainPanel(
        tabsetPanel(
          tabPanel(
            "10. SummaryData",
            tags$h4("10.1. Dataset Summary", style = "color: gray; margin-top: 10px;"),
            DT::DTOutput(ns("infotable2")),
            tags$br(),
            tags$h4("10.2. Dataset", style = "color: gray; margin-top: 10px;"),
            DT::DTOutput(ns("data2")),
            tags$br(),
            tags$h4("10.3. Selected ROI Trace", style = "color: gray; margin-top: 10px;"),
            plotOutput(ns("selected_trace_plot"), height = "420px"),
            tags$br(),
            tags$h4("10.4. Selected ROI Descriptive Summary", style = "color: gray; margin-top: 10px;"),
            DT::DTOutput(ns("selected_trace_summary"))
          ),
          
          tabPanel(
            "11. Peaks",
            plotOutput(ns("event_plot"), height = "520px"),
            tags$br(),
            tags$h4("11.1. AUC Summary", style = "color: gray; margin-top: 10px;"),
            DT::DTOutput(ns("auc_summary_table")),
            tags$br(),
            tags$h4("11.2. Event-Level Metrics", style = "color: gray; margin-top: 10px;"),
            DT::DTOutput(ns("event_table"))
          ),
          
          tabPanel(
            "12. Baseline Sensitivity Analysis",
            div(
              class = "helper-card",
              tags$h4("12.1. Methodological Note: Why Baseline Matters"),
              tags$p(
                "Wavelet ridgewalking detects candidate transients before the baseline is applied. ",
                "However, post-detection metrics such as amplitude above baseline, DeltaF/Baseline, AUC, FWHM, FWHP, and prominence ",
                "can change depending on the selected baseline definition."
              ),
              tags$p(
                tags$strong("Purpose: "),
                "this panel recomputes the wavelet-derived metrics using the same detected ridge candidates under several baseline definitions. ",
                "Use it to evaluate whether the biological interpretation is robust or baseline-sensitive."
              ),
              tags$p(
                tags$strong("Interpretation: "),
                "if the values are similar across baseline methods, the analysis is relatively robust. ",
                "If the values change strongly, the report should explicitly state that the results are baseline-sensitive."
              )
            ),
            tags$br(),
            DT::DTOutput(ns("baseline_sensitivity_table")),
            tags$br(),
            plotOutput(ns("baseline_sensitivity_plot"), height = "900px"),
            tags$br(),
            uiOutput(ns("baseline_sensitivity_interpretation"))
          ),
          
          tabPanel(
            "13. Metrics",
            tags$h4("13.1. Transient Metrics", style = "color: gray; margin-top: 10px;"),
            DT::DTOutput(ns("transient_metrics_table")),
            tags$br(),
            tags$h4("13.2. Trace Metrics", style = "color: gray; margin-top: 10px;"),
            DT::DTOutput(ns("trace_metrics_table"))
          ),
          
          tabPanel(
            "14. Metric Plots",
            plotOutput(ns("summary_plot"), height = "950px")
          ),
          
          tabPanel(
            "15. Wavelet Analysis",
            plotOutput(ns("cwt_plot"), height = "650px")
          )
        )
      )
    )
  )
}

# ------------------------------------------------------------
# Server
# ------------------------------------------------------------

mod_wavelet_ridgewalking_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    
    observeEvent(input$wavelet_info_button, {
      showModal(modalDialog(
        title = "Wavelet Ridgewalking module",
        tags$p(
          "This module follows the same visual organization used in the FFT + Baseline Analysis module: ",
          "data input, ROI selection, method arguments, baseline options, biological filters, and graph display options."
        ),
        tags$p(
          "The wavelet ridgewalking detection logic was preserved. The updates focus on labels, tab names, section names, helper text, and output organization."
        ),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
    
    filedata <- reactive({
      if (as.numeric(input$data_simulate) > 0) {
        return(example_calcium_wavelet_data())
      }
      
      req(input$file)
      ext <- tools::file_ext(input$file$name)
      
      if (exists("load_file")) {
        dat <- load_file(input$file$name, input$file$datapath, ext)
      } else {
        dat <- switch(
          ext,
          csv = vroom::vroom(input$file$datapath, delim = ","),
          tsv = vroom::vroom(input$file$datapath, delim = "\t"),
          json = jsonlite::fromJSON(input$file$datapath),
          stop("Unsupported file format.")
        )
      }
      
      as.data.frame(dat)
    })
    
    selected_trace <- reactive({
      dat <- filedata()
      shiny::validate(shiny::need(ncol(dat) >= 2, "Dataset must contain Time plus at least one ROI."))
      
      roi_col <- as.numeric(input$roi) + 1
      shiny::validate(shiny::need(roi_col <= ncol(dat), "Selected ROI is out of range."))
      
      data.frame(
        Time = as.numeric(dat[[1]]),
        signal = as.numeric(dat[[roi_col]])
      )
    })
    
    data_info <- reactive({
      dat <- filedata()
      shiny::validate(shiny::need(ncol(dat) >= 2, "Dataset must contain Time plus at least one ROI."))
      n_observations <- nrow(dat)
      n_rois <- ncol(dat) - 1
      summary_data <- data.frame(
        Number = c(n_rois, n_observations),
        row.names = c("Region of Interest (ROI)", "Time observations")
      )
      list(
        SummaryData = summary_data,
        data = data.frame(dat, row.names = NULL)
      )
    })
    
    output$data2 <- DT::renderDT({
      DT::datatable(
        data_info()$data,
        options = list(pagingType = "simple"),
        caption = htmltools::tags$caption(htmltools::tags$strong("Dataset:"))
      )
    })
    
    output$infotable2 <- DT::renderDT({
      DT::datatable(
        data_info()$SummaryData,
        options = list(pagingType = "simple", dom = "t"),
        caption = htmltools::tags$caption(htmltools::tags$strong("Dataset Summary:"))
      )
    })
    
    output$selected_trace_plot <- renderPlot({
      df <- selected_trace()
      ggplot(df, aes(x = Time, y = signal)) +
        geom_line(linewidth = 0.6) +
        labs(
          title = "Selected ROI Calcium Trace",
          subtitle = paste0("ROI: ", input$roi),
          x = "Time",
          y = "Signal"
        ) +
        theme_bw(base_size = 13) +
        theme(
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank()
        )
    })
    
    output$selected_trace_summary <- DT::renderDT({
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
          as.character(input$roi),
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
      
      DT::datatable(
        trace_summary,
        rownames = FALSE,
        options = list(pagingType = "simple", dom = "t"),
        caption = htmltools::tags$caption(htmltools::tags$strong("Selected ROI Summary:"))
      )
    })
    
    wavelet_results <- reactive({
      df <- selected_trace()
      
      scale_min <- min(input$scale_min, input$scale_max)
      scale_max <- max(input$scale_min, input$scale_max)
      
      n_scales <- safe_int(input$n_scales, 32L)
      if (!is.finite(n_scales) || n_scales < 8) n_scales <- 32L
      
      scales <- seq(scale_min, scale_max, length.out = n_scales)
      
      cwt <- compute_cwt_ricker(
        time = df$Time,
        signal = df$signal,
        scales = scales
      )
      
      maxima <- find_local_maxima_by_scale(
        amp = cwt$amplitude,
        time = df$Time,
        scales = scales,
        min_coef_quantile = input$coef_quantile
      )
      
      ridges <- build_wavelet_ridges(
        maxima_df = maxima,
        time_tolerance = input$time_tolerance,
        min_ridge_length = input$min_ridge_length
      )
      
      events <- merge_close_wavelet_events(
        events = ridges,
        min_distance_time = input$min_event_distance
      )
      
      metrics <- compute_wavelet_event_metrics(
        time = df$Time,
        signal = df$signal,
        events = events,
        event_window = input$event_window,
        baseline_method = input$baseline_method,
        initial_fraction = input$initial_f0_fraction,
        minimal_fraction = input$minimal_f0_fraction,
        rolling_window = input$rolling_window,
        rolling_percentile = input$rolling_percentile,
        event_specific_fft_f = input$event_specific_fft_f,
        local_pre_event_window = input$local_pre_event_window,
        local_f0_statistic = input$local_f0_statistic,
        min_amp_above_baseline = input$min_amp_above_baseline,
        accepted_min_ridge_length = input$accepted_min_ridge_length,
        min_cwt_amplitude = input$min_cwt_amplitude
      )
      
      sigmoid_fft_fraction <- safe_num(input$sigmoid_fft_fraction_wavelet, 0.20)
      sigmoid_fft_fraction <- max(0.01, min(1, sigmoid_fft_fraction))
      sigmoid_keep_mean <- isTRUE(input$sigmoid_keep_mean_wavelet)
      sigmoid_fft_res <- wavelet_fft_lowpass_smooth_trace(
        time = df$Time,
        signal = df$signal,
        f = sigmoid_fft_fraction,
        keep_mean = sigmoid_keep_mean
      )
      sigmoid_smoothed_trace <- sigmoid_fft_res$data
      sigmoid_curves <- data.frame()
      sigmoids_active <- isTRUE(input$use_single_sigmoid_wavelet) ||
        isTRUE(input$use_double_sigmoid_wavelet)
      
      if (sigmoids_active) {
        sigmoid_out <- append_wavelet_sigmoid_metrics(
          metrics = metrics,
          time = sigmoid_smoothed_trace$Time,
          signal = sigmoid_smoothed_trace$signal,
          use_single_sigmoid = isTRUE(input$use_single_sigmoid_wavelet),
          use_double_sigmoid = isTRUE(input$use_double_sigmoid_wavelet),
          pre_event_window = input$sigmoid_pre_event_window_wavelet,
          single_left_pad = input$single_sigmoid_left_pad_wavelet,
          single_right_pad = input$single_sigmoid_right_pad_wavelet,
          double_left_pad = input$double_sigmoid_left_pad_wavelet,
          double_right_pad = input$double_sigmoid_right_pad_wavelet,
          fft_fraction = sigmoid_fft_fraction
        )
        metrics <- sigmoid_out$metrics
        sigmoid_curves <- sigmoid_out$fitted_curves
      }
      
      # Keep the static trace-start-to-onset baseline internally consistent.
      # This is especially important after single-sigmoid onset refinement:
      # the black dashed selected-baseline line, blue amplitude segment,
      # Event_AUC, and exported baseline columns must all use the same value.
      if (identical(as.character(input$baseline_method), "f0_onset_mean")) {
        metrics <- sync_wavelet_onset_mean_baseline_with_metrics(
          metrics = metrics,
          time = df$Time,
          signal = df$signal
        )
      }
      
      if (sigmoids_active || identical(as.character(input$baseline_method), "f0_onset_mean")) {
        metrics <- update_wavelet_event_auc_from_metrics(
          metrics = metrics,
          time = df$Time,
          signal = df$signal
        )
      }
      
      list(
        trace = df,
        cwt = cwt,
        maxima = maxima,
        ridges = ridges,
        events = events,
        metrics = metrics,
        sigmoid_fft_res = sigmoid_fft_res,
        sigmoid_smoothed_trace = sigmoid_smoothed_trace,
        sigmoid_curves = sigmoid_curves,
        sigmoid_fft_fraction = sigmoid_fft_fraction,
        sigmoids_active = sigmoids_active
      )
    })
    
    filtered_metrics_for_summary <- reactive({
      metrics <- wavelet_results()$metrics
      
      if (is.null(metrics) || nrow(metrics) == 0) return(metrics)
      
      if (input$summary_event_set == "accepted") {
        metrics <- metrics[metrics$Event_Status == "Accepted", , drop = FALSE]
      }
      
      if (input$summary_event_set == "accepted_low") {
        metrics <- metrics[
          metrics$Event_Status == "Accepted" |
            grepl("^Low confidence", metrics$Event_Status),
          ,
          drop = FALSE
        ]
      }
      
      metrics
    })
    
    selected_baseline_info <- reactive({
      df <- selected_trace()
      
      time_onset_for_baseline <- NA_real_
      if (identical(as.character(input$baseline_method), "f0_onset_mean")) {
        metrics_for_onset <- wavelet_results()$metrics
        
        if (!is.null(metrics_for_onset) && nrow(metrics_for_onset) > 0) {
          # Prefer the synchronized Baseline column when available. This ensures
          # the black dashed selected-baseline line is exactly the same value as
          # the baseline used by the blue amplitude segments and the exports.
          if ("Baseline" %in% names(metrics_for_onset)) {
            baseline_candidates <- suppressWarnings(as.numeric(metrics_for_onset$Baseline))
            baseline_candidates <- baseline_candidates[is.finite(baseline_candidates)]
            if (length(baseline_candidates) > 0) {
              onset_candidates <- numeric(0)
              if ("Transient_Occurrence_Time" %in% names(metrics_for_onset)) {
                onset_candidates <- suppressWarnings(as.numeric(metrics_for_onset$Transient_Occurrence_Time))
              } else if ("Onset_Time" %in% names(metrics_for_onset)) {
                onset_candidates <- suppressWarnings(as.numeric(metrics_for_onset$Onset_Time))
              }
              onset_candidates <- onset_candidates[is.finite(onset_candidates)]
              time_onset_for_baseline <- if (length(onset_candidates) > 0) min(onset_candidates, na.rm = TRUE) else NA_real_
              
              return(list(
                method = "Constant baseline: mean from trace start to Time_Onset",
                type = "static",
                scalar = baseline_candidates[1],
                trace = NULL,
                time_onset = time_onset_for_baseline,
                source = "Synchronized metrics Baseline column"
              ))
            }
          }
          
          if ("Transient_Occurrence_Time" %in% names(metrics_for_onset)) {
            onset_candidates <- suppressWarnings(as.numeric(metrics_for_onset$Transient_Occurrence_Time))
            onset_candidates <- onset_candidates[is.finite(onset_candidates)]
            if (length(onset_candidates) > 0) {
              time_onset_for_baseline <- min(onset_candidates, na.rm = TRUE)
            }
          }
        }
      }
      
      compute_wavelet_baseline_info(
        time = df$Time,
        signal = df$signal,
        baseline_method = input$baseline_method,
        initial_fraction = input$initial_f0_fraction,
        minimal_fraction = input$minimal_f0_fraction,
        rolling_window = input$rolling_window,
        rolling_percentile = input$rolling_percentile,
        time_onset = time_onset_for_baseline
      )
    })
    
    selected_baseline_value <- reactive({
      baseline_info <- selected_baseline_info()
      metrics <- wavelet_results()$metrics
      df <- selected_trace()
      
      if (!is.null(baseline_info$type) && baseline_info$type == "event_specific") {
        if (!is.null(metrics) && nrow(metrics) > 0 && "Baseline" %in% names(metrics)) {
          return(stats::median(metrics$Baseline, na.rm = TRUE))
        }
        return(stats::median(df$signal, na.rm = TRUE))
      }
      
      baseline_info$scalar
    })
    
    output$cwt_plot <- renderPlot({
      x <- wavelet_results()
      
      amp_df <- expand.grid(
        Scale = x$cwt$scales,
        Time = x$cwt$time
      )
      amp_df$Amplitude <- as.vector(x$cwt$amplitude)
      
      p <- ggplot2::ggplot(amp_df, ggplot2::aes(x = Time, y = Scale, fill = Amplitude)) +
        ggplot2::geom_raster(interpolate = TRUE) +
        ggplot2::labs(
          title = "Wavelet Analysis: Continuous Wavelet Transform Amplitude Map",
          subtitle = "Local maxima indicate ridge candidates used for event detection",
          x = "Time",
          y = "Scale",
          fill = "CWT amplitude"
        ) +
        ggplot2::theme_minimal(base_size = 13)
      
      if (!is.null(x$maxima) && nrow(x$maxima) > 0) {
        p <- p + ggplot2::geom_point(
          data = x$maxima,
          ggplot2::aes(x = Time, y = scale),
          inherit.aes = FALSE,
          size = 0.7
        )
      }
      
      p
    })
    
    wavelet_event_plot <- reactive({
      x <- wavelet_results()
      df <- x$trace
      all_metrics <- x$metrics
      metrics <- filtered_metrics_for_summary()
      baseline_info <- selected_baseline_info()
      baseline <- selected_baseline_value()
      
      global_auc <- compute_global_auc_above_baseline(df$Time, df$signal, baseline_info)
      event_auc_total <- if (!is.null(metrics) && nrow(metrics) > 0) {
        sum(metrics$Event_AUC, na.rm = TRUE)
      } else {
        0
      }
      
      # Green area:
      # - For global/dynamic baseline: total AUC above the selected baseline.
      # - For event-specific baseline: local left-right FFT-minima interval for each transient.
      if (!is.null(baseline_info$type) && baseline_info$type == "event_specific") {
        global_auc_df <- build_event_specific_f0_interval_plot_data(
          time = df$Time,
          signal = df$signal,
          metrics = metrics
        )
      } else {
        global_auc_df <- build_global_auc_plot_data(
          time = df$Time,
          signal = df$signal,
          baseline_info = baseline_info
        )
      }
      
      event_auc_df <- build_event_auc_plot_data(
        time = df$Time,
        signal = df$signal,
        metrics = metrics
      )
      
      baseline_line_df <- build_baseline_line_plot_data(
        time = df$Time,
        baseline_info = baseline_info,
        metrics = metrics
      )
      
      p <- ggplot2::ggplot(df, ggplot2::aes(x = Time, y = signal))
      
      # Total AUC between baseline and curve across the full trace.
      # Total AUC / green area
      if (!is.null(global_auc_df) && nrow(global_auc_df) > 0) {
        
        if (!is.null(baseline_info$type) && baseline_info$type == "event_specific") {
          
          p <- p +
            ggplot2::geom_ribbon(
              data = global_auc_df,
              ggplot2::aes(
                x = Time,
                ymin = Baseline,
                ymax = AUC_Upper,
                group = Event_ID
              ),
              inherit.aes = FALSE,
              fill = "#E1FFE1",
              alpha = 0.85
            )
          
        } else {
          
          p <- p +
            ggplot2::geom_ribbon(
              data = global_auc_df,
              ggplot2::aes(
                x = Time,
                ymin = Baseline,
                ymax = AUC_Upper
              ),
              inherit.aes = FALSE,
              fill = "#E1FFE1",
              alpha = 0.85
            )
        }
      }
      
      # Event-specific AUC between event baseline and curve.
      if (!is.null(event_auc_df) && nrow(event_auc_df) > 0) {
        p <- p +
          ggplot2::geom_ribbon(
            data = event_auc_df,
            ggplot2::aes(
              x = Time,
              ymin = Baseline,
              ymax = AUC_Upper,
              group = Event_ID
            ),
            inherit.aes = FALSE,
            fill = "#FFD166",
            alpha = 0.60
          )
      }
      
      p <- p +
        ggplot2::geom_line(linewidth = 0.75)
      
      # The FFT-smoothed trace is used internally as the input for the
      # optional single/double sigmoid fits, but it is intentionally not drawn
      # on the Peaks plot. This keeps the plot focused on the original calcium
      # trace and the fitted sigmoid model curves.
      
      # Optional single/double sigmoid fitted curves.
      if (!is.null(x$sigmoid_curves) && nrow(x$sigmoid_curves) > 0) {
        sigmoid_curves_plot <- x$sigmoid_curves
        if (!is.null(metrics) && nrow(metrics) > 0 && "Event_ID" %in% names(metrics)) {
          sigmoid_curves_plot <- sigmoid_curves_plot[
            sigmoid_curves_plot$Event_ID %in% metrics$Event_ID,
            ,
            drop = FALSE
          ]
        }
        
        single_curves <- sigmoid_curves_plot[
          sigmoid_curves_plot$Model == "Single sigmoid",
          ,
          drop = FALSE
        ]
        double_curves <- sigmoid_curves_plot[
          sigmoid_curves_plot$Model == "Double sigmoid",
          ,
          drop = FALSE
        ]
        
        if (nrow(single_curves) > 0) {
          p <- p +
            ggplot2::geom_line(
              data = single_curves,
              ggplot2::aes(x = Time, y = Fitted, group = Event_ID),
              inherit.aes = FALSE,
              linewidth = 1.05,
              color = "#0072B2"
            )
        }
        
        if (nrow(double_curves) > 0) {
          p <- p +
            ggplot2::geom_line(
              data = double_curves,
              ggplot2::aes(x = Time, y = Fitted, group = Event_ID),
              inherit.aes = FALSE,
              linewidth = 1.05,
              color = "#7B3294"
            )
        }
      }
      
      # Draw the selected baseline:
      # - static baseline: one horizontal dashed line,
      # - smooth baseline: one time-varying dashed curve,
      # - event-specific baseline: one dashed segment per transient.
      if (!is.null(baseline_line_df) && nrow(baseline_line_df) > 0) {
        p <- p +
          ggplot2::geom_line(
            data = baseline_line_df,
            ggplot2::aes(
              x = Time,
              y = Baseline,
              group = Baseline_Group
            ),
            inherit.aes = FALSE,
            linetype = 2,
            linewidth = 0.75
          )
      }
      
      p <- p +
        ggplot2::labs(
          title = "Peaks: Wavelet Ridgewalking Detected Events",
          subtitle = paste0(
            "Yellow = Event_AUC | Light green = selected baseline area | Black dashed line = selected baseline | Blue vertical over red = amplitude | Red vertical = prominence | Blue horizontal = FWHM | Optional red horizontal = FWHP. ",
            "Total global AUC = ", round(global_auc, 4),
            " | Total event AUC = ", round(event_auc_total, 4)
          ),
          x = "Time",
          y = "Calcium Signal"
        ) +
        ggplot2::theme_classic(base_size = 14)
      
      # Show vertical event lines.
      # Accepted event lines are always shown.
      # Low-confidence and rejected event lines are optional.
      if (!is.null(all_metrics) && nrow(all_metrics) > 0) {
        
        accepted_lines <- all_metrics[
          all_metrics$Event_Status == "Accepted",
          ,
          drop = FALSE
        ]
        
        nonaccepted_lines <- all_metrics[
          all_metrics$Event_Status %in% c(
            "Low confidence: short ridge",
            "Rejected: below baseline"
          ),
          ,
          drop = FALSE
        ]
        
        if (nrow(accepted_lines) > 0) {
          p <- p +
            ggplot2::geom_vline(
              data = accepted_lines,
              ggplot2::aes(xintercept = Wavelet_Event_Time, linetype = Event_Status),
              linewidth = 0.5
            )
        }
        
        if (as.numeric(input$show_nonaccepted_event_lines) > 0 && nrow(nonaccepted_lines) > 0) {
          p <- p +
            ggplot2::geom_vline(
              data = nonaccepted_lines,
              ggplot2::aes(xintercept = Wavelet_Event_Time, linetype = Event_Status),
              linewidth = 0.5
            )
        }
      }
      
      # Mark selected events and draw amplitude/prominence annotations.
      if (!is.null(metrics) && nrow(metrics) > 0) {
        
        # Small horizontal offset so amplitude and prominence do not overlap exactly.
        x_range <- max(df$Time, na.rm = TRUE) - min(df$Time, na.rm = TRUE)
        x_offset <- ifelse(is.finite(x_range) && x_range > 0, 0.006 * x_range, 0.2)
        
        metrics$Prominence_X <- metrics$Peak_Time + x_offset
        
        p <- p +
          # FWHM: full width at half maximum (yellow horizontal segment).
          ggplot2::geom_segment(
            data = metrics[is.finite(metrics$FWHM_Left_Time) & is.finite(metrics$FWHM_Right_Time), , drop = FALSE],
            ggplot2::aes(
              x = FWHM_Left_Time,
              xend = FWHM_Right_Time,
              y = FWHM_Level,
              yend = FWHM_Level
            ),
            inherit.aes = FALSE,
            linewidth = 2.0,
            color = "blue"
          ) +
          # Prominence: peak amplitude minus the FFT-smoothed left valley before the peak.
          # This layer is drawn first.
          ggplot2::geom_segment(
            data = metrics,
            ggplot2::aes(
              x = Prominence_X,
              xend = Prominence_X,
              y = Prominence_Reference,
              yend = Peak_Amplitude
            ),
            inherit.aes = FALSE,
            linewidth = 1.15,
            color = "red",
            linetype = "dotdash"
          ) +
          # Amplitude above baseline: peak amplitude minus baseline.
          # This layer is drawn after prominence so the blue line appears on top.
          ggplot2::geom_segment(
            data = metrics,
            ggplot2::aes(
              x = Peak_Time,
              xend = Peak_Time,
              y = Baseline,
              yend = Peak_Amplitude
            ),
            inherit.aes = FALSE,
            linewidth = 1.25,
            color = "blue",
            linetype = "dashed"
          ) +
          
          # Peak points.
          ggplot2::geom_point(
            data = metrics,
            ggplot2::aes(x = Peak_Time, y = Peak_Amplitude, shape = Event_Status),
            size = 2.6,
            inherit.aes = FALSE
          ) +
          # Onset/baseline points.
          ggplot2::geom_point(
            data = metrics,
            ggplot2::aes(x = Onset_Time, y = Baseline),
            size = 1.8,
            inherit.aes = FALSE
          )
        
        # Optional sigmoid onset/offset annotations.
        if ("SingleSigmoid_Onset_Time" %in% names(metrics)) {
          single_onsets <- metrics[
            is.finite(metrics$SingleSigmoid_Onset_Time),
            ,
            drop = FALSE
          ]
          if (nrow(single_onsets) > 0) {
            p <- p +
              ggplot2::geom_vline(
                data = single_onsets,
                ggplot2::aes(xintercept = SingleSigmoid_Onset_Time),
                inherit.aes = FALSE,
                linetype = "dotted",
                linewidth = 0.85,
                color = "#D55E00"
              )
          }
        }
        
        if ("DoubleSigmoid_Offset_Time" %in% names(metrics)) {
          double_offsets <- metrics[
            is.finite(metrics$DoubleSigmoid_Offset_Time),
            ,
            drop = FALSE
          ]
          if (nrow(double_offsets) > 0) {
            p <- p +
              ggplot2::geom_vline(
                data = double_offsets,
                ggplot2::aes(xintercept = DoubleSigmoid_Offset_Time),
                inherit.aes = FALSE,
                linetype = "dotted",
                linewidth = 0.85,
                color = "#CC79A7"
              )
          }
        }
        
        # Optional FWHP display controlled by the user.
        if (as.numeric(input$show_fwhp) > 0) {
          p <- p +
            ggplot2::geom_segment(
              data = metrics[
                is.finite(metrics$FWHP_Left_Time) &
                  is.finite(metrics$FWHP_Right_Time),
                ,
                drop = FALSE
              ],
              ggplot2::aes(
                x = FWHP_Left_Time,
                xend = FWHP_Right_Time,
                y = FWHP_Level,
                yend = FWHP_Level
              ),
              inherit.aes = FALSE,
              linewidth = 1.8,
              color = "red"
            )
        }
      }
      
      p
    })
    
    output$event_plot <- renderPlot({
      wavelet_event_plot()
    })
    
    output$auc_summary_table <- DT::renderDT({
      metrics <- filtered_metrics_for_summary()
      df <- selected_trace()
      baseline_info <- selected_baseline_info()
      
      auc_df <- build_auc_summary_table(
        metrics = metrics,
        time = df$Time,
        signal = df$signal,
        baseline_info = baseline_info
      )
      
      numeric_cols <- sapply(auc_df, is.numeric)
      auc_df[numeric_cols] <- lapply(auc_df[numeric_cols], function(z) round(z, 4))
      
      DT::datatable(
        auc_df,
        options = list(
          dom = "t",
          scrollX = TRUE
        )
      )
    })
    
    output$baseline_sensitivity_table <- DT::renderDT({
      df <- selected_trace()
      results <- wavelet_results()
      
      sensitivity_df <- build_wavelet_baseline_sensitivity_table(
        time = df$Time,
        signal = df$signal,
        events = results$events,
        event_window = input$event_window,
        initial_fraction = input$initial_f0_fraction,
        minimal_fraction = input$minimal_f0_fraction,
        rolling_window = input$rolling_window,
        rolling_percentile = input$rolling_percentile,
        event_specific_fft_f = input$event_specific_fft_f,
        local_pre_event_window = input$local_pre_event_window,
        local_f0_statistic = input$local_f0_statistic,
        min_amp_above_baseline = input$min_amp_above_baseline,
        accepted_min_ridge_length = input$accepted_min_ridge_length,
        min_cwt_amplitude = input$min_cwt_amplitude,
        summary_event_set = input$summary_event_set
      )
      
      numeric_cols <- sapply(sensitivity_df, is.numeric)
      sensitivity_df[numeric_cols] <- lapply(
        sensitivity_df[numeric_cols],
        function(z) round(z, 4)
      )
      
      DT::datatable(
        sensitivity_df,
        rownames = FALSE,
        options = list(
          scrollX = TRUE,
          pageLength = 10
        ),
        caption = htmltools::tags$caption(
          htmltools::tags$strong(
            "Baseline Sensitivity Analysis: comparison of wavelet-derived metrics across baseline definitions"
          )
        )
      )
    })
    
    
    output$baseline_sensitivity_plot <- renderPlot({
      df <- selected_trace()
      results <- wavelet_results()
      
      sensitivity_df <- build_wavelet_baseline_sensitivity_table(
        time = df$Time,
        signal = df$signal,
        events = results$events,
        event_window = input$event_window,
        initial_fraction = input$initial_f0_fraction,
        minimal_fraction = input$minimal_f0_fraction,
        rolling_window = input$rolling_window,
        rolling_percentile = input$rolling_percentile,
        event_specific_fft_f = input$event_specific_fft_f,
        local_pre_event_window = input$local_pre_event_window,
        local_f0_statistic = input$local_f0_statistic,
        min_amp_above_baseline = input$min_amp_above_baseline,
        accepted_min_ridge_length = input$accepted_min_ridge_length,
        min_cwt_amplitude = input$min_cwt_amplitude,
        summary_event_set = input$summary_event_set
      )
      
      metric_cols <- c(
        "Number_of_Wavelet_Events",
        "Frequency",
        "Mean_Amplitude_Above_Baseline",
        "Mean_DeltaF_over_Baseline_Peak",
        "Mean_Prominence",
        "Mean_FWHP",
        "Mean_FWHM",
        "Mean_Event_AUC",
        "Total_Event_AUC",
        "Total_AUC_Global_Above_Baseline",
        "Mean_CWT_Amplitude",
        "Mean_Ridge_Length"
      )
      
      metric_cols <- metric_cols[metric_cols %in% names(sensitivity_df)]
      
      shiny::validate(
        shiny::need(length(metric_cols) > 0, "No numeric metrics are available for the baseline sensitivity plot.")
      )
      
      plot_long <- do.call(
        rbind,
        lapply(metric_cols, function(m) {
          data.frame(
            Baseline_Method = sensitivity_df$Baseline_Method,
            Metric = m,
            Value = suppressWarnings(as.numeric(sensitivity_df[[m]])),
            stringsAsFactors = FALSE
          )
        })
      )
      
      plot_long <- plot_long[is.finite(plot_long$Value), , drop = FALSE]
      
      shiny::validate(
        shiny::need(nrow(plot_long) > 0, "No finite metric values are available for plotting.")
      )
      
      plot_long$Baseline_Method <- factor(
        plot_long$Baseline_Method,
        levels = sensitivity_df$Baseline_Method
      )
      
      ggplot(plot_long, aes(x = Baseline_Method, y = Value)) +
        geom_col() +
        facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
        labs(
          title = "Sensitivity of wavelet-derived calcium transient metrics to baseline definition",
          subtitle = "Comparison across baseline definitions using the same wavelet-detected candidate events",
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
      df <- selected_trace()
      results <- wavelet_results()
      
      sensitivity_df <- build_wavelet_baseline_sensitivity_table(
        time = df$Time,
        signal = df$signal,
        events = results$events,
        event_window = input$event_window,
        initial_fraction = input$initial_f0_fraction,
        minimal_fraction = input$minimal_f0_fraction,
        rolling_window = input$rolling_window,
        rolling_percentile = input$rolling_percentile,
        event_specific_fft_f = input$event_specific_fft_f,
        local_pre_event_window = input$local_pre_event_window,
        local_f0_statistic = input$local_f0_statistic,
        min_amp_above_baseline = input$min_amp_above_baseline,
        accepted_min_ridge_length = input$accepted_min_ridge_length,
        min_cwt_amplitude = input$min_cwt_amplitude,
        summary_event_set = input$summary_event_set
      )
      
      cv_fun <- function(x) {
        x <- suppressWarnings(as.numeric(x))
        x <- x[is.finite(x)]
        if (length(x) < 2 || mean(abs(x), na.rm = TRUE) == 0) return(NA_real_)
        100 * stats::sd(x, na.rm = TRUE) / mean(abs(x), na.rm = TRUE)
      }
      
      # Keep the same diagnostic logic used in the FFT baseline sensitivity panel:
      # amplitude, AUC and FWHM are the most baseline-dependent metrics.
      cv_amp <- if ("Mean_Amplitude_Above_Baseline" %in% names(sensitivity_df)) {
        cv_fun(sensitivity_df$Mean_Amplitude_Above_Baseline)
      } else NA_real_
      
      cv_auc <- if ("Total_Event_AUC" %in% names(sensitivity_df)) {
        cv_fun(sensitivity_df$Total_Event_AUC)
      } else if ("AUC" %in% names(sensitivity_df)) {
        cv_fun(sensitivity_df$AUC)
      } else NA_real_
      
      cv_fwhm <- if ("Mean_FWHM" %in% names(sensitivity_df)) {
        cv_fun(sensitivity_df$Mean_FWHM)
      } else NA_real_
      
      cv_vector <- c(cv_amp, cv_auc, cv_fwhm)
      max_cv <- if (all(!is.finite(cv_vector))) NA_real_ else max(cv_vector, na.rm = TRUE)
      
      interpretation <- if (!is.finite(max_cv)) {
        "The sensitivity index could not be computed because some metrics are unavailable."
      } else if (max_cv < 10) {
        "The selected calcium trace appears relatively robust to the baseline definition. Biological interpretation is less affected because amplitude, AUC, and FWHM change only slightly across baseline definitions."
      } else if (max_cv < 30) {
        "The selected calcium trace shows moderate baseline sensitivity. Biological interpretation should report the selected baseline definition because amplitude, AUC, or FWHM change noticeably across baseline definitions."
      } else {
        "The selected calcium trace shows strong baseline sensitivity. Biological interpretation should be made carefully because amplitude, AUC, or FWHM change substantially depending on the baseline definition."
      }
      
      fmt_cv <- function(x) {
        if (!is.finite(x)) return("not available")
        paste0(round(x, 2), "%")
      }
      
      tags$div(
        class = "helper-card",
        tags$h4("Automatic interpretation"),
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
          tags$li(paste0("Amplitude CV: ", fmt_cv(cv_amp))),
          tags$li(paste0("AUC CV: ", fmt_cv(cv_auc))),
          tags$li(paste0("FWHM CV: ", fmt_cv(cv_fwhm)))
        )
      )
    })
    
    output$event_table <- DT::renderDT({
      metrics <- wavelet_results()$metrics
      
      if (is.null(metrics) || nrow(metrics) == 0) {
        return(DT::datatable(
          data.frame(Message = "No wavelet events were detected with the selected parameters."),
          options = list(dom = "t")
        ))
      }
      
      numeric_cols <- sapply(metrics, is.numeric)
      metrics[numeric_cols] <- lapply(metrics[numeric_cols], function(z) round(z, 4))
      
      DT::datatable(
        metrics,
        options = list(
          scrollX = TRUE,
          pageLength = 10
        )
      )
    })
    
    output$summary_table <- DT::renderDT({
      metrics <- filtered_metrics_for_summary()
      df <- selected_trace()
      
      baseline_info <- selected_baseline_info()
      
      summary_df <- summarize_wavelet_metrics(
        metrics = metrics,
        time = df$Time,
        signal = df$signal,
        baseline_info = baseline_info
      )
      
      numeric_cols <- sapply(summary_df, is.numeric)
      summary_df[numeric_cols] <- lapply(summary_df[numeric_cols], function(z) round(z, 4))
      
      DT::datatable(
        summary_df,
        options = list(
          dom = "t",
          scrollX = TRUE
        )
      )
    })
    
    
    output$transient_metrics_table <- DT::renderDT({
      metrics <- filtered_metrics_for_summary()
      
      transient_df <- build_wavelet_transient_metrics_table(metrics)
      
      if (nrow(transient_df) == 0) {
        return(DT::datatable(
          data.frame(Message = "No selected wavelet events available for transient metrics."),
          options = list(dom = "t")
        ))
      }
      
      numeric_cols <- sapply(transient_df, is.numeric)
      transient_df[numeric_cols] <- lapply(
        transient_df[numeric_cols],
        function(z) round(z, 4)
      )
      
      DT::datatable(
        transient_df,
        options = list(
          scrollX = TRUE,
          pageLength = 10
        )
      )
    })
    
    output$trace_metrics_table <- DT::renderDT({
      metrics <- filtered_metrics_for_summary()
      df <- selected_trace()
      baseline_info <- selected_baseline_info()
      
      trace_df <- build_wavelet_trace_metrics_table(
        metrics = metrics,
        time = df$Time,
        signal = df$signal,
        baseline_info = baseline_info
      )
      
      numeric_cols <- sapply(trace_df, is.numeric)
      trace_df[numeric_cols] <- lapply(
        trace_df[numeric_cols],
        function(z) round(z, 4)
      )
      
      DT::datatable(
        trace_df,
        options = list(
          dom = "t",
          scrollX = TRUE
        )
      )
    })
    
    
    output$descargarP <- downloadHandler(
      filename = function() "Trace_Metrics.csv",
      content = function(file) {
        metrics <- filtered_metrics_for_summary()
        df <- selected_trace()
        baseline_info <- selected_baseline_info()
        
        trace_df <- build_wavelet_trace_metrics_table(
          metrics = metrics,
          time = df$Time,
          signal = df$signal,
          baseline_info = baseline_info
        )
        
        write.csv(trace_df, file, row.names = FALSE)
      }
    )
    
    output$descargar <- downloadHandler(
      filename = function() "Transient_Metrics.csv",
      content = function(file) {
        metrics <- filtered_metrics_for_summary()
        transient_df <- build_wavelet_transient_metrics_table(metrics)
        
        # Export the same transient table that is shown in the interface.
        write.csv(transient_df, file, row.names = FALSE)
      }
    )
    
    output$Calcium_Trance_Graph <- downloadHandler(
      filename = function() paste0("wavelet_calcium_trace_", Sys.Date(), ".png"),
      content = function(file) {
        ggplot2::ggsave(
          filename = file,
          plot = wavelet_event_plot(),
          dpi = 300,
          width = 10,
          height = 6
        )
      }
    )
    
    output$summary_plot <- renderPlot({
      metrics <- filtered_metrics_for_summary()
      df <- selected_trace()
      baseline_info <- selected_baseline_info()
      
      summary_df <- summarize_wavelet_metrics(
        metrics = metrics,
        time = df$Time,
        signal = df$signal,
        baseline_info = baseline_info
      )
      
      # ------------------------------------------------------------
      # Metric plots by family
      # ------------------------------------------------------------
      # Rationale:
      # The previous plot placed all summary metrics in one bar chart.
      # That was visually simple but statistically misleading because it
      # mixed quantities with different units and scales, such as event
      # counts, time variables, amplitudes, rates, CWT coefficients, and AUC.
      #
      # The new plot separates metrics into conceptual families. This keeps
      # the original numeric values while making the visualization more
      # interpretable and pedagogically consistent with the app.
      # ------------------------------------------------------------
      
      metric_families <- list(
        "Event detection" = c(
          "Number_of_Wavelet_Events",
          "Frequency"
        ),
        "Time-domain metrics" = c(
          "Mean_Onset",
          "Mean_Peak_Occurrence_Time",
          "Mean_Peak_Rise_Time",
          "Mean_Transient_Occurrence_Time",
          "Mean_FWHP",
          "Mean_FWHM"
        ),
        "Signal amplitude metrics" = c(
          "Mean_Amplitude_Above_Baseline",
          "Mean_DeltaF_over_F0_Peak",
          "Mean_Prominence",
          "Mean_Rise_Rate"
        ),
        "Wavelet ridge metrics" = c(
          "Mean_CWT_Amplitude",
          "Mean_Ridge_Length",
          "Mean_Scale"
        ),
        "Area under the curve" = c(
          "Mean_Event_AUC",
          "Total_Event_AUC",
          "Total_AUC_Global_Above_Baseline"
        ),
        "Optional sigmoid fit metrics" = c(
          "Mean_SingleSigmoid_Onset_Time",
          "Mean_SingleSigmoid_Rise_Time",
          "Mean_SingleSigmoid_Rise_Rate",
          "Mean_SingleSigmoid_Slope_at_t0",
          "Mean_SingleSigmoid_R2",
          "Mean_DoubleSigmoid_Decay_Time",
          "Mean_DoubleSigmoid_Decay_Slope",
          "Mean_DoubleSigmoid_Duration"
        )
      )
      
      pretty_metric_names <- c(
        Number_of_Wavelet_Events = "Number of wavelet events",
        Frequency = "Frequency",
        Mean_Onset = "Mean onset time",
        Mean_Peak_Occurrence_Time = "Mean peak occurrence time",
        Mean_Peak_Rise_Time = "Mean peak rise time",
        Mean_Transient_Occurrence_Time = "Mean transient occurrence time",
        Mean_FWHP = "Mean FWHP",
        Mean_FWHM = "Mean FWHM",
        Mean_Amplitude_Above_Baseline = "Mean amplitude above baseline",
        Mean_DeltaF_over_F0_Peak = "Mean DeltaF/Baseline peak",
        Mean_Prominence = "Mean prominence",
        Mean_Rise_Rate = "Mean rise rate",
        Mean_CWT_Amplitude = "Mean CWT amplitude",
        Mean_Ridge_Length = "Mean ridge length",
        Mean_Scale = "Mean scale",
        Mean_Event_AUC = "Mean event AUC",
        Total_Event_AUC = "Total event AUC",
        Total_AUC_Global_Above_Baseline = "Total global AUC above baseline",
        Mean_SingleSigmoid_Onset_Time = "Mean single-sigmoid onset",
        Mean_SingleSigmoid_Rise_Time = "Mean single-sigmoid rise time",
        Mean_SingleSigmoid_Rise_Rate = "Mean single-sigmoid rise rate",
        Mean_SingleSigmoid_Slope_at_t0 = "Mean single-sigmoid slope at t0",
        Mean_SingleSigmoid_R2 = "Mean single-sigmoid R2",
        Mean_DoubleSigmoid_Decay_Time = "Mean double-sigmoid decay time",
        Mean_DoubleSigmoid_Decay_Slope = "Mean double-sigmoid decay slope",
        Mean_DoubleSigmoid_Duration = "Mean double-sigmoid duration"
      )
      
      build_family_df <- function(family_name, cols) {
        cols <- cols[cols %in% names(summary_df)]
        if (length(cols) == 0) {
          return(data.frame(
            Family = character(0),
            Metric = character(0),
            Metric_Label = character(0),
            Value = numeric(0),
            stringsAsFactors = FALSE
          ))
        }
        
        values <- suppressWarnings(as.numeric(summary_df[1, cols]))
        keep <- is.finite(values)
        cols <- cols[keep]
        values <- values[keep]
        
        if (length(cols) == 0) {
          return(data.frame(
            Family = character(0),
            Metric = character(0),
            Metric_Label = character(0),
            Value = numeric(0),
            stringsAsFactors = FALSE
          ))
        }
        
        labels <- ifelse(
          cols %in% names(pretty_metric_names),
          unname(pretty_metric_names[cols]),
          gsub("_", " ", cols)
        )
        
        data.frame(
          Family = family_name,
          Metric = cols,
          Metric_Label = labels,
          Value = values,
          stringsAsFactors = FALSE
        )
      }
      
      plot_df <- do.call(
        rbind,
        Map(build_family_df, names(metric_families), metric_families)
      )
      
      if (is.null(plot_df) || nrow(plot_df) == 0) {
        return(
          ggplot2::ggplot() +
            ggplot2::annotate(
              "text",
              x = 0,
              y = 0,
              label = "No selected wavelet events available for metric plots.",
              size = 5
            ) +
            ggplot2::theme_void()
        )
      }
      
      make_family_plot <- function(family_name) {
        family_df <- plot_df[plot_df$Family == family_name, , drop = FALSE]
        
        if (nrow(family_df) == 0) {
          return(
            ggplot2::ggplot() +
              ggplot2::annotate(
                "text",
                x = 0,
                y = 0,
                label = paste("No available metrics for", family_name),
                size = 4
              ) +
              ggplot2::theme_void()
          )
        }
        
        family_df <- family_df[order(family_df$Value), , drop = FALSE]
        family_df$Metric_Label <- factor(
          family_df$Metric_Label,
          levels = family_df$Metric_Label
        )
        
        ggplot2::ggplot(
          family_df,
          ggplot2::aes(x = Metric_Label, y = Value)
        ) +
          ggplot2::geom_col(width = 0.70) +
          ggplot2::coord_flip() +
          ggplot2::labs(
            title = family_name,
            x = NULL,
            y = "Value"
          ) +
          ggplot2::theme_classic(base_size = 12) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13),
            axis.text.y = ggplot2::element_text(size = 10),
            axis.title.x = ggplot2::element_text(size = 10),
            plot.margin = ggplot2::margin(t = 5, r = 15, b = 5, l = 5)
          )
      }
      
      plot_list <- lapply(names(metric_families), make_family_plot)
      
      gridExtra::grid.arrange(
        grobs = plot_list,
        ncol = 1,
        top = grid::textGrob(
          paste0(
            "Wavelet Ridgewalking: Summary Metrics by Family\n",
            "Events used in summary metrics: ", input$summary_event_set
          ),
          gp = grid::gpar(fontsize = 15, fontface = "bold")
        )
      )
    })
    
    return(wavelet_results)
  })
}
