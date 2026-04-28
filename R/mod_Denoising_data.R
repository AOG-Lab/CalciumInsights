# Adapted from the standalone FFT-based app.R to work inside a golem app.
# Keep the function names used by app_ui.R: mod_Denoising_data_ui/server.

# app.R
# ============================================================
# FFT-based Denoising Data App
# - Preserves the original denoising workflow
# - Replaces LOESS with FFT low-pass smoothing
# - Adds Fourier explicit reconstruction formula using kept frequencies
# - Uses helper functions saved in the same folder
# ============================================================

library(shiny)
library(shinyjs)
library(DT)
library(ggplot2)
library(latex2exp)
library(prospectr)
library(gridExtra)
library(pracma)
library(vroom)
library(jsonlite)
library(dplyr)

# =========================
# Source helper functions
# =========================
# source("utils_load_file.R")
# source("utils_peaks.R")
# source("utils_Time_of_the_first_peak.R")
# source("utils_prominens2.R")
# source("utils_FWHP2.R")
# source("utils_response_time.R")
# source("utils_right_left_FWHP.R")
# source("utils_AUC2.R")
# source("utils_Savitzky_Golay.R")

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
                              baseline_mode = 1,
                              lim_inf = 0,
                              lim_sup = 20,
                              own_baseline = 0,
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
      if (!is.null(time_start_increasin_peak) &&
          is.data.frame(time_start_increasin_peak) &&
          "Time" %in% names(time_start_increasin_peak) &&
          length(time_start_increasin_peak$Time) > 0 &&
          is.finite(suppressWarnings(as.numeric(time_start_increasin_peak$Time[1])))) {
        Time_One_set <- suppressWarnings(as.numeric(time_start_increasin_peak$Time[1]))
        posicion_Time_One_set <- which.min(abs(time_vals - Time_One_set))
        if (length(posicion_Time_One_set) == 0 || !is.finite(posicion_Time_One_set)) {
          safe_mean(signal_vals)
        } else {
          safe_mean(signal_vals[seq_len(posicion_Time_One_set)])
        }
      } else {
        safe_mean(signal_vals)
      }
    },
    "3" = {
      lim_inf <- safe_scalar(lim_inf, fallback = min(time_vals, na.rm = TRUE))
      lim_sup <- safe_scalar(lim_sup, fallback = max(time_vals, na.rm = TRUE))
      df_filtrado <- data_smoothed[time_vals >= lim_inf & time_vals <= lim_sup, , drop = FALSE]
      if (nrow(df_filtrado) > 0) safe_mean(df_filtrado$signal) else safe_mean(signal_vals)
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
    safe_mean(signal_vals)
  )

  scalar <- safe_scalar(scalar, fallback = safe_mean(signal_vals))
  list(scalar = scalar, trace = baseline_trace, mode = baseline_mode)
}

baseline_at_times <- function(time_points, baseline_info) {
  time_points <- as.numeric(time_points)
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
    baseline_mode = 1,
    lim_inf = 0,
    lim_sup = 20,
    own_baseline = 0,
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
      if (!is.null(time_start_increasin_peak) && "Time" %in% names(time_start_increasin_peak) &&
          length(time_start_increasin_peak$Time) > 0 && is.finite(time_start_increasin_peak$Time[1])) {
        Time_One_set <- time_start_increasin_peak$Time[1]
        posicion_Time_One_set <- which.min(abs(time_vals - Time_One_set))
        mean(signal_vals[1:posicion_Time_One_set], na.rm = TRUE)
      } else {
        mean(signal_vals, na.rm = TRUE)
      }
    },
    "3" = {
      if (!is.finite(lim_inf)) lim_inf <- min(time_vals, na.rm = TRUE)
      if (!is.finite(lim_sup)) lim_sup <- max(time_vals, na.rm = TRUE)
      df_filtrado <- data_smoothed[time_vals >= lim_inf & time_vals <= lim_sup, , drop = FALSE]
      if (nrow(df_filtrado) > 0) mean(df_filtrado$signal, na.rm = TRUE) else mean(signal_vals, na.rm = TRUE)
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
    mean(signal_vals, na.rm = TRUE)
  )
}

extract_peak_metrics_from_smoothed <- function(
    data_smoothed,
    nups,
    ndowns,
    minpeakheight,
    minpeakdistance,
    baseline_mode = 1,
    lim_inf = 0,
    lim_sup = 20,
    own_baseline = 0,
    min_FWHP = 0,
    min_prominence = 0
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
      empty_result$baseline <- compute_baseline_value_from_mode(
        data_smoothed = data_smoothed,
        baseline_mode = baseline_mode,
        lim_inf = lim_inf,
        lim_sup = lim_sup,
        own_baseline = own_baseline,
        time_start_increasin_peak = NULL
      )
      return(empty_result)
    }

    table_peak <- as.data.frame(peaks_found$p_eak)
    table_positions_peaks <- as.data.frame(peaks_found$peak)
    peaks_idx <- table_positions_peaks[, 2]

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

    baseline1 <- compute_baseline_value_from_mode(
      data_smoothed = data_smoothed,
      baseline_mode = baseline_mode,
      lim_inf = lim_inf,
      lim_sup = lim_sup,
      own_baseline = own_baseline,
      time_start_increasin_peak = time_start_increasin_peak
    )

    p_eak_mediun <- c((table_positions_peaks[, 1] + baseline1) / 2)
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
    if (is.finite(baseline1)) {
      df_p <- df_p[df_p$Amplitude > baseline1, , drop = FALSE]
    }

    if (nrow(df_p) > 0) {
      df_p <- df_p[order(df_p$Peak_Occurence_Time), , drop = FALSE]
      mean_peak_vals <- df_p$Amplitude
      if (is.finite(baseline1)) {
        mean_peak_vals <- mean_peak_vals - baseline1
      }
      summary_df <- data.frame(
        n_peaks = nrow(df_p),
        mean_peak = mean(mean_peak_vals, na.rm = TRUE),
        mean_peak_occurrence = mean(df_p$Peak_Occurence_Time, na.rm = TRUE),
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

# =========================
# UI Module
# =========================
mod_Denoising_data_ui <- function(id) {
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

        actionButton(
          ns("param_info_button11"), "Help",
          class = "btn-sm",
          style = "position:absolute; top:0; right:15px; margin:5px;"
        ),

        radioButtons(
          ns("data_simulate"), "Example Data",
          choices = c("Yes" = 1, "No" = 0),
          selected = 1
        ),

        fileInput(
          ns("fileBcsv2"),
          accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".tsv", ".json"),
          label = h5("Dataset")
        ),

        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),

        numericInput(ns("Cell2"), "Region of Interest (ROI):", value = 1, min = 1),
        sliderInput(ns("fft_fraction"), "FFT Low-Frequency Fraction:", min = 0.01, max = 1,
                    value = 0.20, step = 0.01),
        checkboxInput(ns("keep_mean"), "Keep mean (DC component)", TRUE),

        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),

        tags$h4("Find Peaks Function Arguments", style = "color: gray; margin-top: 10px;"),

        fluidRow(
          column(
            width = 6,
            numericInput(ns("minpeakheight2"), "1. Peak Height (min)", value = 0, min = 0, max = 100, step = 0.1),
            numericInput(ns("ndowns2"), "3. Peak Descent", value = 1, min = 0, max = 100)
          ),
          column(
            width = 6,
            numericInput(ns("nups2"), "2. Peak Ascent", value = 1, min = 0, max = 100),
            numericInput(ns("minpeakdistance2"), "4. Min Peak Distance", value = 1, min = 0, max = 100)
          ),
          column(
            width = 6,
            numericInput(ns("min_FWHP"), "5. FWHP (min)", value = 0, min = 0, step = 0.1)
          ),
          column(
            width = 6,
            numericInput(ns("min_prominence"), "6. Prominence (min)", value = 0, min = 0, step = 0.1)
          )
        ),

        div(style = "border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;"),

        selectInput(
          ns("Baseline"), "Baseline:",
          choices = c(
            "Reference Level 0" = 1,
            "Standard definition" = 2,
            "Interval" = 3,
            "Your baseline" = 4,
            "Min" = 5,
            "Rolling percentile baseline" = 6,
            "Fixed-window percentile baseline" = 7
          ),
          selected = 1
        ),

        conditionalPanel(
          condition = "input.Baseline==3", ns = ns,
          fluidRow(
            column(6, textInput(ns("Lim_inf"), "Lim inf:", value = "0")),
            column(6, textInput(ns("Lim_sup"), "Lim sup:", value = "20"))
          )
        ),

        conditionalPanel(
          condition = "input.Baseline==4", ns = ns,
          numericInput(ns("own_baseline"), "Own Baseline:", value = 0, step = 0.1)
        ),

        conditionalPanel(
          condition = "input.Baseline==6 || input.Baseline==7", ns = ns,
          fluidRow(
            column(
              6,
              numericInput(
                ns("rolling_window"),
                "Moving window size:",
                value = 600,
                min = 0.0001,
                step = 1
              )
            ),
            column(
              6,
              numericInput(
                ns("rolling_percentile"),
                "Percentile:",
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
            "Show fixed-window intervals on calcium trace",
            value = TRUE
          )

        ),

        radioButtons(ns("auc2"), "Area Under the Curve (AUC):", choices = c("No" = 1, "Yes" = 2), selected = 1),
        radioButtons(ns("raw_data"), "Raw Data", choices = c("No" = 1, "Yes" = 2), selected = 2),
        radioButtons(ns("FWHM"), "Full Width at Half Maximum:", choices = c("No" = 1, "Yes" = 2), selected = 1),

        downloadButton(ns("descargarP"), "Trace Metrics"),
        downloadButton(ns("descargar"), "Transient Metrics"),
        downloadButton(ns("Calcium_Trance_Graph"), "Calcium Trace Graph")
      ),

      mainPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "SummaryData",
            DTOutput(ns("infotable2")),
            DTOutput(ns("data2"))
          ),
          tabPanel(
            "Peaks",
            tabsetPanel(
              type = "tabs",
              tabPanel(
                "Metrics",
                DTOutput(ns("table_peaks2")),
                DTOutput(ns("table_peaks22"))
              ),
              tabPanel(
                "Metric plots",
                plotOutput(ns("plot_peak3")),
                plotOutput(ns("derivative")),
                plotOutput(ns("plot_raw_smoothed"))
              ),
              tabPanel(
                "Fourier analysis",
                fluidRow(
                  column(
                    width = 6,
                    h4("FFT Summary"),
                    verbatimTextOutput(ns("fft_infoText")),
                    plotOutput(ns("fft_fitPlot"), height = "320px"),
                    plotOutput(ns("fft_residPlot"), height = "260px")
                  ),
                  column(
                    width = 6,
                    h4("Spectrum and Frequency Table"),
                    plotOutput(ns("fft_specPlot"), height = "320px"),
                    DTOutput(ns("fft_freqTable"))
                  )
                ),
                fluidRow(
                  column(
                    width = 12,
                    tags$div(
                      class = "helper-card",
                      h4("Explicit reconstruction formula using the kept FFT frequencies"),
                      tags$p("This formula writes the recovered series explicitly as a finite sum of cosines and sines using only the frequencies retained by the FFT low-pass filter."),
                      div(class = "formula-box", textOutput(ns("fft_formula_text")))
                    )
                  )
                ),
                fluidRow(
                  column(
                    width = 12,
                    h4("Theory and Interpretation"),
                    uiOutput(ns("fft_formula_box"))
                  )
                )
              )
            )
          ),
          tabPanel(
            "FFT selection criteria",
            tags$div(
              class = "helper-card",
              h4("Selection of FFT Low-Frequency Fraction"),
              tags$p("This panel evaluates the candidate fractions f = 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, and 0.40 using four descriptive criteria."),
              tags$p("Specifically, the assessment considers: (1) preservation of the original signal morphology after FFT-based smoothing, (2) the structure and magnitude of the residual signal, (3) the fraction of spectral power retained in the periodogram, and (4) the closeness of each candidate mean peak amplitude (mean_peak) to a weighted reference mean computed across all candidate fractions."),
              tags$p("The weighted reference mean of mean_peak is calculated using the number of detected peaks as weights. Based on Criterion 4, the app proposes a candidate fraction f whose mean_peak is closest to this weighted reference mean. However, the final selection of f remains at the user's discretion.")
            ),
            # tags$div(
            #   class = "helper-card",
            #   h4("Fixed parameters used in this panel"),
            #   uiOutput(ns("fft_grid_current_params"))
            # ),
            uiOutput(ns("fft_grid_recommendation")),
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
mod_Denoising_data_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    observeEvent(input$param_info_button11, {
      showModal(modalDialog(
        title = "Help",
        size = "l",
        HTML("<p style='text-align:justify;'><strong>This module reproduces the original denoising workflow, but the smoothing step is now based on FFT low-pass filtering instead of LOESS.</strong></p>"),
        HTML("<p style='text-align:justify;'><strong>FFT Low-Frequency Fraction:</strong> controls how many low frequencies are preserved. Smaller values produce stronger smoothing; values near 1 preserve more of the original signal.</p>"),
        HTML("<p style='text-align:justify;'><strong>Keep mean (DC component):</strong> keeps the global average level of the signal during FFT reconstruction.</p>"),
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

    peaks_df <- reactive({
      req(filedata()$fileInput)
      data <- filedata()$fileInput
      shiny::validate(shiny::need(ncol(data) >= 2, "Dataset must contain a time column and at least one ROI."))

      time <- as.numeric(data[[1]])
      roi_index <- as.numeric(input$Cell2) + 1
      shiny::validate(shiny::need(roi_index <= ncol(data), "Selected ROI is out of range."))

      data_raw <- data.frame(Time = time, signal = as.numeric(data[[roi_index]]))

      fft_res <- fft_lowpass_smooth_trace(
        time = data_raw$Time,
        signal = data_raw$signal,
        f = input$fft_fraction,
        keep_mean = isTRUE(input$keep_mean)
      )
      df_smoothed <- fft_res$data

      peaks_found <- peaks(
        data = df_smoothed,
        nups = input$nups2,
        ndowns = input$ndowns2,
        minpeakheight = input$minpeakheight2,
        minpeakdistance = input$minpeakdistance2
      )

      list(
        table_peak = peaks_found$p_eak,
        table_positions_peaks = peaks_found$peak,
        data_raw = data_raw,
        df_smoothed = df_smoothed,
        data_matrix = as.matrix(data[, -1, drop = FALSE]),
        fft_res = fft_res
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
      baseline_mode = "1",
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

      df_raw <- data.frame(Time = time, signal = as.numeric(data[[roi_index]]))
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

        peak_res <- extract_peak_metrics_from_smoothed(
          data_smoothed = df_smoothed,
          nups = fft_grid_fixed_params$nups,
          ndowns = fft_grid_fixed_params$ndowns,
          minpeakheight = fft_grid_fixed_params$minpeakheight,
          minpeakdistance = fft_grid_fixed_params$minpeakdistance,
          baseline_mode = fft_grid_fixed_params$baseline_mode,
          lim_inf = lim_inf_num,
          lim_sup = lim_sup_num,
          own_baseline = fft_grid_fixed_params$own_baseline,
          min_FWHP = fft_grid_fixed_params$min_FWHP,
          min_prominence = fft_grid_fixed_params$min_prominence
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
          peak_res = peak_res,
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
      peaks_idx <- table_positions_peaks[, 2]

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
        time_start_increasin_peak = time_start_increasin_peak,
        rolling_window = input$rolling_window,
        rolling_percentile = input$rolling_percentile
      )
      baseline1 <- baseline_info$scalar
      baseline_at_peak <- baseline_at_times(table_peak$posision_peaks, baseline_info)
      table_peak$Baseline_at_peak <- baseline_at_peak

      if (identical(as.character(input$auc2), "2")) {
        if (!is.null(baseline_info$trace)) {
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
             x = "Raw [Delta F/F0]",
             y = "Smoothed [Delta F/F0]") +
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
             y = "Delta F/F0 * s^-1") +
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

      peaks_idx <- table_positions_peaks[, 2]
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
      df_FWHM2 <- df_FWHM2[df_FWHM2$FWHP > input$min_FWHP, , drop = FALSE]
      df_FWHM2 <- df_FWHM2[df_FWHM2$Prominence > input$min_prominence, , drop = FALSE]
      if ("Baseline_at_peak" %in% names(df_FWHM2)) {
        df_FWHM2 <- df_FWHM2[df_FWHM2$Amplitude > df_FWHM2$Baseline_at_peak, , drop = FALSE]
      } else if (is.finite(peaks_plot()$baseline1)) {
        df_FWHM2 <- df_FWHM2[df_FWHM2$Amplitude > peaks_plot()$baseline1, , drop = FALSE]
      }
      list(df_p = df_FWHM2)
    })

    Trance_Graph <- reactive({
      data_smoothed <- peaks_df()$df_smoothed
      data_raw <- peaks_df()$data_raw
      colnames(data_smoothed) <- c("Time", "Sing")
      df_p <- Peaks_Data_Final()$df_p

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
                                      y = Prominence_Midpoint, yend = Prominence_Midpoint),
                     linetype = "solid", linewidth = 1.2, color = "orange") +
        labs(title = "Calcium Trace",
             x = "Time [s]",
             y = "Delta F/F0") +
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

      if (identical(as.character(input$auc2), "2")) {
        Integration_Reference <- peaks_plot()$baseline1
        auc_df <- data_smoothed
        baseline_trace <- peaks_plot()$baseline_trace

        if (!is.null(baseline_trace) && nrow(baseline_trace) > 0) {
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

      if (as.character(input$Baseline) %in% c("6", "7") && !is.null(peaks_plot()$baseline_trace)) {
        gg3 <- gg3 +
          geom_line(
            data = peaks_plot()$baseline_trace,
            aes(x = Time, y = Baseline),
            linetype = "dashed",
            linewidth = 1.2,
            color = "darkgreen",
            inherit.aes = FALSE
          )
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
                y = Amplitude_Midpoint, yend = Amplitude_Midpoint),
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

    output$table_peaks2 <- renderDT({
      df_p <- Peaks_Data_Final()$df_p
      shiny::validate(shiny::need(nrow(df_p) > 0, "No peaks satisfy the selected thresholds."))
      df_p <- subset(df_p, select = c("Amplitude", "Baseline_at_peak", "Peak_Occurence_Time", "Prominence", "FWHP", "FWHM",
                                      "Peak_Rise_Time", "Transient_Ocurrence_Time", "Rise_Rate"))
      df_p$Amplitude <- df_p$Amplitude - df_p$Baseline_at_peak
      datatable(round(df_p, 3), caption = tags$caption(tags$strong("Transient Metrics")))
    })

    output$table_peaks22 <- renderDT({
      df_p <- Peaks_Data_Final()$df_p
      shiny::validate(shiny::need(nrow(df_p) > 0, "No trace metrics are available because no valid peaks were detected."))
      time1 <- min(peaks_df()$data_raw$Time)
      time2 <- max(peaks_df()$data_raw$Time)
      Time_OnSet <- df_p$Transient_Ocurrence_Time[1]
      Frequency <- length(df_p$Amplitude) / (time2 - time1)
      Baseline <- peaks_plot()$baseline1
      number_of_peaks <- length(df_p$Amplitude)
      df_p2 <- data.frame(Time_Onset = Time_OnSet, Frequency = Frequency, Baseline = Baseline,
                          Number_of_Peaks = number_of_peaks)

      if (identical(as.character(input$auc2), "2")) {
        Transient_Metrics <- cbind(df_p2, peaks_plot()$tabla_AUC)
      } else {
        Transient_Metrics <- df_p2
      }

      datatable(round(Transient_Metrics, 3), options = list(pagingType = "simple", dom = "t", autoWidth = TRUE),
                caption = tags$caption(tags$strong("Trace Metrics")))
    })

    output$descargar <- downloadHandler(
      filename = function() "Transient_Metrics.csv",
      content = function(file) {
        df_p <- Peaks_Data_Final()$df_p
        shiny::validate(shiny::need(nrow(df_p) > 0, "No transient metrics available to export."))
        time1 <- min(peaks_df()$data_raw$Time)
        time2 <- max(peaks_df()$data_raw$Time)
        Time_OnSet <- df_p$Transient_Ocurrence_Time[1]
        Frequency <- length(df_p$Amplitude) / (time2 - time1)
        Baseline <- peaks_plot()$baseline1
        number_of_peaks <- length(df_p$Amplitude)
        df_p2 <- data.frame(id = 1, Time_Onset = Time_OnSet, Frequency = Frequency,
                            Baseline = Baseline, Number_of_Peaks = number_of_peaks)
        out <- if (identical(as.character(input$auc2), "2")) cbind(df_p2, peaks_plot()$tabla_AUC) else df_p2
        write.csv(out, file, row.names = FALSE)
      }
    )

    output$descargarP <- downloadHandler(
      filename = function() "Trace_Metrics.csv",
      content = function(file) {
        df_p <- Peaks_Data_Final()$df_p
        shiny::validate(shiny::need(nrow(df_p) > 0, "No trace metrics available to export."))
        df_p$Amplitude <- df_p$Amplitude - df_p$Baseline_at_peak
        df_p$id <- seq_len(nrow(df_p))
        df_p <- subset(df_p, select = c("id", "Amplitude", "Baseline_at_peak", "Peak_Occurence_Time", "Prominence", "FWHP",
                                        "FWHM", "Peak_Rise_Time", "Transient_Ocurrence_Time", "Rise_Rate"))
        write.csv(df_p, file, row.names = FALSE)
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

    output$fft_grid_current_params <- renderUI({
      baseline_label <- switch(
        as.character(fft_grid_fixed_params$baseline_mode),
        "1" = "Reference Level 0",
        "2" = "Standard definition",
        "3" = paste0("Interval [", fft_grid_fixed_params$lim_inf, ", ", fft_grid_fixed_params$lim_sup, "]"),
        "4" = paste0("Your baseline (", fft_grid_fixed_params$own_baseline, ")"),
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
      df <- fft_fraction_grid()$summary_df
      df_show <- df[, c(
        "f", "k_keep", "corr_raw_smooth", "retained_power_frac",
        "n_peaks", "mean_peak"
      )]

      round_cols <- c("f", "corr_raw_smooth", "retained_power_frac", "mean_peak")
      for (nm in round_cols) {
        if (nm %in% names(df_show)) df_show[[nm]] <- round(df_show[[nm]], 4)
      }

      if ("mean_peak" %in% names(df_show)) {
        df_show$mean_peak <- ifelse(is.na(df_show$mean_peak), "—", sprintf("%.4f", as.numeric(df_show$mean_peak)))
      }

      datatable(
        df_show,
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 7),
        caption = tags$caption(tags$strong("Grid assessment for candidate FFT fractions"))
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
          y = "Delta F/F0"
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

