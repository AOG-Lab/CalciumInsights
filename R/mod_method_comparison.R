# This comparison module is retained from the current app source but remains inactive in app_ui/app_server.
# R/mod_method_comparison.R
# ============================================================
# Method Comparison Module
# CalciumInsights
#
# Purpose:
#   This module provides a professional and didactic comparison panel
#   for interpreting the FFT + Baseline Analysis and Wavelet Ridgewalking
#   workflows.
#
# Current version:
#   - Works as a complete standalone explanatory module.
#   - Provides structured interpretation rules.
#   - Provides comparison tables and recommended reporting language.
#   - Includes an optional future-ready server structure to receive
#     reactive outputs from FFT and Wavelet modules.
#
# Future-ready usage:
#   fft_results     <- mod_Denoising_data_fft_server("fft_module")
#   wavelet_results <- mod_wavelet_ridgewalking_server("wavelet_module")
#   mod_method_comparison_server(
#     "comparison_module",
#     fft_results = fft_results,
#     wavelet_results = wavelet_results
#   )
# ============================================================


# ============================================================
# UI
# ============================================================

mod_method_comparison_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tagList(
    
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        .ci-page {
          padding-top: 10px;
        }

        .ci-title {
          font-size: 28px;
          font-weight: 700;
          color: #2f2f2f;
          margin-bottom: 6px;
        }

        .ci-subtitle {
          font-size: 15px;
          color: #555555;
          margin-bottom: 18px;
          line-height: 1.55;
        }

        .ci-card {
          background: #ffffff;
          border: 1px solid #e5e5e5;
          border-radius: 12px;
          padding: 18px 20px;
          margin-bottom: 18px;
          box-shadow: 0 1px 3px rgba(0,0,0,0.04);
        }

        .ci-card-soft {
          background: #f7f9fb;
          border: 1px solid #e5e5e5;
          border-radius: 12px;
          padding: 18px 20px;
          margin-bottom: 18px;
        }

        .ci-card h3,
        .ci-card-soft h3 {
          margin-top: 0;
          margin-bottom: 10px;
          font-size: 20px;
          font-weight: 700;
          color: #333333;
        }

        .ci-card h4,
        .ci-card-soft h4 {
          margin-top: 14px;
          margin-bottom: 8px;
          font-size: 16px;
          font-weight: 700;
          color: #333333;
        }

        .ci-card p,
        .ci-card-soft p {
          font-size: 14.5px;
          line-height: 1.6;
          color: #333333;
        }

        .ci-grid-3 {
          display: grid;
          grid-template-columns: repeat(3, 1fr);
          gap: 14px;
          margin-top: 12px;
        }

        .ci-grid-2 {
          display: grid;
          grid-template-columns: repeat(2, 1fr);
          gap: 14px;
          margin-top: 12px;
        }

        .ci-mini-card {
          background: #ffffff;
          border: 1px solid #e6e6e6;
          border-radius: 10px;
          padding: 14px;
        }

        .ci-mini-card h4 {
          margin: 0 0 8px 0;
          font-size: 15px;
          font-weight: 700;
          color: #2f2f2f;
        }

        .ci-mini-card p {
          margin: 0;
          font-size: 13.5px;
          line-height: 1.55;
          color: #444444;
        }

        .ci-badge {
          display: inline-block;
          padding: 4px 9px;
          border-radius: 999px;
          font-size: 12px;
          font-weight: 700;
          background: #eaf3ff;
          color: #337ab7;
          margin-bottom: 8px;
        }

        .ci-note {
          background: #fff8e8;
          border-left: 5px solid #f0ad4e;
          border-radius: 10px;
          padding: 14px 16px;
          margin-top: 12px;
          font-size: 14px;
          line-height: 1.55;
        }

        .ci-good {
          background: #edf7ed;
          border-left: 5px solid #5cb85c;
          border-radius: 10px;
          padding: 14px 16px;
          margin-top: 12px;
          font-size: 14px;
          line-height: 1.55;
        }

        .ci-warning {
          background: #fff3f3;
          border-left: 5px solid #d9534f;
          border-radius: 10px;
          padding: 14px 16px;
          margin-top: 12px;
          font-size: 14px;
          line-height: 1.55;
        }

        .ci-code {
          background: #f7f7f7;
          border: 1px solid #e5e5e5;
          border-radius: 8px;
          padding: 12px;
          font-size: 13px;
          overflow-x: auto;
          white-space: pre;
        }

        .ci-status-box {
          background: #f7f9fb;
          border: 1px dashed #c9d6e2;
          border-radius: 10px;
          padding: 14px 16px;
          margin-bottom: 15px;
          color: #444444;
          font-size: 14px;
          line-height: 1.55;
        }

        @media (max-width: 900px) {
          .ci-grid-3,
          .ci-grid-2 {
            grid-template-columns: 1fr;
          }
        }
      "))
    ),
    
    shiny::fluidPage(
      class = "ci-page",
      
      shiny::div(
        class = "ci-card",
        shiny::div(class = "ci-badge", "CalciumInsights"),
        shiny::h2(class = "ci-title", "Method Comparison"),
        shiny::p(
          class = "ci-subtitle",
          "This module explains how to interpret results from the FFT + Baseline Analysis and Wavelet Ridgewalking workflows. ",
          "The goal is not to declare one method universally superior, but to evaluate whether calcium transient conclusions are robust across complementary detection strategies."
        )
      ),
      
      shiny::tabsetPanel(
        id = ns("comparison_tabs"),
        type = "tabs",
        
        # ------------------------------------------------------------
        # TAB 1: Overview
        # ------------------------------------------------------------
        shiny::tabPanel(
          title = "Overview",
          value = "overview",
          
          shiny::br(),
          
          shiny::div(
            class = "ci-card-soft",
            shiny::h3("Why compare methods?"),
            shiny::p(
              "Calcium transient analysis depends on several methodological decisions: smoothing, baseline definition, event detection, ",
              "peak filtering, and metric extraction. Different pipelines may detect slightly different events, especially when signals are noisy, ",
              "broad, low-amplitude, or baseline-sensitive."
            ),
            shiny::p(
              "The FFT + Baseline Analysis module focuses on a smoothed signal and peak detection. ",
              "The Wavelet Ridgewalking module detects events as time-scale structures in a continuous wavelet transform map. ",
              "Together, they provide a stronger framework for evaluating robustness."
            )
          ),
          
          shiny::div(
            class = "ci-grid-3",
            
            shiny::div(
              class = "ci-mini-card",
              shiny::h4("1. FFT + Baseline Analysis"),
              shiny::p(
                "Best suited for clear peaks after smoothing. It provides detailed baseline-dependent transient metrics such as amplitude, AUC, FWHM, FWHP, prominence, and frequency."
              )
            ),
            
            shiny::div(
              class = "ci-mini-card",
              shiny::h4("2. Wavelet Ridgewalking"),
              shiny::p(
                "Useful for detecting events as multiscale structures. It can help identify broader or lower-amplitude transients that may not be captured by a fixed peak-detection setting."
              )
            ),
            
            shiny::div(
              class = "ci-mini-card",
              shiny::h4("3. Method Comparison"),
              shiny::p(
                "Used to determine whether the biological interpretation is consistent across methods, or whether conclusions depend strongly on the selected detection strategy."
              )
            )
          ),
          
          shiny::div(
            class = "ci-note",
            shiny::strong("Recommended principle: "),
            "If both methods support the same biological conclusion, the result is more robust. ",
            "If the methods disagree, the disagreement should be reported and interpreted as method sensitivity rather than ignored."
          )
        ),
        
        # ------------------------------------------------------------
        # TAB 2: Conceptual comparison
        # ------------------------------------------------------------
        shiny::tabPanel(
          title = "Comparison matrix",
          value = "comparison_matrix",
          
          shiny::br(),
          
          shiny::div(
            class = "ci-card",
            shiny::h3("Conceptual comparison between workflows"),
            shiny::p(
              "This table summarizes the practical differences between the two detection strategies implemented in CalciumInsights."
            ),
            DT::DTOutput(ns("conceptual_comparison_table"))
          )
        ),
        
        # ------------------------------------------------------------
        # TAB 3: Interpretation guide
        # ------------------------------------------------------------
        shiny::tabPanel(
          title = "Interpretation guide",
          value = "interpretation_guide",
          
          shiny::br(),
          
          shiny::div(
            class = "ci-card",
            shiny::h3("How to interpret agreement or disagreement"),
            shiny::p(
              "Use this guide when reporting whether the FFT + Baseline Analysis and Wavelet Ridgewalking modules produce consistent conclusions."
            ),
            DT::DTOutput(ns("interpretation_table"))
          ),
          
          shiny::div(
            class = "ci-good",
            shiny::strong("Robust conclusion: "),
            "When both methods detect a similar number of events, similar timing, and similar overall signal behavior, the calcium transient pattern can be considered more reliable."
          ),
          
          shiny::div(
            class = "ci-warning",
            shiny::strong("Sensitive conclusion: "),
            "When event counts, timing, amplitude, AUC, or frequency differ strongly between methods, the interpretation should explicitly state that the result is method-sensitive."
          )
        ),
        
        # ------------------------------------------------------------
        # TAB 4: Metrics to compare
        # ------------------------------------------------------------
        shiny::tabPanel(
          title = "Metrics to compare",
          value = "metrics_to_compare",
          
          shiny::br(),
          
          shiny::div(
            class = "ci-card",
            shiny::h3("Recommended metrics for side-by-side comparison"),
            shiny::p(
              "The following groups of metrics are the most useful for comparing FFT-based peak detection and wavelet ridgewalking."
            ),
            DT::DTOutput(ns("recommended_metrics_table"))
          ),
          
          shiny::div(
            class = "ci-note",
            shiny::strong("Important: "),
            "Do not compare all metrics as if they had the same meaning or scale. ",
            "Time metrics, amplitude metrics, frequency metrics, and AUC metrics should be interpreted separately."
          )
        ),
        
        # ------------------------------------------------------------
        # TAB 5: Direct comparison
        # ------------------------------------------------------------
        shiny::tabPanel(
          title = "Direct comparison",
          value = "direct_comparison",
          
          shiny::br(),
          
          shiny::div(
            class = "ci-card",
            shiny::h3("Direct comparison output"),
            shiny::uiOutput(ns("comparison_status")),
            DT::DTOutput(ns("direct_comparison_table")),
            shiny::br(),
            shiny::plotOutput(ns("direct_comparison_plot"), height = "620px")
          )
        ),
        
        # ------------------------------------------------------------
        # TAB 6: Reporting guide
        # ------------------------------------------------------------
        shiny::tabPanel(
          title = "Reporting guide",
          value = "reporting_guide",
          
          shiny::br(),
          
          shiny::div(
            class = "ci-card",
            shiny::h3("Suggested language for scientific reporting"),
            shiny::p(
              "The following reporting statements can be adapted depending on whether both methods agree or disagree."
            ),
            
            shiny::h4("When both methods agree"),
            shiny::div(
              class = "ci-code",
              "Both the FFT-based peak detection workflow and the wavelet ridgewalking workflow identified a similar calcium transient pattern. This agreement suggests that the main conclusion is robust to the event-detection strategy."
            ),
            
            shiny::h4("When FFT detects more events"),
            shiny::div(
              class = "ci-code",
              "The FFT-based workflow detected a larger number of events than the wavelet ridgewalking workflow. This may indicate that the selected peak-detection parameters split complex transients into multiple peaks or that the result is sensitive to amplitude-threshold settings."
            ),
            
            shiny::h4("When Wavelet detects more events"),
            shiny::div(
              class = "ci-code",
              "The wavelet ridgewalking workflow detected additional events not identified by the FFT-based workflow. These events may correspond to broad, low-amplitude, or multiscale calcium transients that are less easily detected using fixed peak-detection parameters."
            ),
            
            shiny::h4("When baseline sensitivity is strong"),
            shiny::div(
              class = "ci-code",
              "Amplitude- and area-based metrics changed substantially across baseline definitions. Therefore, quantitative conclusions based on amplitude, AUC, or FWHM should be interpreted as baseline-sensitive and reported with caution."
            )
          )
        ),
        
        # ------------------------------------------------------------
        # TAB 7: Implementation note
        # ------------------------------------------------------------
        shiny::tabPanel(
          title = "Implementation note",
          value = "implementation_note",
          
          shiny::br(),
          
          shiny::div(
            class = "ci-card",
            shiny::h3("Future-ready implementation"),
            shiny::p(
              "This module currently works as a professional interpretation and comparison panel. ",
              "It is also written so that later you can pass reactive outputs from the FFT and Wavelet modules."
            ),
            
            shiny::h4("Current app.R usage"),
            shiny::div(
              class = "ci-code",
              "mod_method_comparison_server(\"comparison_module\")"
            ),
            
            shiny::h4("Future app.R usage with connected modules"),
            shiny::div(
              class = "ci-code",
              "fft_results <- mod_Denoising_data_fft_server(\"fft_module\")\nwavelet_results <- mod_wavelet_ridgewalking_server(\"wavelet_module\")\n\nmod_method_comparison_server(\n  \"comparison_module\",\n  fft_results = fft_results,\n  wavelet_results = wavelet_results\n)"
            ),
            
            shiny::div(
              class = "ci-note",
              shiny::strong("Next programming step: "),
              "Modify the FFT and Wavelet server modules so each returns a named reactive list containing at least: ",
              shiny::code("summary_metrics"),
              ", ",
              shiny::code("event_table"),
              ", and ",
              shiny::code("selected_trace"),
              "."
            )
          )
        )
      )
    )
  )
}


# ============================================================
# SERVER
# ============================================================

mod_method_comparison_server <- function(id,
                                         fft_results = NULL,
                                         wavelet_results = NULL) {
  
  shiny::moduleServer(id, function(input, output, session) {
    
    # ------------------------------------------------------------
    # Helper: evaluate object only if it is reactive
    # ------------------------------------------------------------
    resolve_object <- function(x) {
      if (is.null(x)) return(NULL)
      
      if (shiny::is.reactive(x)) {
        return(x())
      }
      
      if (is.function(x)) {
        out <- tryCatch(x(), error = function(e) NULL)
        return(out)
      }
      
      x
    }
    
    # ------------------------------------------------------------
    # Helper: extract a data frame of metrics from a possible result list
    # ------------------------------------------------------------
    extract_metrics_df <- function(x) {
      x <- resolve_object(x)
      if (is.null(x)) return(NULL)
      
      if (is.data.frame(x)) return(x)
      
      if (is.list(x)) {
        candidate_names <- c(
          "summary_metrics",
          "summary_table",
          "metrics",
          "metric_table",
          "trace_metrics",
          "baseline_metrics"
        )
        
        for (nm in candidate_names) {
          if (!is.null(x[[nm]]) && is.data.frame(x[[nm]])) {
            return(x[[nm]])
          }
        }
        
        # Fallback: first data frame inside the list
        df_positions <- which(vapply(x, is.data.frame, logical(1)))
        if (length(df_positions) > 0) {
          return(x[[df_positions[1]]])
        }
      }
      
      NULL
    }
    
    # ------------------------------------------------------------
    # Helper: convert metrics to long format
    # ------------------------------------------------------------
    standardize_metrics <- function(df, method_name) {
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
        return(data.frame(
          Method = character(0),
          Metric = character(0),
          Value = numeric(0),
          Family = character(0),
          stringsAsFactors = FALSE
        ))
      }
      
      names(df) <- make.names(names(df))
      
      # Case 1: already long format
      metric_col <- names(df)[tolower(names(df)) %in% c("metric", "metrics", "variable", "name")]
      value_col <- names(df)[tolower(names(df)) %in% c("value", "values", "estimate")]
      
      if (length(metric_col) > 0 && length(value_col) > 0) {
        out <- data.frame(
          Method = method_name,
          Metric = as.character(df[[metric_col[1]]]),
          Value = suppressWarnings(as.numeric(df[[value_col[1]]])),
          stringsAsFactors = FALSE
        )
      } else {
        # Case 2: wide format, use numeric columns
        num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
        
        if (length(num_cols) == 0) {
          return(data.frame(
            Method = character(0),
            Metric = character(0),
            Value = numeric(0),
            Family = character(0),
            stringsAsFactors = FALSE
          ))
        }
        
        # If there are multiple rows, use the first row as summary
        first_row <- df[1, num_cols, drop = FALSE]
        
        out <- data.frame(
          Method = method_name,
          Metric = num_cols,
          Value = suppressWarnings(as.numeric(first_row[1, ])),
          stringsAsFactors = FALSE
        )
      }
      
      out <- out[is.finite(out$Value), , drop = FALSE]
      out$Metric <- clean_metric_names(out$Metric)
      out$Family <- classify_metric_family(out$Metric)
      out
    }
    
    # ------------------------------------------------------------
    # Helper: clean metric names for display
    # ------------------------------------------------------------
    clean_metric_names <- function(x) {
      x <- as.character(x)
      x <- gsub("\\.", "_", x)
      x <- gsub("DeltaF_over_F0", "DeltaF_over_Baseline", x)
      x <- gsub("F0", "Baseline", x)
      x
    }
    
    # ------------------------------------------------------------
    # Helper: classify metric family
    # ------------------------------------------------------------
    classify_metric_family <- function(metric) {
      metric <- clean_metric_names(metric)
      
      out <- rep("Other metrics", length(metric))
      
      out[grepl("Number|Peak|Event|Frequency", metric, ignore.case = TRUE)] <- "Event detection"
      out[grepl("Time|Onset|Rise|FWHM|FWHP|Duration", metric, ignore.case = TRUE)] <- "Time-domain metrics"
      out[grepl("Amplitude|Prominence|DeltaF|Signal|Baseline", metric, ignore.case = TRUE)] <- "Signal amplitude metrics"
      out[grepl("AUC|Area", metric, ignore.case = TRUE)] <- "Area under the curve"
      out[grepl("CWT|Ridge|Scale|Wavelet", metric, ignore.case = TRUE)] <- "Wavelet ridge metrics"
      
      out
    }
    
    # ------------------------------------------------------------
    # Static conceptual comparison table
    # ------------------------------------------------------------
    output$conceptual_comparison_table <- DT::renderDT({
      df <- data.frame(
        Component = c(
          "Primary goal",
          "Signal preparation",
          "Event detection principle",
          "Main strength",
          "Main limitation",
          "Baseline dependence",
          "Best use case",
          "Main outputs"
        ),
        `FFT + Baseline Analysis` = c(
          "Smooth the trace and detect peaks using user-defined peak parameters.",
          "Optional polynomial detrending followed by FFT low-pass smoothing.",
          "Events are detected as local peaks in the smoothed calcium trace.",
          "Direct, intuitive, and highly interpretable peak-based workflow.",
          "Can be sensitive to peak-height, peak-distance, ascent, descent, and baseline choices.",
          "Baseline affects amplitude, AUC, DeltaF/Baseline, FWHM, FWHP, and accepted-event metrics.",
          "Clear transient signals with well-defined peaks.",
          "Peak table, transient metrics, baseline sensitivity, FFT selection criteria, and summary plots."
        ),
        `Wavelet Ridgewalking` = c(
          "Detect events as multiscale time-scale structures.",
          "Uses a continuous wavelet transform approximation across temporal scales.",
          "Events are detected from connected local maxima across wavelet scales.",
          "Can capture broader or low-amplitude events that may not appear as sharp peaks.",
          "More abstract and depends on scale range, ridge length, and CWT amplitude quantile.",
          "Detection is less directly baseline-driven, but post-detection metrics still depend on baseline.",
          "Complex, broad, noisy, or multiscale calcium transients.",
          "Ridge event table, CWT map, event classification, AUC, transient metrics, and baseline sensitivity."
        ),
        check.names = FALSE
      )
      
      DT::datatable(
        df,
        rownames = FALSE,
        options = list(
          pageLength = 8,
          dom = "t",
          scrollX = TRUE
        )
      )
    })
    
    # ------------------------------------------------------------
    # Interpretation table
    # ------------------------------------------------------------
    output$interpretation_table <- DT::renderDT({
      df <- data.frame(
        Situation = c(
          "Both methods detect similar event number and timing",
          "FFT detects more events than Wavelet",
          "Wavelet detects more events than FFT",
          "FFT detects sharper events, Wavelet detects broader events",
          "Amplitude and AUC change strongly across baseline definitions",
          "Timing is stable but amplitude is not",
          "Both methods disagree strongly"
        ),
        Interpretation = c(
          "The calcium event pattern is likely robust to the event-detection pipeline.",
          "The FFT workflow may be splitting complex transients, detecting noise, or being influenced by permissive peak parameters.",
          "The wavelet workflow may be capturing broad, low-amplitude, or multiscale events not detected by fixed peak rules.",
          "The signal may contain events with different temporal widths; both methods may be highlighting complementary biological features.",
          "Quantitative conclusions based on amplitude, AUC, or FWHM should be reported as baseline-sensitive.",
          "Event occurrence may be reliable, but quantitative intensity-based interpretation should be cautious.",
          "The result should be interpreted as method-sensitive; additional parameter exploration or manual inspection is recommended."
        ),
        Recommended_Action = c(
          "Report agreement and use both methods as supporting evidence.",
          "Review peak height, minimum distance, ascent, descent, and baseline settings.",
          "Inspect the CWT map and verify whether additional wavelet events are biologically plausible.",
          "Report both sharp and broad event behavior separately if biologically relevant.",
          "Include baseline sensitivity analysis in the report.",
          "Emphasize timing/frequency metrics more than amplitude/AUC metrics.",
          "Do not overstate the biological conclusion without sensitivity analysis."
        ),
        check.names = FALSE
      )
      
      DT::datatable(
        df,
        rownames = FALSE,
        options = list(
          pageLength = 7,
          dom = "t",
          scrollX = TRUE
        )
      )
    })
    
    # ------------------------------------------------------------
    # Recommended metrics table
    # ------------------------------------------------------------
    output$recommended_metrics_table <- DT::renderDT({
      df <- data.frame(
        Metric_Family = c(
          "Event detection",
          "Event detection",
          "Time-domain metrics",
          "Time-domain metrics",
          "Time-domain metrics",
          "Signal amplitude metrics",
          "Signal amplitude metrics",
          "Area under the curve",
          "Wavelet ridge metrics",
          "Baseline sensitivity"
        ),
        Metric = c(
          "Number of detected events",
          "Frequency",
          "Mean peak occurrence time",
          "Mean onset / rise time",
          "Mean FWHM / FWHP",
          "Mean amplitude above baseline",
          "Mean prominence",
          "Mean event AUC / total AUC",
          "Mean CWT amplitude / ridge length / scale",
          "Coefficient of variation across baseline methods"
        ),
        Why_it_matters = c(
          "Shows whether methods detect similar event counts.",
          "Summarizes event occurrence over time.",
          "Evaluates whether methods identify events at similar temporal locations.",
          "Describes the temporal dynamics of event initiation.",
          "Captures transient width and duration.",
          "Quantifies signal intensity relative to the selected baseline.",
          "Measures how strongly events stand above neighboring valleys.",
          "Summarizes total event burden or integrated calcium activity.",
          "Provides wavelet-specific information about multiscale event structure.",
          "Shows whether conclusions depend strongly on baseline definition."
        ),
        Recommended_Interpretation = c(
          "Compare directly across FFT and Wavelet.",
          "Compare directly, but verify duration/time scale.",
          "Compare only if both methods detect similar event groups.",
          "Interpret alongside event plots.",
          "Do not mix with amplitude metrics in the same scale.",
          "Interpret with baseline sensitivity results.",
          "Useful for comparing local event contrast.",
          "Highly baseline-sensitive; report with caution.",
          "Wavelet-specific; use mainly to characterize wavelet-detected events.",
          "Use as a robustness diagnostic."
        ),
        check.names = FALSE
      )
      
      DT::datatable(
        df,
        rownames = FALSE,
        options = list(
          pageLength = 10,
          dom = "tip",
          scrollX = TRUE
        )
      )
    })
    
    # ------------------------------------------------------------
    # Direct comparison reactive data
    # ------------------------------------------------------------
    comparison_long <- shiny::reactive({
      fft_df <- extract_metrics_df(fft_results)
      wavelet_df <- extract_metrics_df(wavelet_results)
      
      fft_long <- standardize_metrics(fft_df, "FFT + Baseline Analysis")
      wavelet_long <- standardize_metrics(wavelet_df, "Wavelet Ridgewalking")
      
      dplyr::bind_rows(fft_long, wavelet_long)
    })
    
    # ------------------------------------------------------------
    # Direct comparison status
    # ------------------------------------------------------------
    output$comparison_status <- shiny::renderUI({
      dat <- comparison_long()
      
      if (is.null(fft_results) || is.null(wavelet_results) || nrow(dat) == 0) {
        return(
          shiny::div(
            class = "ci-status-box",
            shiny::strong("Current status: "),
            "This panel is ready for direct comparison, but the FFT and Wavelet modules are not yet returning reactive result objects to this module. ",
            "For now, use this tab as a placeholder for the next implementation step. The interpretation guide and metric tables above are fully usable."
          )
        )
      }
      
      shiny::div(
        class = "ci-status-box",
        shiny::strong("Connected mode: "),
        "Reactive outputs were detected. The table and plot below summarize comparable numeric metrics returned by the FFT and Wavelet modules."
      )
    })
    
    # ------------------------------------------------------------
    # Direct comparison table
    # ------------------------------------------------------------
    output$direct_comparison_table <- DT::renderDT({
      dat <- comparison_long()
      
      if (nrow(dat) == 0) {
        dat <- data.frame(
          Method = c("FFT + Baseline Analysis", "Wavelet Ridgewalking"),
          Metric = c("Not connected", "Not connected"),
          Value = c(NA_real_, NA_real_),
          Family = c("Implementation status", "Implementation status"),
          Note = c(
            "Return a reactive summary_metrics table from the FFT module.",
            "Return a reactive summary_metrics table from the Wavelet module."
          ),
          check.names = FALSE
        )
      }
      
      DT::datatable(
        dat,
        rownames = FALSE,
        options = list(
          pageLength = 12,
          scrollX = TRUE
        )
      )
    })
    
    # ------------------------------------------------------------
    # Direct comparison plot
    # ------------------------------------------------------------
    output$direct_comparison_plot <- shiny::renderPlot({
      dat <- comparison_long()
      
      if (nrow(dat) == 0) {
        graphics::plot.new()
        graphics::text(
          x = 0.5,
          y = 0.58,
          labels = "Direct method comparison is ready,\nbut FFT and Wavelet reactive outputs are not connected yet.",
          cex = 1.1
        )
        graphics::text(
          x = 0.5,
          y = 0.42,
          labels = "Use this module as the interpretation and reporting panel\nuntil both analysis modules return summary metrics.",
          cex = 0.9,
          col = "gray40"
        )
        return(invisible(NULL))
      }
      
      dat <- dat[is.finite(dat$Value), , drop = FALSE]
      
      if (nrow(dat) == 0) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "No finite numeric metrics available for plotting.")
        return(invisible(NULL))
      }
      
      # Keep a manageable number of metrics for plotting
      metric_counts <- sort(table(dat$Metric), decreasing = TRUE)
      common_metrics <- names(metric_counts[metric_counts >= 1])
      common_metrics <- head(common_metrics, 18)
      plot_dat <- dat[dat$Metric %in% common_metrics, , drop = FALSE]
      
      ggplot2::ggplot(
        plot_dat,
        ggplot2::aes(
          x = Metric,
          y = Value,
          fill = Method
        )
      ) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::facet_wrap(~ Family, scales = "free_y") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "Direct comparison of returned summary metrics",
          subtitle = "Metrics are grouped by family because they have different units and scales",
          x = "Metric",
          y = "Value"
        ) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          legend.position = "bottom",
          strip.background = ggplot2::element_rect(fill = "gray95"),
          strip.text = ggplot2::element_text(face = "bold"),
          plot.title = ggplot2::element_text(face = "bold")
        )
    })
  })
}