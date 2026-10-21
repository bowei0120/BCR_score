#!/usr/bin/env Rscript
# ROC analysis for BCR scores.
#
# Usage:
#   Rscript ROC.R --table BCR_score_result/score_PRJNA797175.csv \
#       --case UC --control Control --score score_UC \
#       [--label label] [--cutoff 3.258] [--out roc.pdf]
#
# Arguments:
#   --table    Score table (CSV/TSV) with a label column and score column(s).
#              CSV is read with row.names = 1 (sample id as row names).
#   --case     Case label in the label column (e.g. UC, CD)
#   --control  Control label in the label column (e.g. Control, HT)
#   --score    Score column name (e.g. score_UC, score_CD)
#   --label    Label column name (default: label)
#   --cutoff   Optional fixed threshold to report sensitivity/specificity at
#   --out      Optional output plot file; extension picks the device
#              (.pdf / .png / .jpg / .svg). Omit to plot interactively.
#
# Prints AUC (95% CI), best Youden cutoff, and metrics at --cutoff (if given).

suppressPackageStartupMessages(library(pROC))

parse_args <- function(argv) {
  opts <- list(label = "label", cutoff = NA_real_, out = NULL)
  i <- 1
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[i])
    if (!key %in% c("table", "case", "control", "score", "label", "cutoff", "out")) {
      stop("Unknown argument: ", argv[i])
    }
    if (i + 1 > length(argv)) stop("Missing value for --", key)
    val <- argv[i + 1]
    if (key == "cutoff") val <- as.numeric(val)
    opts[[key]] <- val
    i <- i + 2
  }
  for (req in c("table", "case", "control", "score")) {
    if (is.null(opts[[req]])) stop("Missing required argument: --", req)
  }
  opts
}

calculate_roc <- function(exp_m, class_m, case_name, control_name, score_name) {
  select_column <- class_m %in% c(case_name, control_name)

  exp_m <- exp_m[select_column, score_name]
  class_m <- class_m[select_column]
  class_m <- class_m %in% case_name

  roc_curve <- roc(class_m, exp_m, quiet = TRUE)
  return(roc_curve)
}

plot_roc_ci <- function(roc_curve) {
  # AUC 95% CI
  ci_auc <- ci.auc(roc_curve)

  # ROC curve 95% CI
  ci_roc <- ci.se(
    roc_curve,
    specificities = seq(0, 1, 0.01),
    conf.level = 0.95
  )

  # ROC + CI
  plot(
    roc_curve,
    col = "#0072B5FF",
    lwd = 2,
    print.auc = FALSE
  )

  plot(
    ci_roc,
    type = "shape",
    col = "#0072B520",
    border = NA
  )

  # ROC curves
  lines(
    roc_curve$specificities,
    roc_curve$sensitivities,
    col = "#0072B5FF",
    lwd = 2
  )

  # AUC + 95% CI
  text(
    0.35, 0.10,
    labels = sprintf("AUC = %.3f\n(95%% CI %.3f-%.3f)",
                     as.numeric(auc(roc_curve)), ci_auc[1], ci_auc[3]),
    cex = 2,
    col = "#0072B5FF"
  )
}

get_best_cutoff <- function(roc_curve) {
  best <- coords(
    roc_curve,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  data.frame(
    cutoff = best$threshold,
    sensitivity = best$sensitivity,
    specificity = best$specificity
  )
}

get_metrics_at_cutoff <- function(roc_curve, cutoff) {
  result <- coords(
    roc_curve,
    x = cutoff,
    input = "threshold",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  data.frame(
    cutoff = result$threshold,
    sensitivity = result$sensitivity,
    specificity = result$specificity
  )
}

open_device <- function(out) {
  ext <- tolower(tools::file_ext(out))
  switch(ext,
    pdf = pdf(out),
    png = png(out, width = 1600, height = 1600, res = 300),
    jpg = , jpeg = jpeg(out, width = 1600, height = 1600, res = 300),
    svg = svg(out),
    stop("Unsupported output extension: .", ext, " (use pdf/png/jpg/svg)")
  )
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))

  score_df <- if (grepl("\\.csv$", opts$table, ignore.case = TRUE)) {
    read.csv(opts$table, row.names = 1)
  } else {
    read.table(opts$table, header = TRUE, sep = "\t", row.names = 1)
  }

  for (col in c(opts$label, opts$score)) {
    if (!col %in% colnames(score_df)) {
      stop("Column '", col, "' not found in ", opts$table,
           " (available: ", paste(colnames(score_df), collapse = ", "), ")")
    }
  }

  roc_curve <- calculate_roc(score_df, score_df[[opts$label]],
                             opts$case, opts$control, opts$score)

  cat(sprintf("%s vs %s (%s):\n", opts$case, opts$control, opts$score))
  cat(sprintf("  AUC = %.3f (95%% CI %.3f-%.3f)\n",
              as.numeric(auc(roc_curve)),
              as.numeric(ci.auc(roc_curve)[1]),
              as.numeric(ci.auc(roc_curve)[3])))

  best <- get_best_cutoff(roc_curve)
  cat(sprintf("  best cutoff (Youden) = %.3f  sens = %.3f  spec = %.3f\n",
              best$cutoff, best$sensitivity, best$specificity))

  if (!is.na(opts$cutoff)) {
    at_cut <- get_metrics_at_cutoff(roc_curve, opts$cutoff)
    cat(sprintf("  at cutoff %.3f:            sens = %.3f  spec = %.3f\n",
                at_cut$cutoff, at_cut$sensitivity, at_cut$specificity))
  }

  if (!is.null(opts$out)) {
    open_device(opts$out)
    plot_roc_ci(roc_curve)
    dev.off()
    cat("  plot written to ", opts$out, "\n", sep = "")
  } else if (interactive()) {
    plot_roc_ci(roc_curve)
  }
}

main()
