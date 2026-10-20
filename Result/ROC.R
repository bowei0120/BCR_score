calculate_roc  <- function(exp_m, class_m, case_name, control_name, score_name){
  select_column = class_m %in% c(case_name, control_name)

  exp_m = exp_m[select_column, score_name]
  class_m = class_m[select_column]
  class_m = class_m %in% case_name
  
  roc_curve <- roc(class_m, exp_m)
  return(roc_curve)
}

plot_roc_ci <- function(roc_curve){
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
    labels = sprintf(
      "AUC = %.3f\n(95%% CI %.3f–%.3f)",
      as.numeric(auc(roc_curve)),
      ci_auc[1],
      ci_auc[3]
    ),
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

### Usage ###
# Read file
score_df = read.csv("score_PRJNA797175.csv", row.names = 1)

# Generate ROC curve
roc_curve = calculate_roc(score_df, score_df$label, "Control", c("UC"), "score_UC")

# Get best cutoff
get_best_cutoff(roc_curve)

# Get metrics at a speciallized cutoff
get_metrics_at_cutoff(roc_curve, 3.258)

# Generate ROC plot
plot_roc_ci(roc_curve)