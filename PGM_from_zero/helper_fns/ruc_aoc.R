library(pROC)

auc_roc_analysis <- function(true_edge, predicted_edge_prob, plot = FALSE) {
  roc_curve <- roc(as.numeric(true_edge), as.numeric(predicted_edge_prob), quiet = TRUE)
  if (plot) {
    with(
      roc_curve, 
      plot(
        1-specificities, 
        sensitivities, 
        type = 'l',
        main = "ROC for QPGM",
        xlab = "False Positive Rate",
        ylab = "True Positive Rate",
        xlim = c(0,1), ylim = c(0,1),
        lwd = 2
      )
    )
    axis(1, axTicks(1), labels=F)
    lines(c(0, 1), c(0, 1), col="grey")
  }
  return(roc_curve$auc)
}
