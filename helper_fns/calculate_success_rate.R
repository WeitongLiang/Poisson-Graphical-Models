calculate_success_rate <- function(true_adj, pred_adj) {
  # Only compare upper triangular parts to avoid double-counting undirected edges
  upper_true <- true_adj[upper.tri(true_adj)]
  upper_pred <- pred_adj[upper.tri(pred_adj)]
  
  # Calculate proportion of matching edges
  sum(upper_true == upper_pred) / length(upper_true)
}