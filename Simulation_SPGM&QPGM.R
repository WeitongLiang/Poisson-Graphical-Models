library(igraph)
library(pROC)
library(PRROC)
library(Matrix)
library(reshape2)
library(ggplot2)

set.seed(42)

## ===== Set work directory =====
setwd("C:/Users/30604/Desktop/Duke/STAT_841/Poisson-Graphical-Models")
purrr::walk(list.files("./helper_fns/", pattern = "*.R$", full.names=TRUE), source,.GlobalEnv)
library(igraph)
library(Matrix)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(pROC)

simulate_pgm_data <- function(adj_matrix, n, omega, burn_in = 1e4, rate_cap = 50) {
  p <- nrow(adj_matrix)
  X <- matrix(0, nrow = n, ncol = p)
  theta <- adj_matrix * omega
  x_current <- rpois(p, lambda = 1)
  for (iter in 1:(burn_in + n)) {
    for (s in 1:p) {
      neighbors <- which(theta[s, ] != 0)
      eta <- sum(theta[s, neighbors] * x_current[neighbors])
      lambda <- min(exp(eta), rate_cap)
      x_current[s] <- rpois(1, lambda = lambda)
    }
    if (iter > burn_in) X[iter - burn_in, ] <- x_current
  }
  return(X)
}

plot_theta_matrix <- function(mat, title) {
  df <- melt(abs(mat))
  ggplot(df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "white", high = "black") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 14) +
    theme(axis.text = element_text(size = 10),
          axis.ticks = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
}

# Sweep settings
p_list <- c(4, 9, 16)
n_list <- c(10000)
lambda_list <- c(0.03)
omega <- 0.1

# Output storage
plot_list <- list()
metric_df <- data.frame()
plot_i <- 1
set.seed(2025)

for (p in p_list) {
  for (n in n_list) {
    for (lambda_val in lambda_list) {
      cat(sprintf("Running p=%d, n=%d, λ=%.3f\n", p, n, lambda_val))
      
      # generate graph
      g <- sample_pa(p, power = 1.2, directed = FALSE)
      A_true <- as.matrix(as_adjacency_matrix(g))
      theta_true <- A_true * omega
      
      # simulate data
      X_sim <- simulate_pgm_data(A_true, n = p*n, omega = omega, rate_cap = 50)
      
      # fit model
      THETA_spgm <- fit_Quadratic_Poisson_Graphical_Model_optim(X_sim, lambda = lambda_val)
      theta_est <- (THETA_spgm + t(THETA_spgm)) / 2
      
      # upper triangle comparison
      true_edges <- as.vector(A_true[upper.tri(A_true)])
      score_edges <- as.vector(theta_est[upper.tri(theta_est)])
      bin_est_edges <- (score_edges > 1e-3) * 1
      
      tp <- sum(true_edges == 1 & bin_est_edges == 1)
      fp <- sum(true_edges == 0 & bin_est_edges == 1)
      fn <- sum(true_edges == 1 & bin_est_edges == 0)
      
      precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
      recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
      f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
      
      auc_val <- tryCatch({
        roc_obj <- pROC::roc(response = true_edges, predictor = score_edges, quiet = TRUE)
        pROC::auc(roc_obj)
      }, error = function(e) NA)
      
      # Titles and plots
      title_true <- paste0("True θ\n(p=", p, ", n=", n, ")")
      title_est <- sprintf("λ=%.3f\nTP=%d  FP=%d\nF1=%.2f  AUC=%.2f", 
                           lambda_val, tp, fp, f1, auc_val)
      
      g1 <- plot_theta_matrix(theta_true, title_true)
      g2 <- plot_theta_matrix(theta_est, title_est)
      plot_list[[plot_i]] <- g1
      plot_list[[plot_i + 1]] <- g2
      plot_i <- plot_i + 2
      
      # Store metrics
      metric_df <- rbind(metric_df, data.frame(
        p = p,
        n = n,
        lambda = lambda_val,
        TP = tp,
        FP = fp,
        FN = fn,
        Precision = round(precision, 3),
        Recall = round(recall, 3),
        F1 = round(f1, 3),
        AUC = round(as.numeric(auc_val), 3)
      ))
    }
  }
}

# Save visual output
grDevices::pdf("spgm_visual_grid.pdf", width = 12, height = 14)
grid.arrange(grobs = plot_list, ncol = 3)
dev.off()

# Save metrics as CSV
write.csv(metric_df, "qpgm_metrics.csv", row.names = FALSE)
cat("\n✅ All done! Saved PDF and CSV.\n")
