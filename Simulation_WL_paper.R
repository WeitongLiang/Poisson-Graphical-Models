# ====== Environment set-up =======
library(XMRF)
library(pROC)
library(verification)
library(glmnet)       # For penalized GLMs
library(igraph)       # Network generation and analysis
library(ggplot2)
library(gridExtra)
library(grid)

## ===== Set seed =====
set.seed(1)

## ===== Set work directory =====
setwd("~/Documents/Duke University/Courses/STA 841 Categotrical Data Analysis/Poisson-Graphical-Models")
purrr::walk(list.files("./PGM_from_zero/helper_fns/", pattern = "*.R$", full.names=TRUE), source,.GlobalEnv)

# ===== Simulation =====
results_list <- list()  # storage model performance result
N_range <- c(5*(4:25), 200*(1:4), 1000*(1:10))  # range of number of samples
Success_range <- numeric(length = length(N_range))  # stores the success rate
p_values <- c(4, 9, 36)
M <- 512 # number of runings per combination of (p, n)
omega <- -0.1         # Edge weight (negative for Poisson PGM)

for (p in p_values) {
  message(paste("\nRunning experiments for p =", p))
  # Generate 4-nearest neighbor lattice graph
  adj_matrix <- make_lattice(length = sqrt(p), dim = 2, nei = 1) |> 
    as_adjacency_matrix(sparse = FALSE)
  diag(adj_matrix) <- 0 # Remove self-loops
  
  # Simulate data
  X_all <- simulate_pgm_data(adj_matrix, n = 2*max(N_range), omega)
  
  for (N in 1:length(N_range)) {
    n <- N_range[N] # number of samples
    tmp_success <- numeric(length = M)
    for (m in 1:M){
      # show current progress
      cat(paste("Calculating success rate in progress (n = ", n, "): ", floor(100*(m/M)), "%", collapse=""),"\r")
      flush.console()
      
      # data for the m-th simulation
      X <- X_all[sample(1:nrow(X_all), n, replace = FALSE), ]
      # fit edges based on Original Poisson Graphical Model
      THETA <- fit_Original_Poisson_Graphical_Model(X = X, lambda = 2 * sqrt(log(p) / n))
      
      # store results
      tmp_success[m] <- calculate_success_rate(adj_matrix, THETA)
    }
    Success_range[N] <- mean(tmp_success)
  }
  
  # Store results for this p
  results_list[[as.character(p)]] <- list(
    p = p,
    Success_range = Success_range
  )
}

# ===== Draw Figures =====
result_df <- data.frame(
  n = N_range,
  "4" = results_list[['4']]$Success_range,
  "9" = results_list[['9']]$Success_range,
  "36" = results_list[['36']]$Success_range
) %>% 
  setNames(c("n",'p_4','p_9','p_36'))

ggplot(result_df, aes(x = n)) +
  geom_line(aes(y = p_4, color = "4"), linewidth = 1.2) +
  geom_line(aes(y = p_9, color = "9"), linewidth = 1.2) +
  geom_line(aes(y = p_36, color = "36"), linewidth = 1.2) +
  scale_color_manual(
    name = "p",  # Legend title
    values = c("4" = "#E69F00", "9" = "#56B4E9", "36" = "#009E73"),  # Colors
    breaks = c("4", "9", "36")  # Ensures legend order: 4, 9, 36
  ) +
  labs(
    title = "Relationship Between n and Success Rate",
    subtitle = "Comparing different values of dimension p",
    x = "Sample Size (n)",
    y = "Success Rate",
    caption = "Simplified simulation replicating Yang, Eunho, et al. (2012)."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "gray50"),
    axis.title = element_text(face = "bold"),
    plot.caption = element_text(color = "gray50", hjust = 1)
  )

# saveRDS(result_df, file = "./PGM_from_zero/RESULT_SuccessRate_vs_SampleSize/SuccessRate_vs_n.RData")

## ===== Refit model under p=36, n = 10000 =====
n = 10000
# OR set n = 3000
X <- X_all[sample(1:nrow(X_all), n, replace = FALSE), ]
THETA <- fit_Original_Poisson_Graphical_Model(X = X, lambda = 2 * sqrt(log(p) / n))

create_matrix_plot <- function(matrix, title) {
  df <- reshape2::melt(matrix)
  df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))  # Reverse y-axis levels
  
  ggplot(df, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "black", linewidth = 0.2) +
    scale_fill_gradient(low = "white", high = "black", guide = "none") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(title = title, x = "", y = "") +
    coord_fixed() +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA),
          axis.text = element_blank(),
          panel.grid = element_blank())
}

# Generate plots
p1 <- create_matrix_plot(THETA, "Fitted Adjacency Matrix by OPGM")
p2 <- create_matrix_plot(adj_matrix, "True Adjacency Matrix")

# Arrange with global caption
combined <- grid.arrange(
  p1, p2, 
  ncol = 2,
  bottom = textGrob("p = 36, n = 10000", 
                    gp = gpar(fontface = "plain", fontsize = 15),
                    hjust = 0.4, x = 0.5)
)
