library(rstan)
library(extraDistr)

# 参数设置
N <- 3       # 每个样本的节点数
n <- 100     # 样本数
J <- matrix(c(0, 0.2, -0.1, 
             0.2, 0, 0.3, 
             -0.1, 0.3, 0), nrow = N)  # 真实相互作用矩阵（对称）
h <- c(0.5, -0.2, 0.1)                  # 真实偏置项

# 模拟数据（Gibbs采样）
X <- matrix(1, n, N)  # 初始化所有X >= 1
for (k in 1:n) {
  for (i in 1:N) {
    lambda <- exp(h[i] + J[i, -i] %*% X[k, -i])
    X[k, i] <- rtpois(1, lambda = lambda, a = 1, b = Inf)
  }
}

stan_data <- list(n = n, N = N, X = X)
fit <- stan(
  file = "/Users/yanyan/Documents/Duke University/Courses/STA 841 Categotrical Data Analysis/Poisson-Graphical-Models/Truncated_Poisson_Graphical_Model/TPGM.stan",
  data = stan_data,
  chains = 4,
  iter = 2000
)

# 输出结果
print(fit, pars = c("J_sym", "h"))

# 提取后验均值
posterior <- extract(fit)
J_est <- apply(posterior$J_sym, c(2, 3), mean)  # 估计的相互作用矩阵
h_est <- colMeans(posterior$h)              # 估计的偏置项

cat("真实相互作用矩阵 J:\n"); print(J)
cat("估计相互作用矩阵 J_est:\n"); print(J_est)

cat("True h:\n", h_true, "\n")
cat("True J:\n"); print(J_true)