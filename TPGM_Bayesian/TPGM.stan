functions {
  real truncated_poisson_lpmf(int x, real lambda) {
    return poisson_lpmf(x | lambda) - poisson_lccdf(0 | lambda);
  }
}

data {
  int<lower=0> n;  // 样本数
  int<lower=1> N;  // 每个样本的节点数
  int<lower=1> X[n, N];  // 观测数据（X >= 1）
}

parameters {
  matrix[N, N] J_raw;  // 非对称的原始参数
  vector[N] h;         // 偏置项
}

transformed parameters {
  matrix[N, N] J_sym;  // 对称的相互作用矩阵
  for (i in 1:N) {
    for (j in 1:N) {
      if (i < j) {
        J_sym[i, j] = J_raw[i, j];
        J_sym[j, i] = J_sym[i, j];  // 强制对称
      } else if (i == j) {
        J_sym[i, j] = 0;  // 对角线为0
      }
    }
  }
}

model {
  // 先验
  h ~ normal(0, 1);
  for (i in 1:N) {
    for (j in (i+1):N) {
      J_raw[i, j] ~ normal(0, 0.5);  // 仅上三角部分有先验
    }
  }

  // 伪似然
  for (k in 1:n) {
    for (i in 1:N) {
      real lambda = exp(h[i] + dot_product(J_sym[i, ], to_vector(X[k, ])) - J_sym[i, i] * X[k, i]);
      X[k, i] ~ truncated_poisson(lambda);
    }
  }
}