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
set.seed(42)

## ===== Set work directory =====
setwd("~/Documents/Duke University/Courses/STA 841 Categotrical Data Analysis/Poisson-Graphical-Models")
purrr::walk(list.files("./PGM_from_zero/helper_fns/", pattern = "*.R$", full.names=TRUE), source,.GlobalEnv)

data("brcadat")
brca <- t(processSeq(t(brcadat), PercentGenes=1)) # n = 445, p = 353

n <- ncol(brca)
p <- nrow(brca)

X <- t(as.matrix(brca))

# ===== Fit QPGM and SPGM with optim() =====
THETA_qpgm <- fit_Quadratic_Poisson_Graphical_Model_optim(X = X, lambda = 1e3)
THETA_spgm <- fit_Sublinear_Poisson_Graphical_Model_optim(X = X, lambda = 1e3)
























